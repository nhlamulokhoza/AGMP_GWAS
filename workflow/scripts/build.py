"""Build the curated AGMP GWAS table for one GWAS Catalog release.

Streams the association zip and ancestry TSV from the EBI FTP into memory
(with size and CRC checks), selects studies relevant to African and
African-diaspora populations, normalises and rescues variant IDs, explodes,
maps regions, annotates genes, and writes ONE compressed parquet table plus a
small list of variants that could not be mapped. No raw downloads or
intermediate tables are written to disk.

Study inclusion routes (each ancestry row gets the first that applies):
  explicit     broad ancestral category is African / Greater Middle Eastern
  label        broad ancestral category is NR / Other
  country      an African country of recruitment or origin, any ancestry label
  description  the row's own sample description names African/diaspora samples
Only "explicit" rows count as resolved; the rest go to NLP triage.

Output columns are those of agmp_gwas_pipeline_June2026.ipynb; the extra
internal columns are dropped at export.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
import zipfile
from datetime import date
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from annotate import annotate_genes  # noqa: E402
from common import (fetch_bytes, get_logger, http_session, load_yaml,  # noqa: E402
                    read_previous_ids, release_info)
from variant_rescue import Rescuer, parse_positional  # noqa: E402
from vep import consequences  # noqa: E402

log = get_logger("build")

ASSOC_COLUMNS = {
    "SNPS": "id",
    "PUBMEDID": "PUBMEDID",
    "DISEASE/TRAIT": "phenotype",
    "P-VALUE": "p-value",
    "REPORTED GENE(S)": "gene_symbol_mined",
    "MAPPED_GENE": "curated_gene_symbol",
    "CHR_ID": "chromosome",
    "STUDY": "title",
    "DATE": "publication_year",
    "INTERGENIC": "intergenic",
    "STUDY ACCESSION": "id_in_source",
    "SNP_GENE_IDS": "gene_name",
    "INITIAL SAMPLE SIZE": "initial_sample_size",
    "REPLICATION SAMPLE SIZE": "replication_sample_size",
    "CONTEXT": "context",
    # internal, used for variant ID normalisation and rescue
    "STRONGEST SNP-RISK ALLELE": "_risk_allele",
    "SNP_ID_CURRENT": "_snp_current",
    "MERGED": "_merged",
}

# GWAS Catalog headers contain typos; accept either spelling.
ANCESTRY_COLUMNS = {
    "STUDY ACCESSION": "id_in_source",
    "PUBMEDID": "ancestry_pubmedid", "PUBMED ID": "ancestry_pubmedid",
    "FIRST AUTHOR": "first_author",
    "DATE": "ancestry_date",
    "INITIAL SAMPLE DESCRIPTION": "initial_sample_desc",
    "REPLICATION SAMPLE DESCRIPTION": "replication_sample_desc",
    "STAGE": "study_stage",
    "NUMBER OF INDIVDUALS": "num_individuals", "NUMBER OF INDIVIDUALS": "num_individuals",
    "BROAD ANCESTRAL CATEGORY": "broad_ancestral_category",
    "COUNTRY OF ORIGIN": "country_of_origin",
    "COUNTRY OF RECRUITMENT": "origin_of_participants",
    "ADDITONAL ANCESTRY DESCRIPTION": "ethnicity",
    "ADDITIONAL ANCESTRY DESCRIPTION": "ethnicity",
}
REQUIRED_ANCESTRY = {"id_in_source", "broad_ancestral_category", "origin_of_participants",
                     "country_of_origin", "ethnicity", "initial_sample_desc",
                     "replication_sample_desc", "study_stage"}
# Internal columns that must not make otherwise-identical records distinct.
INTERNAL = ["_id_reported", "variant_id_source", "inclusion_route", "rescue_build"]

RS_RE = re.compile(r"^rs\d+$", re.I)
SPLIT_RE = re.compile(r"[;,]")
INTERACTION_RE = re.compile(r"\s[xX]\s")


# --------------------------------------------------------------------------
# Download with integrity checks
# --------------------------------------------------------------------------
def fetch_checked(url: str, session, is_zip: bool, attempts: int = 3):
    for attempt in range(1, attempts + 1):
        info = release_info(url, session)
        data = fetch_bytes(url, log, session)
        problem = None
        if info["bytes"] and len(data) != info["bytes"]:
            problem = f"size {len(data)} != Content-Length {info['bytes']}"
        elif is_zip:
            try:
                with zipfile.ZipFile(io.BytesIO(data)) as zf:
                    bad = zf.testzip()
                if bad:
                    problem = f"CRC error in {bad}"
            except zipfile.BadZipFile as e:
                problem = f"bad zip: {e}"
        if not problem:
            return data, info
        log.warning("download check failed (%s), attempt %d/%d", problem, attempt, attempts)
        time.sleep(10 * attempt)
    raise RuntimeError(f"could not download a valid copy of {url}")


# --------------------------------------------------------------------------
# Ancestry: selection routes
# --------------------------------------------------------------------------
def african_countries(regions_cfg: dict) -> set[str]:
    african = {r for r in regions_cfg["regions"] if "Africa" in r and "American" not in r}
    names = {str(c).strip().lower() for r in african for c in regions_cfg["regions"][r]}
    names |= {k.strip().lower() for k, v in (regions_cfg.get("aliases") or {}).items()
              if v.strip().lower() in names}
    return names


def load_ancestry(raw: bytes, cfg: dict, regions_cfg: dict) -> pl.DataFrame:
    df = pl.read_csv(io.BytesIO(raw), separator="\t", quote_char=None, infer_schema=False)
    log.info("ancestry file: %d rows, %d columns", df.height, df.width)
    keep = {c: ANCESTRY_COLUMNS[c] for c in df.columns if c in ANCESTRY_COLUMNS}
    missing = REQUIRED_ANCESTRY - set(keep.values())
    if missing:
        raise ValueError(f"ancestry file is missing expected columns: {sorted(missing)}; "
                         f"found {df.columns}. Update ANCESTRY_COLUMNS in build.py.")
    df = df.select(list(keep)).rename(keep).fill_null("")

    sel = cfg["selection"]
    explicit_re = re.compile(sel["explicit_ancestry_regex"], re.I)
    label_re = re.compile(sel["unresolved_ancestry_regex"])
    desc_re = re.compile(sel["sample_description_regex"])
    afr = african_countries(regions_cfg)

    def has_african_country(*fields):
        return any(tok.strip().lower() in afr for f in fields for tok in f.split(","))

    def route(bac, stage, init_desc, rep_desc, recruit, origin):
        if explicit_re.search(bac):
            return "explicit"
        if label_re.search(bac):
            return "label"
        if sel["include_african_recruitment"] and has_african_country(recruit, origin):
            return "country"
        own = rep_desc if stage.strip().lower() == "replication" else init_desc
        if sel["include_african_sample_description"] and desc_re.search(own):
            return "description"
        return ""

    cols = ["broad_ancestral_category", "study_stage", "initial_sample_desc",
            "replication_sample_desc", "origin_of_participants", "country_of_origin"]
    keys = df.select(cols).unique()
    keys = keys.with_columns(
        pl.Series("inclusion_route", [route(*r) for r in keys.iter_rows()], dtype=pl.String))
    df = df.join(keys, on=cols, how="left").filter(pl.col("inclusion_route") != "")
    # Sample descriptions describe the whole study, not one ancestry row, so the
    # description route only applies to studies with no labelled row at all.
    # Otherwise every "European" row of a multi-ancestry study would come in.
    labelled = pl.col("inclusion_route").is_in(["explicit", "label"]).any().over("id_in_source")
    df = df.filter((pl.col("inclusion_route") != "description") | ~labelled)
    df = df.with_columns(
        pl.when(pl.col("inclusion_route") == "explicit").then(pl.lit("explicit"))
          .otherwise(pl.lit("unresolved")).alias("ancestry_status"))
    counts = dict(df.group_by("inclusion_route").agg(pl.col("id_in_source").n_unique()).rows())
    log.info("ancestry rows kept: %d from %d studies; studies by route: %s",
             df.height, df["id_in_source"].n_unique(), counts)
    return df, counts


# --------------------------------------------------------------------------
# Associations
# --------------------------------------------------------------------------
def load_associations(raw: bytes, studies: set[str]):
    """Stream the zipped TSV row by row, keeping only the selected studies."""
    t0 = time.time()
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        member = next(n for n in zf.namelist() if n.endswith((".tsv", ".txt")))
        with zf.open(member) as fh:
            reader = csv.reader(io.TextIOWrapper(fh, encoding="utf-8", newline=""),
                                delimiter="\t", quoting=csv.QUOTE_NONE)
            header = next(reader)
            missing = [c for c in ASSOC_COLUMNS if c not in header]
            if missing:
                raise ValueError(f"association file is missing expected columns: {missing}")
            idx = [header.index(c) for c in ASSOC_COLUMNS]
            s_idx = header.index("STUDY ACCESSION")
            cols = [[] for _ in idx]
            total = short = 0
            for row in reader:
                total += 1
                if len(row) < len(header):
                    short += 1
                if len(row) <= s_idx or row[s_idx] not in studies:
                    continue
                for j, i in enumerate(idx):
                    cols[j].append(row[i] if i < len(row) else "")
    df = pl.DataFrame({name: col for name, col in zip(ASSOC_COLUMNS.values(), cols)},
                      schema={n: pl.String for n in ASSOC_COLUMNS.values()})
    log.info("associations: %d rows scanned, %d kept (%.0fs)", total, df.height,
             time.time() - t0)
    if short:
        log.warning("%d association rows had fewer fields than the header", short)
    return df, total


# --------------------------------------------------------------------------
# Variant IDs
# --------------------------------------------------------------------------
def split_tokens(snps: str, risk: str, current: str, merged: str):
    """Yield (token, normalised_id, allele, kind) for one SNPS field."""
    if INTERACTION_RE.search(snps):
        for t in INTERACTION_RE.split(snps):
            yield t.strip(), "", "", "interaction"
        return
    tokens = [t.strip() for t in SPLIT_RE.split(snps) if t.strip()]
    alleles = {}
    for item in SPLIT_RE.split(risk):
        if "-" in item:
            snp, a = item.strip().rsplit("-", 1)
            alleles[snp.strip().lower()] = a.strip().upper()
    for t in tokens:
        allele = alleles.get(t.lower(), "")
        if RS_RE.match(t):
            rid = t.lower()
            if merged == "1" and current and len(tokens) == 1 and rid != f"rs{current}":
                yield t, f"rs{current}", allele, "merged_updated"
            else:
                yield t, rid, allele, "rsid"
        elif "/" in t and all(RS_RE.match(p) for p in t.split("/")):
            yield t, f"rs{current}" if current else t.split("/")[0].lower(), allele, "synonym"
        elif parse_positional(t):
            yield t, "", allele, "positional"
        else:
            yield t, "", allele, "other"


def normalise_ids(df: pl.DataFrame, cfg: dict, base: Path, meta: dict):
    """Explode SNPS into one row per variant and resolve IDs. Returns (df, unmapped)."""
    vcfg = cfg["variants"]
    keys = df.select("id", "_risk_allele", "_snp_current", "_merged").unique()
    rows = []
    for snps, risk, cur, merged in keys.iter_rows():
        for tok, nid, allele, kind in split_tokens(snps, risk, cur, merged):
            rows.append((snps, risk, cur, merged, tok, nid, allele, kind))
    tok_df = pl.DataFrame(rows, orient="row", schema={
        "id": pl.String, "_risk_allele": pl.String, "_snp_current": pl.String,
        "_merged": pl.String, "_id_reported": pl.String, "_nid": pl.String,
        "_allele": pl.String, "variant_id_source": pl.String})
    if not vcfg["use_current_rsid_for_merged"]:
        tok_df = tok_df.with_columns(
            pl.when(pl.col("variant_id_source") == "merged_updated")
              .then(pl.col("_id_reported").str.to_lowercase()).otherwise(pl.col("_nid"))
              .alias("_nid"),
            pl.col("variant_id_source").replace("merged_updated", "rsid"))
    df = df.join(tok_df, on=["id", "_risk_allele", "_snp_current", "_merged"], how="inner")
    df = df.with_columns(pl.lit("").alias("rescue_build"))

    # --- positional rescue -------------------------------------------------
    pos = (df.filter(pl.col("variant_id_source") == "positional")
             .select("PUBMEDID", "_id_reported", "_allele").unique())
    rescued = {}
    if vcfg["rescue_positional"] and pos.height:
        items = []
        for pmid, tok, allele in pos.iter_rows():
            chrom, p, ref, alt = parse_positional(tok)
            items.append({"pmid": pmid, "token": tok, "chrom": chrom, "pos": p,
                          "ref": ref, "alt": alt, "allele": allele})
        log.info("rescuing %d positional IDs from %d papers", len(items),
                 pos["PUBMEDID"].n_unique())
        rescued = Rescuer(vcfg, base, log).resolve(items)
    if rescued:
        res = pl.DataFrame(
            [(pm, tok, rs, b, ",".join(g), r) for (pm, tok), (rs, b, g, r) in rescued.items()],
            orient="row", schema={"PUBMEDID": pl.String, "_id_reported": pl.String,
                                  "_rs": pl.String, "_build": pl.String,
                                  "_genes": pl.String, "_reason": pl.String})
        df = df.join(res, on=["PUBMEDID", "_id_reported"], how="left")
        ok = pl.col("_rs").fill_null("") != ""
        df = df.with_columns(
            pl.when(ok).then(pl.col("_rs")).otherwise(pl.col("_nid")).alias("_nid"),
            pl.when(ok).then(pl.lit("positional_rescued"))
              .otherwise(pl.col("variant_id_source")).alias("variant_id_source"),
            pl.when(ok).then(pl.col("_genes")).otherwise(pl.col("gene_name")).alias("gene_name"),
            pl.when(ok).then(pl.when(pl.col("_genes") != "").then(pl.lit(0.0))
                             .otherwise(pl.lit(1.0)))
              .otherwise(pl.col("intergenic")).alias("intergenic"),
            pl.when(ok & (pl.col("chromosome").fill_null("") == ""))
              .then(pl.col("_id_reported").map_elements(
                  lambda t: parse_positional(t)[0], return_dtype=pl.String))
              .otherwise(pl.col("chromosome")).alias("chromosome"),
            pl.when(ok).then(pl.col("_build")).otherwise(pl.lit("")).alias("rescue_build"),
        )
    else:
        df = df.with_columns(pl.lit(None, pl.String).alias("_reason"))

    # --- what could not be mapped ------------------------------------------
    keep_kinds = ["rsid", "merged_updated", "synonym", "positional_rescued"]
    lost = df.filter(~pl.col("variant_id_source").is_in(keep_kinds))
    unmapped = (lost.group_by("PUBMEDID", "_id_reported", "variant_id_source")
                    .agg(pl.col("_reason").first().fill_null(""),
                         pl.col("_allele").first(),
                         pl.col("id_in_source").n_unique().alias("n_studies"),
                         pl.col("phenotype").first(),
                         pl.col("p-value").first())
                    .rename({"_id_reported": "reported_id", "variant_id_source": "id_type",
                             "_reason": "rescue_status", "_allele": "risk_allele"})
                    .sort("PUBMEDID", "reported_id"))
    kinds = dict(df.group_by("variant_id_source").agg(pl.col("_id_reported").n_unique()).rows())
    reasons = dict(lost.filter(pl.col("_reason").is_not_null())
                       .group_by("_reason").agg(pl.col("_id_reported").n_unique()).rows())
    meta["variant_ids"] = {"unique_ids_by_type": kinds, "positional_not_rescued": reasons}
    log.info("variant IDs by type (unique): %s", kinds)
    if reasons:
        log.info("positional IDs not rescued: %s", reasons)

    if vcfg["keep_unresolved_positional"]:
        keep = pl.col("variant_id_source").is_in(keep_kinds) | (
            pl.col("variant_id_source") == "positional")
        df = df.with_columns(
            pl.when(pl.col("variant_id_source") == "positional")
              .then(pl.col("_id_reported").map_elements(
                  lambda t: "chr{}:{}".format(*parse_positional(t)[:2]), return_dtype=pl.String))
              .otherwise(pl.col("_nid")).alias("_nid"))
    else:
        keep = pl.col("variant_id_source").is_in(keep_kinds)
    df = (df.filter(keep)
            .drop("id", "_nid_tmp", strict=False)
            .rename({"_nid": "id"})
            .drop("_risk_allele", "_snp_current", "_merged", "_allele",
                  "_rs", "_build", "_genes", "_reason", strict=False))
    return df, unmapped


# --------------------------------------------------------------------------
# Regions
# --------------------------------------------------------------------------
def region_mapper(regions_cfg: dict):
    country_to_region = {}
    for region, countries in regions_cfg["regions"].items():
        for c in countries:
            country_to_region[str(c).strip().lower()] = region
    aliases = {k.strip().lower(): v.strip().lower()
               for k, v in (regions_cfg.get("aliases") or {}).items()}
    fallback = [(k.lower(), v) for k, v in regions_cfg["ancestry_fallback"].items()]

    def map_one(origin: str, ancestry: str) -> str:
        found = set()
        for c in (origin or "").split(","):
            key = c.strip().lower()
            key = aliases.get(key, key)
            if key in country_to_region:
                found.add(country_to_region[key])
        if found:
            return ", ".join(sorted(found))
        a = (ancestry or "").lower()
        for k, v in fallback:
            if k in a:
                return v
        return ""

    return map_one, set(country_to_region) | set(aliases)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--unmapped", required=True)
    args = ap.parse_args()

    base = Path.cwd()                                # work directory: caches, inputs
    code = Path(__file__).resolve().parent.parent    # workflow directory: config files
    cfg = load_yaml(args.config)
    regions_cfg = load_yaml(code / cfg["regions"])
    cat = cfg["catalog"]
    t_start = time.time()
    session = http_session()

    # ---- 1. Fetch (into memory, checked) ------------------------------------
    meta = {"release": args.release, "run_date": date.today().isoformat(), "sources": {}}
    raws = {}
    for key in ("ancestry", "associations"):
        local = cat.get(f"local_{key}")
        if local:
            p = Path(local) if Path(local).is_absolute() else base / local
            raws[key] = p.read_bytes()
            meta["sources"][key] = {"path": str(p)}
            log.info("read local %s (%.1f MB)", p, len(raws[key]) / 1e6)
        else:
            url = f"{cat['base_url']}/{cat[key]}"
            raws[key], info = fetch_checked(url, session, is_zip=url.endswith(".zip"))
            if info["release"] != args.release:
                log.warning("%s on the FTP is dated %s, not %s (release changed mid-run?)",
                            key, info["release"], args.release)
            meta["sources"][key] = info

    # ---- 2. Ancestry + associations ---------------------------------------
    ances, route_counts = load_ancestry(raws.pop("ancestry"), cfg, regions_cfg)
    assoc, n_assoc_total = load_associations(
        raws.pop("associations"), set(ances["id_in_source"].unique().to_list()))
    meta["counts"] = {"associations_total": n_assoc_total}
    meta["studies_by_inclusion_route"] = route_counts

    assoc = assoc.with_columns(
        pl.col("publication_year").str.split("-").list.first(),
        pl.col("intergenic").cast(pl.Float64, strict=False).fill_null(0.0),
    )

    # ---- 3. Merge, normalise variant IDs, explode genes --------------------
    df = ances.join(assoc, on="id_in_source", how="inner").filter(pl.col("id") != "")
    log.info("merged: %d rows", df.height)
    df, unmapped = normalise_ids(df, cfg, base, meta)
    log.info("after resolving variant IDs: %d rows", df.height)

    # Annotate genes before exploding, so rescued variants can keep only the
    # NAMED genes they overlap. Ensembl overlap also returns unnamed "novel
    # transcript" genes, which the GWAS Catalog's own mapping does not report;
    # they are kept only when no named gene overlaps (the variant is still genic).
    df = df.with_columns(pl.col("gene_name").fill_null("").str.replace_all(";", ","))
    ensg = sorted({g.strip() for s in df["gene_name"].unique().to_list()
                   for g in s.split(",") if g.strip().startswith("ENSG")})
    genes, ann_status = annotate_genes(ensg, cfg, base, log)
    meta["annotation"] = ann_status
    named = set(genes.filter(pl.col("gene_symbol_external") != "")["gene_name"].to_list())
    rescued = pl.col("variant_id_source") == "positional_rescued"
    lists = df.filter(rescued)["gene_name"].unique().to_list()

    def keep_named(s):
        gs = [g.strip() for g in s.split(",") if g.strip()]
        return ",".join([g for g in gs if g in named] or gs)

    mapping = {s: keep_named(s) for s in lists}
    df = df.with_columns(pl.when(rescued).then(pl.col("gene_name").replace(mapping))
                         .otherwise(pl.col("gene_name")).alias("gene_name"))
    n_dropped = sum(len(s.split(",")) - len(v.split(",")) for s, v in mapping.items() if s)
    log.info("rescued variants: kept named genes only where available "
             "(%d unnamed gene links dropped across %d gene lists)", n_dropped, len(mapping))

    # Rescued variants have no GWAS Catalog MAPPED_GENE, so build one in the
    # Catalog's style: the named genes they fall in ("A, B"), otherwise the
    # nearest named genes on either side ("UP - DOWN").
    vcfg = cfg["variants"]
    sym = dict(genes.select("gene_name", "gene_symbol_external").rows())
    resc = df.filter(rescued).select("gene_name", "rescue_build", "_id_reported").unique()
    mapped, flank_keys = {}, {}
    for gl, b, tok in resc.iter_rows():
        names = [sym[g] for g in gl.split(",") if g and sym.get(g)]
        if names:
            mapped[(gl, b, tok)] = ", ".join(dict.fromkeys(names))
        else:
            chrom, pos = parse_positional(tok)[:2]
            flank_keys[(gl, b, tok)] = (b, chrom, pos)
    flank_pending = 0
    if flank_keys and vcfg.get("flank_genes", True):
        rescuer = Rescuer(vcfg, base, log)
        res = rescuer.flanks(sorted(set(flank_keys.values())), vcfg.get("flank_window", 500000))
        for k, fk in flank_keys.items():
            if res.get(fk) is None:
                flank_pending += 1
            else:
                mapped[k] = res[fk]
        log.info("flanking genes for %d rescued positions without a named gene "
                 "(%d pending)", len(set(flank_keys.values())), flank_pending)
    meta.setdefault("variant_ids", {})["flank_lookup_pending"] = flank_pending
    if mapped:
        mg = pl.DataFrame([(*k, v) for k, v in mapped.items()], orient="row",
                          schema={"gene_name": pl.String, "rescue_build": pl.String,
                                  "_id_reported": pl.String, "_mapped": pl.String})
        df = df.join(mg, on=["gene_name", "rescue_build", "_id_reported"], how="left")
        df = df.with_columns(
            pl.when(rescued & (pl.col("curated_gene_symbol").fill_null("") == ""))
              .then(pl.col("_mapped").fill_null(""))
              .otherwise(pl.col("curated_gene_symbol")).alias("curated_gene_symbol")
        ).drop("_mapped")

    df = (df.with_columns(pl.col("gene_name").str.split(","))
            .explode("gene_name")
            .with_columns(pl.col("gene_name").str.strip_chars()))
    log.info("after exploding gene_name: %d", df.height)

    # ---- 4. Regions (computed once per unique country/ancestry pair) -------
    mapper, known = region_mapper(regions_cfg)
    pairs = df.select("origin_of_participants", "broad_ancestral_category").unique()
    pairs = pairs.with_columns(pl.Series(
        "geographical_region", [mapper(o, a) for o, a in pairs.iter_rows()], dtype=pl.String))
    df = df.join(pairs, on=["origin_of_participants", "broad_ancestral_category"],
                 how="left", nulls_equal=True)
    unknown = sorted({t.strip() for o in pairs["origin_of_participants"].to_list()
                      for t in (o or "").split(",") if t.strip()
                      and t.strip().lower() not in known})
    meta["countries_without_region"] = unknown

    # ---- 5. Gene annotation (looked up above) --------------------------------
    ensg = sorted(g for g in df["gene_name"].unique().to_list() if g and g.startswith("ENSG"))
    genes = genes.filter(pl.col("gene_name").is_in(ensg))
    df = df.join(genes, on="gene_name", how="left")
    sym = pl.col("gene_symbol_external").fill_null("")
    df = df.with_columns(
        pl.when(sym != "").then(sym)
          .when(pl.col("curated_gene_symbol").fill_null("") != "").then(pl.col("curated_gene_symbol"))
          .otherwise(pl.col("gene_symbol_mined").fill_null(""))
          .alias("gene_symbol_external"),
        pl.when(pl.col("gene_name").str.to_lowercase().is_in(["nan", "null", "none"]))
          .then(pl.lit("")).otherwise(pl.col("gene_name")).alias("gene_name"),
    )

    # ---- 5b. Variant consequence for variants without a Catalog CONTEXT ------
    # (all rescued positional variants, plus a few others): Ensembl VEP.
    no_ctx = (pl.col("context").fill_null("") == "") & pl.col("id").str.starts_with("rs")
    vep_filled = vep_pending = 0
    if vcfg.get("vep_context", True):
        ids = df.filter(no_ctx)["id"].unique().to_list()
        if ids:
            csq, vep_pending = consequences(ids, vcfg, base, log)
            before = df.filter(no_ctx).height
            df = df.with_columns(
                pl.when(no_ctx).then(pl.col("id").replace_strict(csq, default=pl.col("context")))
                  .otherwise(pl.col("context")).alias("context"))
            vep_filled = before - df.filter(no_ctx).height
            log.info("context from Ensembl VEP: %d records (%d unique variants looked up, "
                     "%d pending)", vep_filled, len(ids), vep_pending)
    meta["variant_ids"]["context_from_vep"] = vep_filled
    meta["variant_ids"]["vep_lookup_pending"] = vep_pending

    # ---- 6. Metadata and transformations -----------------------------------
    run_tag = date.today().strftime("%d%m%Y")
    first_number = r"(\d+\.?\d*)"
    df = df.with_columns(
        pl.lit("GWAS").alias("study_type"),
        pl.lit(run_tag).alias("data_ac"),
        pl.col("PUBMEDID").alias("publication"),
        pl.lit("SNV").alias("variant_type"),
        pl.lit("PubMed").alias("publication_type"),
        pl.lit("GWAS Catalog").alias("source"),
        # "True"/"False" text, as the June 2026 (pandas) export wrote it
        pl.when(pl.col("origin_of_participants").fill_null("").str.contains(",", literal=True))
          .then(pl.lit("True")).otherwise(pl.lit("False")).alias("mixed_population"),
        pl.when(pl.col("intergenic") == 1.0).then(pl.lit("intergenic"))
          .otherwise(pl.lit("non-intergenic")).alias("intergenic"),
        pl.col("replication_sample_size").fill_null("").str.replace_all(",", "")
          .str.extract(first_number).fill_null(""),
        pl.col("initial_sample_size").fill_null("").str.replace_all(",", "")
          .str.extract(first_number).fill_null(""),
    )
    ig = pl.col("intergenic") == "intergenic"
    df = df.with_columns(
        pl.when(ig).then(pl.lit("unresolved")).otherwise(pl.col("uniprot")).alias("uniprot"),
        pl.when(ig).then(pl.lit("unresolved")).otherwise(pl.col("function")).alias("function"),
    )
    str_cols = [c for c, t in df.schema.items() if t == pl.String]
    df = df.with_columns(pl.col(str_cols).fill_null(""))
    f, u = pl.col("function"), pl.col("uniprot")
    df = df.with_columns(
        pl.when((f != "") & (f != "unresolved") & (u == "")).then(pl.lit("unresolved"))
          .otherwise(u).alias("uniprot"),
        pl.when((u != "") & (u != "unresolved") & (f == "")).then(pl.lit("unresolved"))
          .otherwise(f).alias("function"),
    )

    # ---- 7. Deduplicate keeping the lowest p-value -------------------------
    before = df.height
    df = df.with_columns(pl.col("p-value").cast(pl.Float64, strict=False))
    subset = [c for c in df.columns if c != "p-value" and c not in INTERNAL]
    df = (df.sort("p-value", nulls_last=True, maintain_order=True)
            .unique(subset=subset, keep="first", maintain_order=True))
    log.info("deduplication: %d -> %d", before, df.height)

    # ---- 8. New vs previously imported -----------------------------------
    prev_files = (cfg.get("compare") or {}).get("previous_imports") or cfg.get("previous_imports") or []
    prev = list(read_previous_ids(prev_files, base, log)) if prev_files else []
    meta["compared_with_previous_imports"] = [str(p) for p in prev_files]
    if not prev_files:
        log.info("no previous imports given: every record is marked new")
    df = df.with_columns(
        ~(pl.col("id").is_in(prev) | pl.col("_id_reported").str.to_lowercase().is_in(prev))
        .alias("is_new"))

    # ---- 9. Write -----------------------------------------------------------
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(out, compression="zstd", compression_level=10)
    unmapped.write_csv(args.unmapped, separator="\t")

    eth = df["ethnicity"] != ""
    meta["counts"].update({
        "records": df.height,
        "unique_variants": df["id"].n_unique(),
        "studies": df["id_in_source"].n_unique(),
        "publications": df["PUBMEDID"].n_unique(),
        "intergenic": int((df["intergenic"] == "intergenic").sum()),
        "non_intergenic": int((df["intergenic"] != "intergenic").sum()),
        "new_records": int(df["is_new"].sum()),
        "with_ethnicity": int(eth.sum()),
        "explicit_ancestry_records": int((df["ancestry_status"] == "explicit").sum()),
        "unresolved_ancestry_records": int((df["ancestry_status"] == "unresolved").sum()),
        "rescued_variant_records": int((df["variant_id_source"] == "positional_rescued").sum()),
        "unmapped_variant_ids": unmapped.height,
        "ensg_genes": len(ensg),
        "genes_with_symbol": int(genes.filter(pl.col("gene_symbol_external") != "").height),
        "genes_with_uniprot": int(genes.filter(pl.col("uniprot") != "").height),
    })
    meta["runtime_seconds"] = round(time.time() - t_start)
    Path(args.meta).write_text(json.dumps(meta, indent=2))
    log.info("wrote %s (%.1f MB) in %ds", out, out.stat().st_size / 1e6,
             meta["runtime_seconds"])
    for k, v in meta["counts"].items():
        log.info("  %-28s %s", k, f"{v:,}")
    if unknown:
        log.info("%d recruitment countries have no AGMP region (non-African ones are "
                 "expected; listed in release.json under countries_without_region)",
                 len(unknown))


if __name__ == "__main__":
    main()
