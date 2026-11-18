"""NLP triage of unresolved GWAS Catalog studies for AGMP manual curation.

For every study with an unresolved field, the script gathers text evidence in
three tiers and scores it against an editable lexicon of African countries,
demonyms, regions, diaspora terms, ethnolinguistic groups and cohorts:

  T0  GWAS Catalog sample descriptions and study title (no network)
  T1  PubMed title + abstract (NCBI E-utilities, batched, cached)
  T2  Methods/participants sections of open-access full texts
      (Europe PMC, optional, only for studies still unresolved after T1)

Output is one ranked TSV: the curator queue. High-priority rows have strong
evidence and affect many variants; low-priority rows show no African evidence
and can be checked last (or excluded after a spot check).

Issues flagged per study:
  ancestry   broad ancestral category is only NR / Other
  origin     no African country or diaspora region could be assigned
  ethnicity  African study with no self-reported ethnicity / ethnolinguistic group
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import re
import sqlite3
import sys
import time
import xml.etree.ElementTree as ET
import zlib
from collections import Counter
from datetime import date
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from common import get_logger, http_session, load_yaml  # noqa: E402

log = get_logger("triage")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{}/fullTextXML"
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
METHODS_RE = re.compile(
    r"method|material|participant|subject|population|cohort|sample|patients|"
    r"study design|data source|recruit", re.I)
SENT_RE = re.compile(r"(?<=[.!?;])\s+(?=[A-Z0-9(])")
UNRESOLVED_REGIONS = {"", "NR", "Other", "African unspecified", "Sub-Saharan African"}
POSITIVE = {"african_country", "demonym", "diaspora", "african_region",
            "african_generic", "ethnolinguistic", "cohort"}
PLACE = {"african_country", "demonym"}


# --------------------------------------------------------------------------
# Lexicon matching
# --------------------------------------------------------------------------
class Lexicon:
    def __init__(self, path: Path):
        lex = pl.read_csv(path, separator="\t", infer_schema=False).fill_null("")
        self.info = {}
        pos, neg = [], []
        for t, cat, region, w in lex.iter_rows():
            self.info[t.lower()] = (t, cat, region, int(w or 0))
            (neg if cat == "negative" else pos).append(t)
        self.pos = self._compile(pos)
        self.neg = self._compile(neg)

    @staticmethod
    def _compile(terms):
        alts = []
        for t in sorted(set(terms), key=len, reverse=True):
            esc = re.escape(t)
            # Capitalised terms are proper nouns: match case-sensitively.
            alts.append(esc if t[:1].isupper() or t.isupper() else f"(?i:{esc})")
        return re.compile(r"(?<![\w-])(?:" + "|".join(alts) + r")(?![\w-])")

    def find(self, text: str, source: str):
        """Yield (term, category, region, weight, source, snippet)."""
        if not text:
            return
        # Blank out negative phrases (same length, so offsets still line up)
        # and quote snippets from the original text.
        masked = self.neg.sub(lambda m: " " * len(m.group()), text)
        pos = 0
        for msent in SENT_RE.split(masked):
            start = masked.index(msent, pos)
            pos = start + len(msent)
            sent = text[start:pos]
            for m in self.pos.finditer(msent):
                key = m.group().lower()
                t, cat, region, w = self.info.get(key, (m.group(), "", "", 0))
                snippet = (sent[:m.start()] + "[[" + m.group() + "]]" + sent[m.end():]).strip()
                if len(snippet) > 300:
                    s = max(0, m.start() - 140)
                    snippet = "..." + snippet[s:s + 290] + "..."
                yield t, cat, region, w, source, snippet


# --------------------------------------------------------------------------
# Literature cache (sqlite: abstracts + compressed methods text only)
# --------------------------------------------------------------------------
def open_cache(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE IF NOT EXISTS abstracts (pmid TEXT PRIMARY KEY, "
                "pmcid TEXT, title TEXT, abstract TEXT, fetched TEXT)")
    con.execute("CREATE TABLE IF NOT EXISTS fulltext (pmcid TEXT PRIMARY KEY, "
                "status TEXT, methods BLOB, fetched TEXT)")
    return con


def fetch_abstracts(pmids: list[str], con: sqlite3.Connection) -> dict:
    have = {r[0] for r in con.execute("SELECT pmid FROM abstracts")}
    todo = [p for p in pmids if p and p.isdigit() and p not in have]
    log.info("abstracts: %d needed, %d cached, %d to fetch", len(pmids),
             len(set(pmids) & have), len(todo))
    if todo:
        session = http_session()
        key = os.environ.get("NCBI_API_KEY", "")
        pause = 0.11 if key else 0.35
        params = {"db": "pubmed", "retmode": "xml", "tool": "AGMP_GWAS"}
        if key:
            params["api_key"] = key
        if os.environ.get("NCBI_EMAIL"):
            params["email"] = os.environ["NCBI_EMAIL"]
        today = date.today().isoformat()
        for i in range(0, len(todo), 200):
            batch = todo[i:i + 200]
            r = session.post(EUTILS, data={**params, "id": ",".join(batch)}, timeout=120)
            r.raise_for_status()
            rows = []
            for art in ET.fromstring(r.content).iter("PubmedArticle"):
                pmid = art.findtext(".//MedlineCitation/PMID") or ""
                title = "".join(art.find(".//ArticleTitle").itertext()) \
                    if art.find(".//ArticleTitle") is not None else ""
                abstract = " ".join("".join(a.itertext()) for a in art.iter("AbstractText"))
                pmcid = ""
                for aid in art.iter("ArticleId"):
                    if aid.get("IdType") == "pmc":
                        pmcid = (aid.text or "").strip()
                rows.append((pmid, pmcid, title, abstract, today))
            got = {r[0] for r in rows}
            rows += [(p, "", "", "", today) for p in batch if p not in got]
            con.executemany("INSERT OR REPLACE INTO abstracts VALUES (?,?,?,?,?)", rows)
            con.commit()
            log.info("  PubMed batch %d/%d", i // 200 + 1, (len(todo) + 199) // 200)
            time.sleep(pause)
    q = "SELECT pmid, pmcid, title, abstract FROM abstracts"
    return {r[0]: r[1:] for r in con.execute(q) if r[0] in set(pmids)}


def _methods_text(xml_bytes: bytes) -> str:
    root = ET.fromstring(xml_bytes)
    body = root.find(".//body")
    if body is None:
        return ""
    parts = []
    for sec in body.iter("sec"):
        title = (sec.findtext("title") or "") + " " + (sec.get("sec-type") or "")
        if METHODS_RE.search(title):
            parts.append(" ".join(" ".join(p.itertext()) for p in sec.iter("p")))
    return re.sub(r"\s+", " ", " ".join(parts))


def fetch_fulltexts(pmcids: list[str], con, workers: int) -> dict:
    have = {r[0] for r in con.execute("SELECT pmcid FROM fulltext")}
    todo = [p for p in pmcids if p and p not in have]
    today = date.today().isoformat()
    # Europe PMC answers HTTP 500 (not 404) for closed-access papers, so ask
    # which are open access first (100 IDs per search) and fetch only those.
    search = http_session()
    open_access = set()
    for i in range(0, len(todo), 100):
        batch = todo[i:i + 100]
        r = search.get(EPMC_SEARCH, timeout=120, params={
            "query": " OR ".join(f"PMCID:{p}" for p in batch),
            "resultType": "lite", "format": "json", "pageSize": 1000})
        r.raise_for_status()
        for res in r.json().get("resultList", {}).get("result", []):
            if res.get("isOpenAccess") == "Y" and res.get("pmcid"):
                open_access.add(res["pmcid"])
    closed = [p for p in todo if p not in open_access]
    con.executemany("INSERT OR REPLACE INTO fulltext VALUES (?,?,?,?)",
                    [(p, "not_open_access", b"", today) for p in closed])
    con.commit()
    todo = [p for p in todo if p in open_access]
    log.info("full texts: %d candidates, %d cached, %d closed access, %d to fetch",
             len(pmcids), len(set(pmcids) & have), len(closed), len(todo))
    session = http_session(total_retries=2)

    def get(pmcid):
        try:
            r = session.get(EPMC.format(pmcid), timeout=(15, 120))
            r.raise_for_status()
            return pmcid, "ok", zlib.compress(_methods_text(r.content).encode(), 9)
        except Exception as e:  # keep going; failures are retried next run
            return pmcid, f"error:{type(e).__name__}", b""

    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for n, (pmcid, status, blob) in enumerate(ex.map(get, todo), 1):
            if not status.startswith("error"):
                con.execute("INSERT OR REPLACE INTO fulltext VALUES (?,?,?,?)",
                            (pmcid, status, blob, today))
            if n % 25 == 0:
                con.commit()
                log.info("  full text %d/%d", n, len(todo))
    con.commit()
    out = {}
    for pmcid, status, blob in con.execute("SELECT pmcid, status, methods FROM fulltext"):
        if pmcid in set(pmcids) and status == "ok" and blob:
            out[pmcid] = zlib.decompress(blob).decode()
    return out


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------
def summarise(hits, issues, cfg):
    best = {}
    for t, cat, region, w, src, snip in hits:
        if cat not in POSITIVE:
            continue
        k = t.lower()
        if k not in best or w > best[k][3]:
            best[k] = (t, cat, region, w, src, snip)
    score = sum(v[3] for v in best.values())
    regions = Counter(v[2] for v in best.values() if v[2])
    places = sorted({v[0] for v in best.values() if v[1] in PLACE})
    ethno = sorted({v[0] for v in best.values() if v[1] == "ethnolinguistic"})
    cohorts = sorted({v[0] for v in best.values() if v[1] == "cohort"})
    sources = sorted({v[4] for v in best.values()})

    rank = {"high": 0, "medium": 1, "low": 2}
    prios = []
    if "ancestry" in issues:
        prios.append("high" if score >= cfg["high_score"]
                     else "medium" if score >= cfg["medium_score"] else "low")
    if "origin" in issues:
        prios.append("high" if places else "medium" if regions else "low")
    if "ethnicity" in issues:
        prios.append("high" if ethno else "medium" if places or cohorts else "low")
    priority = min(prios, key=rank.get) if prios else "low"

    seen, evidence = set(), []
    for v in sorted(best.values(), key=lambda v: -v[3]):
        if v[5] not in seen:
            seen.add(v[5])
            evidence.append(f"[{v[4]}] {v[5]}")
        if len(evidence) == 3:
            break
    return {
        "priority": priority, "score": score,
        "candidate_regions": "; ".join(r for r, _ in regions.most_common(3)),
        "candidate_countries": "; ".join(places),
        "candidate_ethnolinguistic": "; ".join(ethno),
        "cohorts": "; ".join(cohorts),
        "evidence_sources": "; ".join(sources),
        "evidence": " | ".join(evidence),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    base = Path.cwd()                                # work directory: caches
    code = Path(__file__).resolve().parent.parent    # workflow directory: lexicon
    cfg = load_yaml(args.config)
    tcfg = cfg["triage"]
    t0 = time.time()

    df = pl.read_parquet(args.table, columns=[
        "id_in_source", "PUBMEDID", "id", "title", "broad_ancestral_category",
        "origin_of_participants", "geographical_region", "ethnicity",
        "ancestry_status", "initial_sample_desc", "replication_sample_desc"])
    # Flag each record, then group by publication: curators resolve a paper
    # once, and one paper can have thousands of GWAS Catalog accessions.
    # Issues are judged per study accession: each variant is repeated once per
    # ancestry row of its study, so an NR row next to an African row adds no
    # new variant and needs no curation.
    explicit = (pl.col("ancestry_status") == "explicit")
    region_unresolved = (pl.col("geographical_region").str.split(",")
                         .list.eval(pl.element().str.strip_chars()
                                    .is_in(list(UNRESOLVED_REGIONS)))
                         .list.all())
    df = df.with_columns(
        explicit.any().over("id_in_source").not_().alias("_anc"),
        region_unresolved.all().over("id_in_source").alias("_origin"),
        (explicit.any().over("id_in_source")
         & (pl.col("ethnicity") == "").all().over("id_in_source")).alias("_eth"),
    )
    join = lambda c: pl.col(c).filter(pl.col(c) != "").unique().sort().str.join(" | ")  # noqa: E731
    affected = lambda flag: pl.col("id").filter(flag).n_unique()  # noqa: E731
    pubs = df.group_by("PUBMEDID").agg(
        pl.col("title").first(),
        pl.col("id_in_source").n_unique().alias("n_studies"),
        pl.col("id_in_source").unique().sort().head(5).str.join(", ").alias("study_accessions"),
        pl.len().alias("n_records"),
        pl.col("id").n_unique().alias("n_variants"),
        affected(pl.col("_anc")).alias("n_variants_ancestry"),
        affected(pl.col("_origin")).alias("n_variants_origin"),
        affected(pl.col("_eth")).alias("n_variants_ethnicity"),
        affected(pl.col("_anc") | pl.col("_origin") | pl.col("_eth")).alias("n_variants_affected"),
        join("broad_ancestral_category").alias("broad_ancestral_category"),
        join("origin_of_participants").alias("origin_of_participants"),
        join("geographical_region").alias("geographical_region"),
        join("initial_sample_desc").alias("initial_sample_desc"),
        join("replication_sample_desc").alias("replication_sample_desc"),
    )

    rows = []
    for r in pubs.iter_rows(named=True):
        r["issues"] = [k for k in ("ancestry", "origin", "ethnicity")
                       if r[f"n_variants_{k}"] > 0]
        if r["issues"]:
            rows.append(r)
    log.info("%d of %d publications have at least one unresolved field",
             len(rows), pubs.height)

    lex = Lexicon(code / tcfg["lexicon"])
    con = open_cache(base / tcfg["cache"])
    abstracts = fetch_abstracts(sorted({r["PUBMEDID"] for r in rows}), con)

    def collect(r, fulltext=None):
        hits = list(lex.find(r["initial_sample_desc"], "catalog"))
        hits += lex.find(r["replication_sample_desc"], "catalog")
        hits += lex.find(r["title"], "catalog")
        pmcid, title, abstract = abstracts.get(r["PUBMEDID"], ("", "", ""))
        hits += lex.find(f"{title} {abstract}", "abstract")
        if fulltext and pmcid in fulltext:
            hits += lex.find(fulltext[pmcid], "fulltext")
        return pmcid, hits

    results = {}
    for r in rows:
        pmcid, hits = collect(r)
        results[r["PUBMEDID"]] = (pmcid, summarise(hits, r["issues"], tcfg))

    if tcfg.get("fulltext"):
        need = sorted((r for r in rows if results[r["PUBMEDID"]][1]["priority"] != "high"
                       and results[r["PUBMEDID"]][0]),
                      key=lambda r: -r["n_variants_affected"])
        need = need[: tcfg["fulltext_max_papers"]]
        texts = fetch_fulltexts([results[r["PUBMEDID"]][0] for r in need], con,
                                tcfg["fulltext_workers"])
        for r in need:
            pmcid, hits = collect(r, texts)
            results[r["PUBMEDID"]] = (pmcid, summarise(hits, r["issues"], tcfg))
    con.close()

    out_rows = []
    for r in rows:
        pmcid, s = results[r["PUBMEDID"]]
        out_rows.append({
            "priority": s["priority"], "issues": "; ".join(r["issues"]),
            "PUBMEDID": r["PUBMEDID"], "pmcid": pmcid,
            "n_variants_affected": r["n_variants_affected"],
            "n_variants_ancestry": r["n_variants_ancestry"],
            "n_variants_origin": r["n_variants_origin"],
            "n_variants_ethnicity": r["n_variants_ethnicity"],
            "n_studies": r["n_studies"], "study_accessions": r["study_accessions"],
            "score": s["score"],
            "candidate_regions": s["candidate_regions"],
            "candidate_countries": s["candidate_countries"],
            "candidate_ethnolinguistic": s["candidate_ethnolinguistic"],
            "cohorts": s["cohorts"], "evidence_sources": s["evidence_sources"],
            "evidence": s["evidence"],
            "broad_ancestral_category": r["broad_ancestral_category"],
            "origin_of_participants": r["origin_of_participants"],
            "geographical_region": r["geographical_region"],
            "title": r["title"],
            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{r['PUBMEDID']}/",
            "curator_decision": "", "curator_notes": "",
        })
    queue = pl.DataFrame(out_rows, infer_schema_length=None)
    queue = (queue.with_columns(pl.col("priority").replace_strict(
                {"high": 0, "medium": 1, "low": 2}, return_dtype=pl.Int8).alias("_p"))
                  .sort(["_p", "n_variants_affected"], descending=[False, True]).drop("_p"))
    queue.write_csv(args.out, separator="\t")

    summary = (queue.with_columns(pl.col("issues").str.split("; "))
                    .explode("issues").group_by("issues", "priority").len()
                    .sort("issues", "priority"))
    log.info("curation queue written: %s (%d publications, %.0fs)", args.out, queue.height,
             time.time() - t0)
    for issue, prio, n in summary.iter_rows():
        log.info("  %-10s %-7s %6d", issue, prio, n)
    Path(args.out).with_suffix(".summary.json").write_text(json.dumps(
        [dict(zip(["issue", "priority", "studies"], x)) for x in summary.iter_rows()],
        indent=1))


if __name__ == "__main__":
    main()
