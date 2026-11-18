"""Gene annotation (Ensembl symbol/description + UniProt accession/function).

Results are cached in a small parquet file, so later releases only query genes
that are new or older than `cache_max_age_days`.

Requests per cold run (~17k genes): ~18 Ensembl POSTs, 1 bulk UniProt
Swiss-Prot download (~7 MB), and a few ID-mapping jobs for genes that only
have unreviewed (TrEMBL) entries.
"""
from __future__ import annotations

import concurrent.futures as cf
import gzip
import json
import re
import time
from datetime import date, timedelta
from pathlib import Path

import polars as pl

from common import http_session

CACHE_SCHEMA = {
    "ensg": pl.String, "symbol": pl.String, "description": pl.String,
    "uniprot": pl.String, "function": pl.String, "fetched": pl.String,
}
SOURCE_RE = re.compile(r"\s*\[Source:.*\]$")


ENSEMBL = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}


def _ensembl_batch(session, ids, server="GRCh38"):
    r = session.post(
        f"{ENSEMBL[server]}/lookup/id",
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        json={"ids": ids}, timeout=120)
    r.raise_for_status()
    out = {}
    for gid, info in (r.json() or {}).items():
        if info:
            out[gid] = (info.get("display_name") or "",
                        SOURCE_RE.sub("", info.get("description") or ""))
    return out


def ensembl_lookup(ids, batch, workers, log):
    """Returns (results, failed_ids). A failed batch never stops the run."""
    session = http_session()
    batches = [ids[i:i + batch] for i in range(0, len(ids), batch)]
    out, failed = {}, []

    def safe(b):
        try:
            return b, _ensembl_batch(session, b), None
        except Exception as e:  # noqa: BLE001
            return b, {}, e

    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for i, (b, res, err) in enumerate(ex.map(safe, batches), 1):
            if err is not None:
                failed += b
                log.warning("  Ensembl batch %d/%d failed: %s", i, len(batches), err)
            else:
                out.update(res)
                log.info("  Ensembl batch %d/%d", i, len(batches))
    # IDs unknown on GRCh38 (retired genes; some rescued variants come from
    # GRCh37 papers): try the GRCh37 archive server.
    missing = [g for g in ids if g not in out and g not in set(failed)]
    if missing:
        found = 0
        for i in range(0, len(missing), batch):
            try:
                res = _ensembl_batch(session, missing[i:i + batch], "GRCh37")
                out.update(res)
                found += len(res)
            except Exception as e:  # noqa: BLE001
                log.warning("  Ensembl GRCh37 lookup failed: %s", e)
        log.info("  Ensembl GRCh37 fallback: %d of %d IDs found", found, len(missing))
    return out, failed


def _truncate(text, n):
    return text if len(text) <= n else text[: n - 3] + "..."


def swissprot_by_ensg(max_chars, log):
    """One streamed request: every reviewed human UniProt entry with its ENSG links."""
    session = http_session()
    url = ("https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=json"
           "&fields=accession,xref_ensembl,cc_function,protein_name"
           "&query=organism_id:9606+AND+reviewed:true")
    t0 = time.time()
    r = session.get(url, timeout=(30, 900))
    r.raise_for_status()
    data = json.loads(gzip.decompress(r.content))
    out = {}
    for e in data.get("results", []):
        acc = e.get("primaryAccession", "")
        func = ""
        for c in e.get("comments", []):
            if c.get("commentType") == "FUNCTION" and c.get("texts"):
                func = c["texts"][0].get("value", "")
                break
        if not func:
            func = (e.get("proteinDescription", {}).get("recommendedName", {})
                    .get("fullName", {}).get("value", ""))
        for x in e.get("uniProtKBCrossReferences", []):
            for p in x.get("properties", []):
                if p.get("key") == "GeneId":
                    ensg = p["value"].split(".")[0]
                    out.setdefault(ensg, (acc, _truncate(func, max_chars)))
    log.info("  UniProt Swiss-Prot: %d entries, %d ENSG links (%.0fs)",
             len(data.get("results", [])), len(out), time.time() - t0)
    return out


def _idmapping_job(session, ids):
    r = session.post("https://rest.uniprot.org/idmapping/run",
                     data={"ids": ",".join(ids), "from": "Ensembl", "to": "UniProtKB"},
                     timeout=120)
    r.raise_for_status()
    job = r.json()["jobId"]
    for _ in range(200):
        time.sleep(3)
        st = session.get(f"https://rest.uniprot.org/idmapping/status/{job}",
                         timeout=60, allow_redirects=False)
        if st.status_code == 303:
            break  # finished: status redirects to the results
        body = st.json() if st.content else {}
        status = body.get("jobStatus")
        if status == "ERROR" or st.status_code >= 400:
            raise RuntimeError(f"ID mapping job {job} failed: {str(body)[:200]}")
        if status not in ("RUNNING", "NEW", "QUEUED"):
            break
    else:
        raise TimeoutError(f"ID mapping job {job} did not finish in 10 minutes")
    res = session.get(
        f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job}",
        params={"format": "tsv", "fields": "accession,reviewed"}, timeout=300)
    res.raise_for_status()
    out = {}
    for line in res.text.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        ensg, acc, reviewed = parts[0].split(".")[0], parts[1], parts[2]
        if reviewed == "reviewed" or ensg not in out:
            out[ensg] = acc
    return out


def search_by_ensembl(ids, log, batch=40):
    """Fallback when ID mapping is down: UniProt search on Ensembl gene xrefs.

    Returns (mapping, failed_ids); reviewed entries are preferred.
    """
    session = http_session()
    out, failed = {}, []
    for i in range(0, len(ids), batch):
        b = ids[i:i + batch]
        query = "(" + " OR ".join(f"xref:ensembl-{g}" for g in b) + ") AND organism_id:9606"
        try:
            r = session.get("https://rest.uniprot.org/uniprotkb/stream", timeout=(30, 300),
                            params={"query": query, "format": "json",
                                    "fields": "accession,reviewed,xref_ensembl"})
            r.raise_for_status()
        except Exception as e:  # noqa: BLE001
            failed += b
            log.warning("  UniProt search fallback batch failed: %s", e)
            continue
        wanted = set(b)
        for e in r.json().get("results", []):
            acc = e.get("primaryAccession", "")
            reviewed = "Swiss-Prot" in e.get("entryType", "")
            for x in e.get("uniProtKBCrossReferences", []):
                for p in x.get("properties", []):
                    g = p["value"].split(".")[0] if p.get("key") == "GeneId" else ""
                    if g in wanted and (reviewed or g not in out):
                        out[g] = acc
    log.info("  UniProt search fallback: %d of %d genes mapped", len(out), len(ids))
    return out, failed


def idmapping(ids, batch, log):
    """Returns (mapping, failed_ids). A failed job never stops the run."""
    session = http_session()
    batches = [ids[i:i + batch] for i in range(0, len(ids), batch)]
    out, failed = {}, []

    def safe(b):
        try:
            return b, _idmapping_job(session, b), None
        except Exception as e:  # noqa: BLE001
            return b, {}, e

    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        for i, (b, res, err) in enumerate(ex.map(safe, batches), 1):
            if err is not None:
                failed += b
                log.warning("  UniProt ID mapping job %d/%d failed: %s", i, len(batches), err)
            else:
                out.update(res)
                log.info("  UniProt ID mapping job %d/%d", i, len(batches))
    return out, failed


def annotate_genes(ensg_ids, cfg, base: Path, log) -> pl.DataFrame:
    """Return gene_name, gene_symbol_external, function, uniprot for each ENSG id."""
    acfg = cfg["annotation"]
    cache_path = base / acfg["cache"]
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        cache = pl.read_parquet(cache_path)
    else:
        cache = pl.DataFrame(schema=CACHE_SCHEMA)

    cutoff = (date.today() - timedelta(days=acfg["cache_max_age_days"])).isoformat()
    # Entries with no symbol and no description are re-tried every run (cheap;
    # they may be found later, e.g. via the GRCh37 fallback).
    known = (pl.col("symbol") != "") | (pl.col("description") != "")
    fresh = set(cache.filter((pl.col("fetched") >= cutoff) & known)["ensg"].to_list())
    todo = sorted(set(ensg_ids) - fresh)
    log.info("genes: %d needed, %d cached, %d to query", len(ensg_ids),
             len(set(ensg_ids) & fresh), len(todo))

    status = {"queried": len(todo), "ensembl_failed": 0, "uniprot_failed": 0,
              "used_stale_cache": 0}
    if todo:
        ens, ens_failed = ensembl_lookup(todo, acfg["ensembl_batch"],
                                         acfg["ensembl_workers"], log)
        try:
            sp = swissprot_by_ensg(acfg["function_max_chars"], log)
            sp_ok = True
        except Exception as e:  # noqa: BLE001
            log.warning("  UniProt Swiss-Prot download failed (%s); using ID mapping only", e)
            sp, sp_ok = {}, False
        rest = [g for g in todo if g not in sp]
        mapped, map_failed = idmapping(rest, acfg["idmapping_batch"], log) if rest else ({}, [])
        if map_failed:  # ID mapping service down: fall back to UniProt search
            extra, map_failed = search_by_ensembl(map_failed, log)
            mapped.update(extra)
        failed = set(ens_failed) | set(map_failed)
        if not sp_ok:
            failed |= set(rest)
        status["ensembl_failed"] = len(ens_failed)
        status["uniprot_failed"] = len(set(map_failed) | (set(rest) if not sp_ok else set()))

        old = {r["ensg"]: r for r in cache.iter_rows(named=True)}
        today = date.today().isoformat()
        rows = []
        for g in todo:
            if g in failed and g in old:
                # keep the previous (stale) annotation and retry next run
                rows.append({**old[g], "fetched": ""})
                status["used_stale_cache"] += 1
                continue
            symbol, desc = ens.get(g, ("", ""))
            if g in sp:
                acc, func = sp[g]
            else:
                acc, func = mapped.get(g, ""), ""
            rows.append({"ensg": g, "symbol": symbol, "description": desc,
                         "uniprot": acc, "function": func or desc,
                         "fetched": "" if g in failed else today})
        new = pl.DataFrame(rows, schema=CACHE_SCHEMA)
        cache = pl.concat([cache.filter(~pl.col("ensg").is_in(todo)), new])
        tmp = cache_path.with_suffix(".tmp")
        cache.write_parquet(tmp, compression="zstd")
        tmp.replace(cache_path)  # atomic: a crash never leaves a half-written cache
        log.info("gene cache updated: %s (%d genes)", cache_path, cache.height)
        if failed:
            log.warning("%d genes could not be (fully) annotated this run; they are retried "
                        "next run", len(failed))

    genes = (cache.filter(pl.col("ensg").is_in(list(ensg_ids)))
             .select(pl.col("ensg").alias("gene_name"),
                     pl.col("symbol").alias("gene_symbol_external"),
                     "function", "uniprot"))
    return genes, status
