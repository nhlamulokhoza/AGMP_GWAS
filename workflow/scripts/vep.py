"""Variant consequence (the `context` column) from Ensembl VEP, for variants
the GWAS Catalog has no CONTEXT for.

Rescued positional variants (and a few others) arrive without a GWAS Catalog
CONTEXT. Ensembl VEP gives the same Sequence Ontology terms the Catalog uses
(intron_variant, missense_variant, ...): its most severe consequence on the
primary GRCh38 assembly is used, across all transcripts (VEP's "pick" option is
NOT used: it reports the most severe consequence of one chosen transcript only).
As in the GWAS Catalog, a variant outside every gene is "intergenic_variant"
even when VEP calls it upstream/downstream (within 5 kb of a gene). On 600
variants that do have a Catalog CONTEXT, this agreed exactly with it for 91%.
Results are cached in sqlite, so each rsID is queried once ever; rsIDs Ensembl
does not know are cached as "" (no answer).
"""
from __future__ import annotations

import concurrent.futures as cf
import json
import sqlite3
import threading
import time
from datetime import date
from pathlib import Path

from common import http_session

URL = "https://rest.ensembl.org/vep/human/id"
PRIMARY = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}
CATALOG_STYLE = {"upstream_gene_variant": "intergenic_variant",
                 "downstream_gene_variant": "intergenic_variant"}


def _query(session, ids):
    r = session.post(URL,
                     headers={"Content-Type": "application/json", "Accept": "application/json"},
                     data=json.dumps({"ids": ids}), timeout=(15, 300))
    r.raise_for_status()
    out = {}
    for e in r.json():
        vid, csq = e.get("id") or e.get("input"), e.get("most_severe_consequence") or ""
        if not vid or not csq:
            continue
        csq = CATALOG_STYLE.get(csq, csq)
        # the primary assembly wins over patches and alternative haplotypes
        if e.get("seq_region_name") in PRIMARY or vid not in out:
            out[vid] = csq
    return out


def consequences(rsids, vcfg: dict, base: Path, log) -> tuple[dict, int]:
    """Return ({rsid: consequence}, n_pending) for the given rsIDs."""
    rsids = sorted({r for r in rsids if r and r.startswith("rs")})
    path = base / vcfg["rescue_cache"]
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path, check_same_thread=False)
    con.execute("CREATE TABLE IF NOT EXISTS vep (rsid TEXT PRIMARY KEY, consequence TEXT, "
                "fetched TEXT)")
    lock = threading.Lock()
    found = {}
    for i in range(0, len(rsids), 900):  # sqlite variable limit
        chunk = rsids[i:i + 900]
        q = f"SELECT rsid, consequence FROM vep WHERE rsid IN ({','.join('?' * len(chunk))})"
        found.update(con.execute(q, chunk).fetchall())
    todo = [r for r in rsids if r not in found]
    cap = vcfg.get("vep_max_new_lookups", 100000)
    if len(todo) > cap:
        log.warning("VEP: %d rsIDs to look up, capped at vep_max_new_lookups=%d; the rest "
                    "are looked up on the next run", len(todo), cap)
        todo = todo[:cap]
    batch = vcfg.get("vep_batch", 200)
    batches = [todo[i:i + batch] for i in range(0, len(todo), batch)]
    session = http_session(pool=8)
    failed, t0, done = 0, time.time(), 0

    def run(ids):
        nonlocal failed, done
        try:
            res = _query(session, ids)
        except Exception as exc:  # noqa: BLE001
            failed += len(ids)
            log.warning("VEP batch of %d failed (%s); retried next run", len(ids), exc)
            return
        rows = [(r, res.get(r, ""), date.today().isoformat()) for r in ids]
        with lock:
            con.executemany("INSERT OR REPLACE INTO vep VALUES (?,?,?)", rows)
            con.commit()
            found.update({r: c for r, c, _ in rows})
            done += 1
            if done % 25 == 0 or done == len(batches):
                log.info("  VEP: %d/%d batches (%.0fs)", done, len(batches), time.time() - t0)

    with cf.ThreadPoolExecutor(max_workers=vcfg.get("vep_workers", 4)) as ex:
        list(ex.map(run, batches))
    con.close()
    pending = len(rsids) - len(found)
    return {r: c for r, c in found.items() if c}, pending
