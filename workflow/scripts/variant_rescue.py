"""Resolve positional variant IDs (e.g. chr7:107194965) to dbSNP rsIDs.

Large recent studies (whole-genome sequencing, biobanks) often report rare
variants by position. The GWAS Catalog leaves these unmapped, so the June 2026
notebook dropped them (~50,000 unique IDs in the 2026-11-02 release).

The lookup is conservative, because dbSNP is so dense that a position alone
often matches some variant by chance:

1. Per publication, a sample of positions is looked up on BOTH GRCh38 and
   GRCh37 (Ensembl REST /overlap). The build whose variants match the
   reported risk allele clearly more often is taken as the paper's build.
   Papers with no clear build are left unresolved.
2. Each position is accepted only if exactly one rsID starts at that position
   AND lists the reported risk allele.
3. Genes overlapping the position are returned too, so rescued variants get a
   gene and a correct intergenic flag.

Lookups are cached in sqlite, so each position is queried once ever.
"""
from __future__ import annotations

import concurrent.futures as cf
import json
import re
import sqlite3
import threading
import time
from collections import defaultdict
from datetime import date
from pathlib import Path

from common import http_session

SERVERS = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
POS_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT|M)[:_](\d+)(?:[:_]([ACGT]+)[:_/]([ACGT]+))?$",
                    re.I)


def parse_positional(token: str):
    """Return (chrom, pos, ref, alt) for a positional ID, else None."""
    m = POS_RE.match(token.strip())
    if not m:
        return None
    chrom = m.group(1).upper()
    chrom = "MT" if chrom == "M" else chrom
    return chrom, int(m.group(2)), (m.group(3) or "").upper(), (m.group(4) or "").upper()


class InvalidRegion(Exception):
    """Ensembl rejected the position itself (HTTP 400), e.g. beyond the end of
    the chromosome on this build or an unknown chromosome. This is permanent,
    so it is cached as "no variants here" rather than retried every run."""


def _get_json(session, url, params, timeout):
    r = session.get(url, params=params, headers={"Content-Type": "application/json"},
                    timeout=timeout)
    if r.status_code == 400:
        raise InvalidRegion(r.text[:200])
    r.raise_for_status()
    return r.json()


class RateLimiter:
    def __init__(self, per_second: float):
        self.interval = 1.0 / per_second
        self.lock = threading.Lock()
        self.next = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            if self.next > now:
                time.sleep(self.next - now)
            self.next = max(now, self.next) + self.interval


class Rescuer:
    def __init__(self, vcfg: dict, base: Path, log):
        self.cfg = vcfg
        self.log = log
        path = base / vcfg["rescue_cache"]
        path.parent.mkdir(parents=True, exist_ok=True)
        self.con = sqlite3.connect(path, check_same_thread=False)
        self.con.execute("CREATE TABLE IF NOT EXISTS overlap (build TEXT, chrom TEXT, pos INT, "
                         "variants TEXT, genes TEXT, fetched TEXT, PRIMARY KEY(build, chrom, pos))")
        self.db_lock = threading.Lock()
        self.session = http_session(pool=vcfg["rescue_workers"] * 2)
        self.limiters = {b: RateLimiter(vcfg["rescue_rate_per_sec"]) for b in SERVERS}
        self.new_lookups = 0
        self.cap_hit = False
        self.errors = 0
        self.invalid = 0
        self.error_kinds = defaultdict(int)
        self.unsaved = 0

    def _failed(self, exc):
        self.errors += 1
        code = getattr(getattr(exc, "response", None), "status_code", None)
        self.error_kinds[f"HTTP {code}" if code else type(exc).__name__] += 1

    # ---------------------------------------------------------------- lookups
    def _cached(self, build, chrom, pos):
        with self.db_lock:
            row = self.con.execute("SELECT variants, genes FROM overlap WHERE build=? AND "
                                   "chrom=? AND pos=?", (build, chrom, pos)).fetchone()
        return (json.loads(row[0]), json.loads(row[1])) if row else None

    def _fetch(self, build, chrom, pos):
        self.limiters[build].wait()
        try:
            data = _get_json(self.session,
                             f"{SERVERS[build]}/overlap/region/human/{chrom}:{pos}-{pos}",
                             [("feature", "variation"), ("feature", "gene")], (15, 60))
        except InvalidRegion:
            self.invalid += 1
            return [], []  # cached: never matches, never retried
        variants, genes = [], []
        for e in data:
            if e.get("feature_type") == "variation":
                vid = str(e.get("id", ""))
                if vid.startswith("rs") and e.get("start") == pos:
                    variants.append([vid, e.get("alleles") or []])
            elif e.get("feature_type") == "gene" and str(e.get("id", "")).startswith("ENSG"):
                genes.append(e["id"])
        return variants, genes

    def lookup(self, build, chrom, pos):
        hit = self._cached(build, chrom, pos)
        if hit is not None:
            return hit
        if self.new_lookups >= self.cfg["rescue_max_new_lookups"]:
            self.cap_hit = True
            return None
        self.new_lookups += 1
        try:
            variants, genes = self._fetch(build, chrom, pos)
        except Exception as exc:  # noqa: BLE001
            self._failed(exc)
            return None  # not cached: retried next run
        with self.db_lock:
            self.con.execute("INSERT OR REPLACE INTO overlap VALUES (?,?,?,?,?,?)",
                             (build, chrom, pos, json.dumps(variants), json.dumps(genes),
                              date.today().isoformat()))
            self.unsaved += 1
            if self.unsaved >= 200:  # an interrupted run keeps its progress
                self.con.commit()
                self.unsaved = 0
        return variants, genes

    def _many(self, build, keys):
        with cf.ThreadPoolExecutor(max_workers=self.cfg["rescue_workers"]) as ex:
            res = list(ex.map(lambda k: self.lookup(build, k[0], k[1]), keys))
        with self.db_lock:
            self.con.commit()
        return dict(zip(keys, res))

    # --------------------------------------------------------------- matching
    @staticmethod
    def match(hit, ref, alt, allele):
        """Unique rsID at the position that carries the reported allele(s)."""
        if hit is None:
            return None
        variants, _ = hit
        wanted = {a for a in (ref, alt, allele) if a and a != "?"}
        if not wanted:
            return None
        ok = [vid for vid, alleles in variants if wanted <= set(alleles)]
        return ok[0] if len(set(ok)) == 1 else None

    # ---------------------------------------------------------- flanking genes
    def _flank_fetch(self, build, chrom, pos, window):
        self.limiters[build].wait()
        start = max(1, pos - window)
        try:
            data = _get_json(self.session,
                             f"{SERVERS[build]}/overlap/region/human/{chrom}:{start}-{pos + window}",
                             {"feature": "gene"}, (15, 90))
        except InvalidRegion:
            self.invalid += 1
            return "", ""
        up = down = None
        for g in data:
            name = g.get("external_name") or ""
            if not name or name.startswith("ENSG"):
                continue  # unnamed genes are not reported, as in GWAS Catalog MAPPED_GENE
            s, e = g.get("start", 0), g.get("end", 0)
            if e < pos and (up is None or e > up[1]):
                up = (name, e)
            elif s > pos and (down is None or s < down[1]):
                down = (name, s)
        return (up[0] if up else ""), (down[0] if down else "")

    def flanks(self, keys, window):
        """Nearest named gene up- and downstream for (build, chrom, pos) keys.

        Returns {key: "UP - DOWN"} in the GWAS Catalog MAPPED_GENE style.
        """
        self.con.execute("CREATE TABLE IF NOT EXISTS flank (build TEXT, chrom TEXT, pos INT, "
                         "up TEXT, down TEXT, fetched TEXT, PRIMARY KEY(build, chrom, pos))")

        def one(key):
            build, chrom, pos = key
            with self.db_lock:
                row = self.con.execute("SELECT up, down FROM flank WHERE build=? AND chrom=? "
                                       "AND pos=?", key).fetchone()
            if row is None:
                if self.new_lookups >= self.cfg["rescue_max_new_lookups"]:
                    self.cap_hit = True
                    return key, None
                self.new_lookups += 1
                try:
                    row = self._flank_fetch(build, chrom, pos, window)
                except Exception as exc:  # noqa: BLE001
                    self._failed(exc)
                    return key, None
                with self.db_lock:
                    self.con.execute("INSERT OR REPLACE INTO flank VALUES (?,?,?,?,?,?)",
                                     (*key, *row, date.today().isoformat()))
                    self.unsaved += 1
                    if self.unsaved >= 200:
                        self.con.commit()
                        self.unsaved = 0
            up, down = row
            return key, " - ".join(x for x in (up, down) if x)

        with cf.ThreadPoolExecutor(max_workers=self.cfg["rescue_workers"]) as ex:
            out = dict(ex.map(one, keys))
        with self.db_lock:
            self.con.commit()
        return out

    # ------------------------------------------------------------------- main
    def resolve(self, items):
        """items: list of dicts with pmid, token, chrom, pos, ref, alt, allele.

        Returns {(pmid, token): (rsid, build, [genes], reason)}.
        """
        by_pub = defaultdict(list)
        for it in items:
            by_pub[it["pmid"]].append(it)
        out = {}
        n_sample = self.cfg["build_detection_sample"]
        min_rate = self.cfg["build_min_match_rate"]
        t0 = time.time()
        pubs = sorted(by_pub, key=lambda p: -len(by_pub[p]))
        for i, pmid in enumerate(pubs, 1):
            its = sorted(by_pub[pmid], key=lambda x: (x["chrom"], x["pos"]))
            usable = [x for x in its if (x["allele"] and x["allele"] != "?") or x["alt"]]
            for x in its:
                if x not in usable:
                    out[(pmid, x["token"])] = ("", "", [], "no_allele")
            if not usable:
                continue
            step = max(1, len(usable) // n_sample)
            sample = usable[::step][:n_sample]
            keys = [(x["chrom"], x["pos"]) for x in sample]
            rates, incomplete = {}, False
            for build in SERVERS:
                hits = self._many(build, keys)
                incomplete |= any(h is None for h in hits.values())
                got = sum(1 for x in sample if self.match(hits[(x["chrom"], x["pos"])],
                                                          x["ref"], x["alt"], x["allele"]))
                rates[build] = got / len(sample)
            if incomplete:  # cap reached or lookups failed: decide on a later run
                for x in usable:
                    out[(pmid, x["token"])] = ("", "", [], "lookup_pending")
                continue
            best = max(rates, key=rates.get)
            other = min(rates, key=rates.get)
            if rates[best] < min_rate or rates[best] < 1.5 * rates[other]:
                for x in usable:
                    out[(pmid, x["token"])] = ("", "", [], "build_undetermined")
                continue
            hits = self._many(best, [(x["chrom"], x["pos"]) for x in usable])
            for x in usable:
                hit = hits[(x["chrom"], x["pos"])]
                rs = self.match(hit, x["ref"], x["alt"], x["allele"])
                if hit is None:
                    reason = "lookup_pending"
                elif rs:
                    reason = "rescued"
                else:
                    reason = "no_unique_match"
                out[(pmid, x["token"])] = (rs or "", best, hit[1] if hit and rs else [], reason)
            if i % 25 == 0 or i == len(pubs):
                self.log.info("  variant rescue: %d/%d papers, %d new lookups (%.0fs)",
                              i, len(pubs), self.new_lookups, time.time() - t0)
        if self.cap_hit:
            self.log.warning("variant rescue stopped at rescue_max_new_lookups=%d; the rest "
                             "are looked up on the next run", self.cfg["rescue_max_new_lookups"])
        if self.invalid:
            self.log.info("variant rescue: %d positions rejected by Ensembl as invalid for the "
                          "build (cached as no match)", self.invalid)
        if self.errors:
            self.log.warning("variant rescue: %d lookups failed (retried next run): %s",
                             self.errors, dict(self.error_kinds))
        return out
