"""Quality checks on the curated table. Export runs only if every FAIL check passes.

Writes qc.tsv (check, value, threshold, status) and exits non-zero on failure,
so a truncated download, a changed GWAS Catalog format or an annotation outage
cannot silently produce a bad AGMP import.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from common import get_logger, load_yaml  # noqa: E402
from export import OUTPUT_COLUMNS  # noqa: E402

log = get_logger("validate")
INTERNAL = {"_id_reported", "variant_id_source", "inclusion_route", "rescue_build"}


def previous_release(path: str):
    """The release to compare with: a results/<release> folder or its release.json."""
    p = Path(path)
    if p.is_dir():
        p = p / "release.json"
    if not p.exists():
        raise SystemExit(f"compare.previous_release: {p} not found. Give the folder of an "
                         "earlier release (results/<date>) or its release.json.")
    return json.loads(p.read_text())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    cfg = load_yaml(args.config)
    q = cfg["qc"]
    meta = json.loads(Path(args.meta).read_text())
    df = pl.read_parquet(args.table)
    n = df.height
    checks = []

    def check(name, value, threshold, ok, level="FAIL"):
        checks.append({"check": name, "value": str(value), "threshold": str(threshold),
                       "status": "PASS" if ok else level})

    missing = [c for c in OUTPUT_COLUMNS if c not in df.columns]
    check("output columns present", f"{len(OUTPUT_COLUMNS) - len(missing)}/{len(OUTPUT_COLUMNS)}",
          "all", not missing)
    total = meta["counts"]["associations_total"]
    check("associations in source file", total, f">= {q['min_associations_total']}",
          total >= q["min_associations_total"])
    check("curated records", n, f">= {q['min_records']}", n >= q["min_records"])

    id_re = (r"^(rs\d+|chr[0-9XYMT]+:\d+)$" if cfg["variants"]["keep_unresolved_positional"]
             else r"^rs\d+$")
    bad_ids = df.filter(~pl.col("id").str.contains(id_re)).height
    check("variant IDs well formed", bad_ids, "0 malformed", bad_ids == 0)

    subset = [c for c in df.columns if c not in INTERNAL and c != "is_new"]
    dups = n - df.select(subset).unique().height
    check("duplicate records", dups, "0", dups == 0)

    bad_ig = df.filter(~pl.col("intergenic").is_in(["intergenic", "non-intergenic"])).height
    check("intergenic labels valid", bad_ig, "0 invalid", bad_ig == 0)

    no_region = df.filter(pl.col("geographical_region") == "").height
    check("records with a geographical region", f"{n - no_region}/{n}", "all", no_region == 0)

    num_p = df["p-value"].is_not_null().mean() if n else 0
    check("numeric p-values", f"{num_p:.4f}", f">= {q['min_numeric_pvalue_fraction']}",
          num_p >= q["min_numeric_pvalue_fraction"])

    prev_path = (cfg.get("compare") or {}).get("previous_release")
    prev = previous_release(prev_path) if prev_path else None
    if prev:
        p_n = prev["counts"]["records"]
        drop = (p_n - n) / p_n if p_n else 0
        check(f"record change vs {prev['release']}", f"{p_n:,} -> {n:,} ({-drop:+.1%})",
              f"drop <= {q['max_record_drop_fraction']:.0%}",
              drop <= q["max_record_drop_fraction"])

    with_sym = df.filter(pl.col("gene_symbol_external") != "").height
    frac = with_sym / max(n, 1)
    check("records with a gene symbol", f"{with_sym:,}/{n:,} ({frac:.1%})",
          f">= {q.get('min_gene_symbol_fraction', 0.97):.0%}", frac >= q.get("min_gene_symbol_fraction", 0.97),
          "WARN")

    with_ctx = df.filter(pl.col("context").fill_null("") != "").height
    check("records with a variant consequence (context)", f"{with_ctx:,}/{n:,} ({with_ctx / max(n, 1):.1%})",
          f">= {q.get('min_context_fraction', 0.95):.0%}",
          with_ctx / max(n, 1) >= q.get("min_context_fraction", 0.95), "WARN")
    ann = meta.get("annotation", {})
    failed = ann.get("ensembl_failed", 0) + ann.get("uniprot_failed", 0)
    check("gene annotation lookups failed", failed, "0 (retried next run)", failed == 0, "WARN")
    vids = meta.get("variant_ids", {})
    pending = (vids.get("positional_not_rescued", {}).get("lookup_pending", 0)
               + vids.get("flank_lookup_pending", 0))
    check("positional IDs awaiting lookup", pending, "0 (continued next run)", pending == 0, "WARN")
    eth = meta["counts"]["with_ethnicity"]
    check("records with self-reported ethnicity", f"{eth:,} ({eth / max(n, 1):.1%})", "info",
          True, "INFO")
    for route, k in sorted(meta.get("studies_by_inclusion_route", {}).items()):
        check(f"studies included via '{route}'", k, "info", True, "INFO")
    for kind, k in sorted(meta.get("variant_ids", {}).get("unique_ids_by_type", {}).items()):
        check(f"unique variant IDs of type '{kind}'", k, "info", True, "INFO")

    out = pl.DataFrame(checks)
    out.write_csv(args.out, separator="\t")
    for c in checks:
        log.info("%-5s %-45s %s", c["status"], c["check"], c["value"])
    failures = [c for c in checks if c["status"] == "FAIL"]
    if failures:
        log.error("%d QC check(s) failed; export skipped. See %s", len(failures), args.out)
        sys.exit(1)


if __name__ == "__main__":
    main()
