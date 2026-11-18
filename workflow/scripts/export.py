"""Write the AGMP import files from the curated parquet table.

Each output file is a filter of the one parquet table. Every file can be
written as CSV (optionally compressed), Excel, or both: export.csv and
export.excel in the settings list which files go in which format.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from columns import NAMES  # noqa: E402
from common import get_logger, load_yaml  # noqa: E402

log = get_logger("export")

OUTPUT_COLUMNS = NAMES  # the 30 standard AGMP columns, in order (scripts/columns.py)
ETHNICITY_COLUMNS = [
    "id_in_source", "PUBMEDID", "phenotype", "title",
    "broad_ancestral_category", "ethnicity",
    "country_of_origin", "origin_of_participants",
    "geographical_region", "initial_sample_size",
    "replication_sample_size",
]

IG = pl.col("intergenic") == "intergenic"
NEW = pl.col("is_new")
ETH = pl.col("ethnicity") != ""
SPLITS = {
    "intergenic": NEW & IG,
    "non_intergenic": NEW & ~IG,
    "combined": NEW,
    "combined_with_old": pl.lit(True),
    "all_intergenic": IG,
    "all_non_intergenic": ~IG,
    "all_intergenic_with_ethnicity": IG & ETH,
    "all_intergenic_without_ethnicity": IG & ~ETH,
    "all_non_intergenic_with_ethnicity": ~IG & ETH,
    "all_non_intergenic_without_ethnicity": ~IG & ~ETH,
    "combined_with_old_with_ethnicity": ETH,
    "combined_with_old_without_ethnicity": ~ETH,
}

# What each output file holds (used by the report, launcher and docs)
ALL_SPLITS_HELP = {
    "intergenic": "New intergenic records (needs --previous)",
    "non_intergenic": "New non-intergenic records (needs --previous)",
    "combined": "All new records (needs --previous)",
    "combined_with_old": "All records",
    "all_intergenic": "All intergenic records",
    "all_non_intergenic": "All non-intergenic records",
    "all_intergenic_with_ethnicity": "Intergenic records with an ethnicity description",
    "all_intergenic_without_ethnicity": "Intergenic records without one",
    "all_non_intergenic_with_ethnicity": "Non-intergenic records with an ethnicity description",
    "all_non_intergenic_without_ethnicity": "Non-intergenic records without one",
    "combined_with_old_with_ethnicity": "All records with an ethnicity description",
    "combined_with_old_without_ethnicity": "All records without one",
    "ethnicity_report": "One row per study population: ancestry, ethnicity, countries, sample sizes",
}


def write(df: pl.DataFrame, stem: Path, fmt: str):
    if fmt == "csv":
        path = stem.with_name(stem.name + ".csv")
        df.write_csv(path)
    elif fmt == "csv.gz":
        import gzip
        path = stem.with_name(stem.name + ".csv.gz")
        with gzip.open(path, "wb", compresslevel=6) as fh:
            df.write_csv(fh)
    elif fmt == "xlsx":
        path = stem.with_name(stem.name + ".xlsx")
        df.write_excel(path, autofit=False)
    else:
        raise ValueError(f"unknown export format: {fmt}")
    return path


NEW_ONLY = {"intergenic", "non_intergenic", "combined"}
ALL_FILES = list(SPLITS) + ["ethnicity_report"]
EXCEL_MAX_ROWS = 1_048_575  # Excel's limit, less the header row


def file_list(value, what: str) -> list[str]:
    """export.csv / export.excel: "all", "none", or a list of output files.

    A list may start from "all" and remove files with a minus sign:
    "all,-intergenic,-combined" is every file except those two.
    """
    if value in (None, "none", [], ["none"]):
        return []
    if value in ("all", ["all"]):
        return ALL_FILES
    if isinstance(value, str):
        value = [x.strip() for x in value.split(",") if x.strip()]
    if "all" in value:
        drop = {x[1:] for x in value if x.startswith("-")}
        unknown = [x for x in drop if x not in ALL_FILES]
        if unknown:
            raise SystemExit(f"Unknown {what} file(s): {', '.join(unknown)}\n"
                             f"Choose from: {', '.join(ALL_FILES)}")
        return [f for f in ALL_FILES if f not in drop]
    unknown = [x for x in value if x not in ALL_FILES]
    if unknown:
        raise SystemExit(f"Unknown {what} file(s): {', '.join(unknown)}\n"
                         f"Choose from: {', '.join(ALL_FILES)}")
    return [f for f in ALL_FILES if f in value]


def chosen_columns(cfg: dict) -> list[str]:
    """export.columns: "all" (default) or a list, kept in the standard AGMP order."""
    want = cfg.get("columns", "all")
    if want in (None, "all", ["all"]):
        return OUTPUT_COLUMNS
    if isinstance(want, str):
        want = [c.strip() for c in want.split(",") if c.strip()]
    unknown = [c for c in want if c not in OUTPUT_COLUMNS]
    if unknown:
        raise SystemExit(f"Unknown export column(s): {', '.join(unknown)}\n"
                         f"Available: {', '.join(OUTPUT_COLUMNS)}")
    return [c for c in OUTPUT_COLUMNS if c in want]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args()
    full = load_yaml(args.config)
    cfg = full["export"]
    columns = chosen_columns(cfg)
    as_csv = file_list(cfg.get("csv", "all"), "CSV")
    as_excel = file_list(cfg.get("excel", "none"), "Excel")
    gzip_csv = bool(cfg.get("gzip", False))
    compared = bool((full.get("compare") or {}).get("previous_imports")
                    or full.get("previous_imports"))

    lf = pl.scan_parquet(args.table)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for old in outdir.glob("*"):  # a re-export never leaves stale files behind
        if old.is_file():
            old.unlink()
    tag = Path(args.table).stem  # e.g. agmp_gwas_2026-11-02
    manifest = []
    log.info("columns (%d): %s", len(columns), ", ".join(columns))

    for name in ALL_FILES:
        formats = (["csv.gz" if gzip_csv else "csv"] if name in as_csv else []) + \
                  (["xlsx"] if name in as_excel else [])
        if not formats:
            continue
        if name in NEW_ONLY and not compared:
            log.info("  %-40s skipped: no previous imports to compare with", name)
            continue
        if name == "ethnicity_report":
            df = lf.select(ETHNICITY_COLUMNS).unique(maintain_order=True).collect()
        else:
            df = lf.filter(SPLITS[name]).select(columns).collect()
        for fmt in formats:
            if fmt == "xlsx" and df.height > EXCEL_MAX_ROWS:
                raise SystemExit(f"{name} has {df.height:,} rows, more than Excel allows "
                                 f"({EXCEL_MAX_ROWS:,}). Write it as CSV instead.")
            path = write(df, outdir / f"{tag}_{name}", fmt)
            manifest.append({"split": name, "format": fmt, "rows": df.height,
                             "file": path.name, "mb": round(path.stat().st_size / 1e6, 2)})
            log.info("  %-40s %9s rows  %s", name, f"{df.height:,}", path.name)

    if not manifest:
        log.warning("no output files were chosen (export.csv and export.excel are both empty)")
    pl.DataFrame(manifest, schema={"split": pl.String, "format": pl.String, "rows": pl.Int64,
                                   "file": pl.String, "mb": pl.Float64}
                 ).write_csv(args.manifest, separator="\t")


if __name__ == "__main__":
    main()
