"""Interactive HTML report for one release (results/<release>/report.html).

One self-contained file: no internet needed to open it, nothing to install.
The table is condensed into counts (a "cube" over a few dimensions), so the
charts can filter each other in the browser without shipping 400,000 rows.
A small random sample of rows is included to preview the column choice.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from columns import COLUMNS, NAMES, RESOURCES  # noqa: E402
from common import get_logger, load_yaml  # noqa: E402
from export import ALL_SPLITS_HELP, NEW_ONLY  # noqa: E402

log = get_logger("report")
TEMPLATE = Path(__file__).resolve().parent.parent / "templates" / "report.html"
TOP_CONTEXTS = 12
SAMPLE_ROWS = 60


def s(col):
    return pl.col(col).cast(pl.Utf8).fill_null("")


def main():
    ap = argparse.ArgumentParser()
    for a in ["--config", "--table", "--meta", "--qc", "--manifest", "--out"]:
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    cfg = load_yaml(args.config)
    meta = json.loads(Path(args.meta).read_text())
    out_dir = Path(args.out).parent
    df = pl.read_parquet(args.table)
    n = df.height

    # ---- dimensions for the cross-filtering charts ------------------------
    ctx_counts = df.group_by(s("context").alias("c")).len().sort("len", descending=True)
    top_ctx = [c for c in ctx_counts["c"].to_list() if c][:TOP_CONTEXTS]
    d = df.select(
        s("geographical_region").alias("region"),
        pl.when(s("context") == "").then(pl.lit("(none)"))
          .when(s("context").is_in(top_ctx)).then(s("context"))
          .otherwise(pl.lit("Other")).alias("context"),
        s("inclusion_route").alias("route"),
        s("ancestry_status").alias("status"),
        s("variant_id_source").alias("idsrc"),
        s("intergenic").alias("intergenic"),
        pl.when(s("ethnicity") != "").then(pl.lit("Reported"))
          .otherwise(pl.lit("Not reported")).alias("eth"),
        pl.when(s("publication_year").str.contains(r"^\d{4}$"))
          .then(s("publication_year")).otherwise(pl.lit("unknown")).alias("year"),
    )
    dims = ["region", "context", "route", "status", "idsrc", "intergenic", "eth", "year"]
    cube = d.group_by(dims).len()
    levels = {}
    for k in dims:
        vals = sorted(cube[k].unique().to_list())
        levels[k] = vals
        cube = cube.with_columns(pl.col(k).replace_strict({v: i for i, v in enumerate(vals)},
                                                          return_dtype=pl.Int32))
    # region values are lists ("West Africa, NR"): charts count each region named
    region_names = sorted({r.strip() for v in levels["region"] for r in v.split(",") if r.strip()})
    region_combos = [[region_names.index(r.strip()) for r in v.split(",") if r.strip()]
                     for v in levels["region"]]

    # ---- column completeness ----------------------------------------------
    filled = df.select([(s(c) != "").sum().alias(c) for c in NAMES]).row(0)
    rname = {r["id"]: r["name"] for r in RESOURCES}
    completeness = [{"col": c["name"], "group": c["group"], "desc": c["desc"],
                     "source": "; ".join(f"{rname[k]} ({v})" for k, v in c["sources"].items()),
                     "steps": c["steps"],
                     "filled": int(f), "pct": round(100 * f / max(n, 1), 2)}
                    for c, f in zip(COLUMNS, filled)]

    # ---- other panels -----------------------------------------------------
    phen = (df.group_by("phenotype")
              .agg(pl.len().alias("records"), pl.col("id_in_source").n_unique().alias("studies"),
                   pl.col("id").n_unique().alias("variants"))
              .sort("records", descending=True).head(25))
    sample = (df.sample(min(SAMPLE_ROWS, n), seed=7).select(NAMES)
                .with_columns([s(c) for c in NAMES]))
    qc = pl.read_csv(args.qc, separator="\t", infer_schema=False).to_dicts()
    manifest = (pl.read_csv(args.manifest, separator="\t", infer_schema=False).to_dicts()
                if Path(args.manifest).stat().st_size else [])
    qfile = out_dir / "curation_queue.summary.json"
    queue = json.loads(qfile.read_text()) if qfile.exists() else None
    settings_file = out_dir / "run_settings.yaml"
    c = meta.get("counts", {})
    sym = df.filter(s("gene_symbol_external") != "").height
    ctx = df.filter(s("context") != "").height
    exp_cols = cfg.get("export", {}).get("columns", "all")

    data = {
        "release": meta.get("release"),
        "run_date": meta.get("run_date"),
        "runtime_seconds": meta.get("runtime_seconds"),
        "sources": meta.get("sources", {}),
        "tiles": {
            "records": n,
            "variants": int(df["id"].n_unique()),
            "studies": int(df["id_in_source"].n_unique()),
            "publications": int(df["PUBMEDID"].n_unique()),
            "gene_symbol_pct": round(100 * sym / max(n, 1), 1),
            "context_pct": round(100 * ctx / max(n, 1), 1),
            "ethnicity_pct": round(100 * c.get("with_ethnicity", 0) / max(n, 1), 1),
            "rescued_records": c.get("rescued_variant_records", 0),
            "new_records": c.get("new_records", n),
        },
        "compared": meta.get("compared_with_previous_imports") or [],
        "qc": qc,
        "levels": levels,
        "region_names": region_names,
        "region_combos": region_combos,
        "cube": cube.select(dims + ["len"]).rows(),
        "dims": dims,
        "completeness": completeness,
        "phenotypes": phen.rows(),
        "sample": {"columns": NAMES, "rows": sample.rows()},
        "manifest": manifest,
        "queue": queue,
        "variant_ids": meta.get("variant_ids", {}),
        "routes": meta.get("studies_by_inclusion_route", {}),
        "export_columns": NAMES if exp_cols in (None, "all", ["all"]) else exp_cols,
        "splits_help": ALL_SPLITS_HELP,
        "new_only": sorted(NEW_ONLY),
        "settings_yaml": settings_file.read_text() if settings_file.exists() else "",
    }
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    html = TEMPLATE.read_text().replace("/*__AGMP_DATA__*/null", payload)
    Path(args.out).write_text(html)
    log.info("report: %s (%.1f MB, %d cube rows)", args.out, len(html) / 1e6, cube.height)


if __name__ == "__main__":
    main()
