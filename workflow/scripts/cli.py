"""AGMP_GWAS command line.

Runs inside agmp_gwas_container (the `agmp_gwas` launcher starts it there).
It turns simple options into a settings file and runs the Snakemake
workflow, so nobody needs to know Snakemake or edit YAML.

    agmp_gwas run                       latest GWAS Catalog release, default settings
    agmp_gwas run --previous old.csv    also mark variants already in AGMP as not new
    agmp_gwas run --columns id,PUBMEDID,phenotype,p-value
    agmp_gwas columns                   what each column holds and where it comes from
    agmp_gwas check                     check the software and internet access

The workflow code lives in the container; everything a run writes (results/,
cache/, logs/) goes to the work directory, which is the current directory.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml

WORKFLOW = Path(__file__).resolve().parent.parent   # code: Snakefile, config, scripts
WORK = Path.cwd()                                   # data: results, caches, settings
sys.path.insert(0, str(WORKFLOW / "scripts"))
from columns import COLUMNS, NAMES, RESOURCES  # noqa: E402
from export import ALL_FILES, NEW_ONLY  # noqa: E402


def csv_list(text: str) -> list[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def abspath(p: str) -> str:
    return str(Path(p).expanduser().resolve())


# --------------------------------------------------------------------- run
def build_settings(a) -> dict:
    """Only what the user changed; everything else comes from workflow/config/config.yaml."""
    s: dict = {}

    def put(path, value):
        d = s
        for k in path[:-1]:
            d = d.setdefault(k, {})
        d[path[-1]] = value

    if a.release:
        put(["catalog", "release"], a.release)
    if a.local_associations or a.local_ancestry:
        if not (a.local_associations and a.local_ancestry and a.release and a.release != "auto"):
            sys.exit("--local-associations and --local-ancestry go together, with --release DATE")
        put(["catalog", "local_associations"], abspath(a.local_associations))
        put(["catalog", "local_ancestry"], abspath(a.local_ancestry))
    if a.results_dir:
        put(["results_dir"], abspath(a.results_dir))
    # comparison with earlier data: off unless asked for
    if a.previous:
        missing = [p for p in a.previous if not Path(p).expanduser().exists()]
        if missing:
            sys.exit(f"--previous file not found: {', '.join(missing)}")
        put(["compare", "previous_imports"], [abspath(p) for p in a.previous])
    elif a.no_previous:
        put(["compare", "previous_imports"], [])
    if a.compare_release:
        p = Path(a.compare_release).expanduser()
        if not (p.is_dir() and (p / "release.json").exists()) and not (p.is_file() and p.name.endswith(".json")):
            sys.exit(f"--compare-release: {p} is not an earlier release folder (it needs a "
                     "release.json) or a release.json file")
        put(["compare", "previous_release"], abspath(a.compare_release))
    # export
    if a.columns:
        cols = NAMES if a.columns == "all" else csv_list(a.columns)
        bad = [c for c in cols if c not in NAMES]
        if bad:
            sys.exit(f"Unknown column(s): {', '.join(bad)}. See: agmp_gwas columns")
        put(["export", "columns"], "all" if cols == NAMES else cols)
    for flag, key in ((a.csv, "csv"), (a.excel, "excel")):
        if flag is None:
            continue
        if flag in ("all", "none"):
            put(["export", key], flag)
            continue
        files = csv_list(flag)
        bad = [x.lstrip("-") for x in files if x != "all" and x.lstrip("-") not in ALL_FILES]
        if bad:
            sys.exit(f"Unknown file(s) for --{key}: {', '.join(bad)}.\n"
                     f"Choose from: {', '.join(ALL_FILES)} (or all / none)")
        put(["export", key], files)
    if a.gzip:
        put(["export", "gzip"], True)
    # variant handling
    if a.no_rescue:
        put(["variants", "rescue_positional"], False)
    if a.no_flanking:
        put(["variants", "flank_genes"], False)
    if a.no_vep:
        put(["variants", "vep_context"], False)
    if a.keep_unresolved_positional:
        put(["variants", "keep_unresolved_positional"], True)
    # study selection
    if a.no_country_route:
        put(["selection", "include_african_recruitment"], False)
    if a.no_description_route:
        put(["selection", "include_african_sample_description"], False)
    # triage
    if a.no_triage:
        put(["triage", "enabled"], False)
    if a.no_fulltext:
        put(["triage", "fulltext"], False)
    return s


def cmd_run(a):
    if a.settings and not Path(a.settings).exists():
        sys.exit(f"--settings file not found: {a.settings}")
    settings = build_settings(a)
    run_dir = WORK / ".agmp_gwas"
    run_dir.mkdir(exist_ok=True)
    sfile = run_dir / "settings.yaml"
    sfile.write_text("# Written by agmp_gwas run on " + datetime.now().isoformat(timespec="seconds")
                     + "\n" + yaml.safe_dump(settings, sort_keys=False) if settings else
                     "# Written by agmp_gwas run: no changes to the default settings\n{}\n")
    cores = a.cores or int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or min(os.cpu_count() or 1, 6)
    cmd = [sys.executable, "-m", "snakemake", "--snakefile", str(WORKFLOW / "Snakefile"),
           "--directory", str(WORK), "--cores", str(cores), "--rerun-incomplete",
           "--quiet", "rules",
           "--configfile", str(WORKFLOW / "config" / "config.yaml")]
    if a.settings:
        cmd.append(abspath(a.settings))
    cmd.append(str(sfile))
    if a.dry_run:
        cmd.append("--dry-run")
    if a.rebuild:
        cmd += ["--forcerun", "build"]
    if a.unlock:
        cmd.append("--unlock")
    cmd += a.snakemake_args or []
    print("Settings changed from the defaults:",
          ("\n" + yaml.safe_dump(settings, sort_keys=False)) if settings else "none", file=sys.stderr)
    print("Work directory:", WORK, file=sys.stderr)
    code = subprocess.call(cmd, cwd=WORK)
    if code == 0 and not (a.dry_run or a.unlock):
        res = Path(settings.get("results_dir") or (WORK / "results"))
        reports = sorted(res.glob("*/report.html"), key=lambda p: p.stat().st_mtime)
        if reports:
            print(f"\nDone. Results: {reports[-1].parent}\nReport:  {reports[-1]}", file=sys.stderr)
    sys.exit(code)


# ----------------------------------------------------------------- columns
def cmd_columns(a):
    names = {r["id"]: r["name"] for r in RESOURCES}
    for c in COLUMNS:
        if a.name and c["name"] != a.name:
            continue
        print(f"\n{c['name']}  [{c['group']}]\n  {c['desc']}")
        for rid, field in c["sources"].items():
            print(f"  from: {names[rid]}: {field}")
        for i, step in enumerate(c["steps"], 1):
            print(f"  {i}. {step}")
    if not a.name:
        print("\nChoose columns with: agmp_gwas run --columns " + ",".join(NAMES[:4]) + ",...")
        print("Output files (--csv / --excel):", ", ".join(ALL_FILES))
        print("  written only with --previous:", ", ".join(sorted(NEW_ONLY)))


# ------------------------------------------------------------------- check
def cmd_check(_a):
    import importlib

    ok = True
    print(f"Python      {sys.version.split()[0]}  ({sys.executable})")
    for mod in ["polars", "pyarrow", "requests", "yaml", "snakemake", "xlsxwriter"]:
        try:
            m = importlib.import_module(mod)
            print(f"{mod:<11} {getattr(m, '__version__', 'ok')}")
        except Exception as e:  # noqa: BLE001
            ok = False
            print(f"{mod:<11} MISSING ({e})")
    try:
        import polars as pl
        pl.DataFrame({"x": [1, 2]}).select(pl.col("x").sum())
        print("polars runs on this CPU: yes")
    except Exception as e:  # noqa: BLE001
        ok = False
        print(f"polars runs on this CPU: NO ({e})")
    import requests
    for name, url in [("GWAS Catalog FTP", "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"),
                      ("Ensembl REST", "https://rest.ensembl.org/info/ping?content-type=application/json"),
                      ("UniProt", "https://rest.uniprot.org/uniprotkb/P69905.json"),
                      ("PubMed", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"),
                      ("Europe PMC", "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=x&format=json")]:
        try:
            r = requests.get(url, timeout=20)
            print(f"{name:<17} reachable (HTTP {r.status_code})")
        except Exception as e:  # noqa: BLE001
            ok = False
            print(f"{name:<17} NOT reachable ({type(e).__name__})")
    print(f"Work directory: {WORK}")
    print(f"Free disk there: {shutil.disk_usage(WORK).free / 1e9:.1f} GB")
    ver = WORKFLOW.parent / "VERSION"
    print(f"AGMP_GWAS version: {ver.read_text().strip() if ver.exists() else 'unknown'}")
    print("\nAll good." if ok else "\nSome checks failed (see above).")
    sys.exit(0 if ok else 1)


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(prog="agmp_gwas", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the pipeline")
    g = r.add_argument_group("basics")
    g.add_argument("--cores", type=int, help="CPU cores (default: up to 6, or the SLURM allocation)")
    g.add_argument("--release", help="GWAS Catalog release date, or 'auto' for the latest")
    g.add_argument("--results-dir", help="where results go (default: results/ in the work directory)")
    g.add_argument("--settings", help="a YAML file of settings to change (advanced)")
    g.add_argument("--dry-run", action="store_true", help="show what would run, run nothing")
    g.add_argument("--rebuild", action="store_true",
                   help="rebuild this release even if results exist (cached lookups are kept)")
    g.add_argument("--unlock", action="store_true", help="clear a lock left by an interrupted run")
    g = r.add_argument_group("optional: compare with earlier data (give the path to it)")
    g.add_argument("--previous", action="append", metavar="FILE",
                   help="file already imported into AGMP (CSV/parquet with an 'id' column); "
                        "repeat for several files")
    g.add_argument("--no-previous", action="store_true",
                   help="ignore any previous_imports in the settings")
    g.add_argument("--compare-release", metavar="PATH",
                   help="QC: stop if records dropped sharply against this earlier release "
                        "(its results/<date> folder or release.json)")
    g = r.add_argument_group("output files")
    g.add_argument("--columns", help="'all' or a comma-separated list (see: agmp_gwas columns)")
    g.add_argument("--csv", metavar="FILES",
                   help="files to write as CSV: all (default), none, a comma-separated list, "
                        "or all minus some (all,-combined)")
    g.add_argument("--excel", metavar="FILES",
                   help="files to write as Excel: none (default), all, or a comma-separated list")
    g.add_argument("--gzip", action="store_true", help="compress the CSV files (.csv.gz)")
    g = r.add_argument_group("what to include")
    g.add_argument("--no-rescue", action="store_true", help="do not convert positional IDs to rsIDs")
    g.add_argument("--no-flanking", action="store_true",
                   help="do not add nearest genes for rescued intergenic variants")
    g.add_argument("--no-vep", action="store_true", help="do not fill missing context with Ensembl VEP")
    g.add_argument("--keep-unresolved-positional", action="store_true",
                   help="keep positional IDs that could not be converted, as chrN:pos")
    g.add_argument("--no-country-route", action="store_true",
                   help="skip studies included only because they recruited in Africa")
    g.add_argument("--no-description-route", action="store_true",
                   help="skip studies included only because their sample description names African samples")
    g.add_argument("--no-triage", action="store_true", help="skip the curation queue")
    g.add_argument("--no-fulltext", action="store_true",
                   help="curation queue: abstracts only, no open-access full text")
    g = r.add_argument_group("offline use")
    g.add_argument("--local-associations", metavar="ZIP", help="local copy of the associations zip")
    g.add_argument("--local-ancestry", metavar="TSV", help="local copy of the ancestry file")
    r.add_argument("snakemake_args", nargs=argparse.REMAINDER,
                   help="after --, extra options passed to Snakemake (advanced)")
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("columns", help="what each output column holds and where it comes from")
    c.add_argument("name", nargs="?", help="one column, e.g. geographical_region")
    c.set_defaults(func=cmd_columns)
    sub.add_parser("check", help="check software and internet access").set_defaults(func=cmd_check)
    sub.add_parser("docs", help="regenerate the documentation's data from the code").set_defaults(
        func=lambda _a: __import__("docs_data").main())
    sub.add_parser("version", help="print the version").set_defaults(
        func=lambda _a: print((WORKFLOW.parent / "VERSION").read_text().strip()))

    a = ap.parse_args()
    if getattr(a, "snakemake_args", None) and a.snakemake_args[0] == "--":
        a.snakemake_args = a.snakemake_args[1:]
    a.func(a)


if __name__ == "__main__":
    main()
