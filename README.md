# AGMP_GWAS

### [Open the documentation →](https://nhlamulokhoza.github.io/AGMP_GWAS/)

GWAS Catalog curation for the African Genomic Medicine Portal (AGMP).

AGMP_GWAS reads each GWAS Catalog release, keeps every association relevant to
African populations and the diaspora, resolves and annotates the variants,
checks the result, and writes the files AGMP imports. It also writes an
interactive report and a ranked queue of studies for manual curation. The whole
workflow runs inside one container, `agmp_gwas_container`.

**Documentation:** https://nhlamulokhoza.github.io/AGMP_GWAS/ (also in `docs/`). It covers installation,
running on an HPC, a command builder with check boxes, every output column (with
the resources it comes from and each transformation applied), and the methods.

## Quick start

```bash
git clone https://github.com/nhlamulokhoza/AGMP_GWAS.git
cd AGMP_GWAS
./agmp_gwas setup      # once: builds agmp_gwas_container (about 5 minutes)
./agmp_gwas check      # confirms the software and internet access work
./agmp_gwas run        # the latest GWAS Catalog release
```

Results go to `results/<release>/`; open `report.html` there. On a SLURM
cluster, submit a job instead of running on the login node:

```bash
./agmp_gwas run --slurm
```

Comparing with earlier data is off unless you ask for it:
`--previous FILE` marks variants already imported into AGMP as not new, and
`--compare-release results/<date>` checks the record count against an earlier
release. Both take the path to the earlier data.
Each output file can be written as CSV, Excel or both:
`--excel ethnicity_report` adds an Excel copy of one file, `--csv all,-combined`
leaves one out. `./agmp_gwas run --help` lists every option, and
`./agmp_gwas columns` describes every column.

## Layout

```
agmp_gwas                    launcher: the only command you use
container/                   recipe for agmp_gwas_container, exact package versions
workflow/                    Snakefile, scripts, default settings (copied into the container)
docs/                        documentation
results/ cache/ logs/        made by runs (not part of the code)
```

The container and the launcher are a complete install: to run elsewhere, copy
`agmp_gwas` and `agmp_gwas_container.sif`.

## Output

By default the export files have the same 30 columns, in the same order, as the
June 2026 AGMP import. `workflow/scripts/columns.py` describes each one: its
sources and every transformation. The documentation, the report and
`./agmp_gwas columns` are all generated from it.
