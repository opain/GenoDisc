# Internal results review

Scripts for producing an internal scientific review of the GenoDisc validation
results from `results_package.rds`: a per-trait summary of upstream signal
(SNP-level loci, fine-mapping, TWAS, PWAS, SMR, MAGMA gene + tissue) and the
downstream drug-repurposing analyses (TWAS-GSEA directional and
non-directional, CMAP, MAGMA DrugTargetor, GCSC) at drug and ATC level, plus
cross-trait overlap analyses.

Outputs (tables, figures, markdown report) are written to a configurable
location outside the repo so the repo stays git-clean. The default location
is `/work/results/internal_results_review/`.

## Folder layout (in this repo)

```
internal_results_review/
├── README.md                               (this file)
├── scripts/
│   ├── _paths.R                            shared path resolution (sourced by every script)
│   ├── inspect_results_package.R           top-level RDS structure dump
│   ├── extract_twas_gsea_results.R         per-trait TWAS-GSEA, MAGMA, GCSC, CMAP tables
│   ├── extract_upstream_summaries.R        per-trait loci/gene/tissue counts and gene lists
│   ├── make_figures.R                      drug-repurposing figures
│   ├── make_upstream_figures.R             upstream-signal figures
│   ├── qc_checks.R                         direction-invariant + panel-concordance + bias checks
│   └── explore_top_results.R               quick per-trait top-hits printer (interactive)
├── tests/
│   └── test_direction_sign_conventions.R   unit test for the TWAS-GSEA sign conventions (no data needed)
└── config_template/
    ├── example_config.yaml                 example packaging config used for the May 2026 review run
    └── example_gwas_list.txt               accompanying gwas_list with all 9 traits
```

## Outputs (written by the scripts)

By default the scripts write to `/work/results/internal_results_review/` with:

```
/work/results/internal_results_review/
├── tables/      (~30 MB; ~35 TSVs of per-trait results, cross-trait overlap, gene x trait matrices)
├── figures/     (~5 MB; PNGs)
└── report/      (markdown reports: genodisc_twas_gsea_internal_report.md, qc_checks.md, etc.)
```

The CMAP per-signature TSV (`cmap_drug_all_traits.tsv`) is ~260 MB and is
*not* shipped with the report by default — it is regenerable by re-running
`extract_twas_gsea_results.R` if needed.

## Configuration

The scripts read two environment variables at runtime, with sensible defaults:

| Variable | Default | What it is |
|---|---|---|
| `GENODISC_RDS` | `/work/results/results_package.rds` | Input RDS produced by `scripts/package_results.R`. |
| `GENODISC_OUTDIR` | `/work/results/internal_results_review` | Output base directory; the scripts create `tables/`, `figures/`, `report/` inside this. |
| `GENODISC_VERBOSE` | `0` | Set to `1` to log resolved paths at startup. |

Example: run against a different RDS and write to a sandbox:

```bash
GENODISC_RDS=/path/to/alternate.rds \
GENODISC_OUTDIR=/tmp/sandbox_run \
Rscript scripts/extract_twas_gsea_results.R
```

The path resolution is implemented in `scripts/_paths.R`. Each script sources
it via a small self-locating idiom that works both when run via `Rscript` and
when sourced interactively.

## How to run

In order:

```bash
cd /path/to/GenoDisc/pipeline/misc/internal_results_review

Rscript scripts/inspect_results_package.R
Rscript scripts/extract_twas_gsea_results.R
Rscript scripts/extract_upstream_summaries.R
Rscript scripts/qc_checks.R
Rscript scripts/make_figures.R
Rscript scripts/make_upstream_figures.R
```

`explore_top_results.R` is an interactive printer; run it via `Rscript` to
dump per-trait top hits to stdout, or `source()` it in an R session.

## R package dependencies

`data.table`, `ggplot2`, `stringr`, `scales`, `gridExtra` (optional for
multi-panel figures), `yaml`, `optparse` (only used by the sign-convention
unit test). Install with:

```r
install.packages(c("data.table","ggplot2","stringr","scales","gridExtra","yaml","optparse"))
```

## Verification

Run the sign-convention unit test (no data required, ~1 sec):

```bash
Rscript tests/test_direction_sign_conventions.R
```

It checks the Wilcoxon HL convention for both the formula-form (DrugTargetor
ATC) and two-vector-form (CMAP MOA) tests, and confirms that the
`Direction` / `Reversal_Z` recipes used by the pipeline produce the expected
values across drug-level, ATC-level, MOA-level and non-directional cases.

## How the columns are interpreted

Every table produced by these scripts has a consistent direction-of-effect
convention:

- `Reversal_Z > 0` always means the candidate-therapeutic direction (the
  drug or class **opposes** the disease TWAS signature).
- `Direction` is a categorical column with values `"Opposes disease"`,
  `"Matches disease"`, or `NA` (non-directional / one-sided tests).
- The raw `Estimate` column has different sign meanings across tables (drug
  vs ATC vs MOA Wilcoxon) — see `pipeline/scripts/functions/package_results_functions.R`
  and `pipeline/scripts/format_twas_gsea_*.R` for the underlying conventions
  and the corresponding fallback recipes in `read_twas_gsea_*`.

## Provenance

These scripts were developed for the May 2026 GenoDisc validation review
covering nine neuropsychiatric and neurodegenerative traits (ADHD, ALS, ALZ,
ASD, BIP, MDD, MIG, PRK, SCZ). The example config in `config_template/` is
the one that was used for that re-run.

Threshold conventions used throughout (matching the pipeline's reading
functions):

- Genome-wide significance: **P < 5×10⁻⁸** (applied to clumped and COJO
  tables; note the pipeline's `.GW.clump.clean.csv` and `.GW.cojo.clean.csv`
  files actually retain signals down to P<10⁻⁵)
- TWAS / PWAS high-confidence: FDR<0.05 AND `COLOC_logical = TRUE`
- SMR high-confidence: FDR<0.05 AND `p_HEIDI > 0.01`
- MAGMA, tissue MAGMA: FDR<0.05; tissue "retained" = retained after the
  conditional analysis the pipeline applies on top of FDR.
