# output_explorer_embed memory profiling

Measures peak R-process RSS + `/tmp` footprint of the `output_explorer_embed`
Shiny app under multi-GWAS `.tar.gz` bundles of varying size, so we can pick an
appropriate shinyapps.io plan tier or bundle-size cap.

**Approach.** Synthesises N-GWAS bundles by duplicating the shipped ALS example
(`data/als_bundle.tar.gz`) N times, uploads each to a locally-running instance
of the app via headless Firefox, walks every tab, and samples the R process's
`VmRSS` from `/proc/<pid>/status` at 2 Hz.

## Prerequisites

- R with `shiny`, `data.table`, `jsonlite`, `ggplot2`, `DT`, `ggplot2`, and the
  full list of app deps (see the app's `renv.lock` / package footer).
- Python with `selenium`.
- Firefox + geckodriver (matched versions).
- `curl`, `tar`, `awk`, `du` (standard).

Point env vars at non-`PATH` binaries as needed:

```
export RSCRIPT=/path/to/Rscript
export PYTHON=/path/to/python
export FIREFOX_BIN=/path/to/firefox
export GECKODRIVER_BIN=/path/to/geckodriver
# Optional, for conda-based Firefox that needs the env's lib dir on LD_LIBRARY_PATH:
export PROFILE_LD_LIBRARY_PATH=/path/to/env/lib
```

## Run

```
bash run_all.sh
```

Defaults to `N ∈ {1, 2, 4, 8, 16}`. Override with `N_LIST="1 4 16"` env.
Total wall time: ~10 minutes on modest hardware.

Bundles are cached in `bundles/` — re-runs skip already-built bundles.

Results land in `results/N<NN>/` per run (samples CSV, checkpoints JSONL,
per-tab screenshots, driver + app logs) plus `results/summary.md` +
`results/peak_rss_vs_N.png` from the aggregator.

## Files

- `scale_bundle.R` — synthesize an N-GWAS `.tar.gz` from a single-GWAS source.
- `profile_sampler.sh` — polls `/proc/<pid>/status` VmRSS + a scoped TMPDIR footprint at 2 Hz.
- `profile_run.py` — Selenium driver: upload bundle, walk every tab and sub-tab, screenshot each.
- `run_one.sh` — orchestrator for one N (starts app, sampler, driver; tears down).
- `run_all.sh` — scales bundles + runs the whole series + aggregates.
- `summary.R` — aggregates per-N CSVs into `summary.md` + `peak_rss_vs_N.png`.

## Interpreting results

The summary reports peak RSS at four moments:

1. **baseline** — RSS in the first 5 s (Shiny + libs, no bundle loaded).
2. **after upload** — RSS right after `data-app-ready="1"` flips (bundle opened,
   `comparison_long()` reactive complete). This is the *unavoidable* cost of
   loading the bundle.
3. **after tab walk** — RSS after every tab + sub-tab has rendered at least
   once. Adds DT tables, ggplot render caches, per-GWAS `mol_assoc/*`, `tx/*`,
   `snp_assoc` block reads via `gd_read()`.
4. **peak** — the max RSS seen at any point during the run.

The linear fit of peak-RSS-vs-N is projected onto shinyapps.io plan sizes
(1 GB Basic, 3 GB Standard, 8 GB Standard-L / Pro), leaving ~100 MB headroom.

## Caveats

- Synthetic bundles duplicate the same source GWAS, so per-GWAS block sizes
  are identical. Real bundles with heavier per-GWAS `snp_assoc` will scale
  worse — treat the reported "safe N" as a lower bound.
- The ALS example bundle has no cross-GWAS LDSC / rG matrix; real bundles do,
  and that adds O(N²) memory that this harness does not measure.
- Measurements are Linux-specific (`/proc/<pid>/status`).
