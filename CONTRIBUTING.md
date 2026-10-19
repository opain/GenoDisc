# Contributing to the GenoDisc pipeline

This guide shows **how to add a new optional analysis** to the Snakemake pipeline
(`pipeline/`), using the GCTB **SBayesRC** analysis (SNP-h² + polygenicity) as the worked
example. The same recipe applies to any new per-GWAS analysis.

AI assistants: see [`AGENTS.md`](AGENTS.md) for the condensed checklist.

## Layout

- `pipeline/Snakefile` — includes the per-topic rule files.
- `pipeline/rules/*.smk` — rules. `dependencies.smk` stages software/reference data into
  `resdir`; `sumstat_qc.smk` defines the GWAS-list DataFrames; `report.smk` assembles the
  `package_results` target.
- `pipeline/scripts/*.R` — worker scripts (one per rule), plus
  `scripts/functions/{utils_functions.R,package_results_functions.R}` and `scripts/reader.R`.
- `pipeline/config.yaml` — the base config (every flag defaults here; user/test configs extend it).
- `pipeline/shiny_apps/output_explorer_embed/` — the results viewer (reads the packaged bundle).

## Conventions (learned from clump/cojo/finemap)

1. **Config flag.** Add `<name>: F` to `config.yaml` with a comment. Read it in rules via
   `config["<name>"] == "T"` and in R via `read_param(config, '<name>', return_obj = F) == "T"`.
2. **Per-population LD reference.** `sumstat_qc.smk` builds `gwas_list_df` and `gwas_list_df_eur`.
   Pick the reference by the GWAS `population` column with a `params` lambda interpolated into the
   shell path (copy `rule clump`/`rule cojo`), **not** finemap.R's hardcoded EUR. Subset the GWAS
   list to the populations an analysis supports (e.g. `gwas_list_df_sbrc` = EUR/EAS/AFR).
3. **Reference staging.** Add a `download_*`/`install_*` rule in `dependencies.smk` (model binaries
   on `rule download_gcta`; large tarballs on `rule download_ldsc_scores`; GitHub R packages on the
   SHA-pinned `install_genoutils` pattern). Append its output to `resource_inputs` under
   `if config["<name>"] == "T":`. Share large references with sibling pipelines (e.g. GenoPred)
   rather than re-downloading where possible.
4. **Rule + worker script.** New `rules/<name>.smk`, `include`d in `Snakefile`. The rule invokes
   `scripts/<name>.R --pipeline_dir {workflow.basedir} --gwas {wildcards.gwas}
   --config_file {params.config_file} ...`; declare `conda:`, `benchmark:`, `log:`, and realistic
   `resources:`. The worker sources `scripts/functions/utils_functions.R` for `read_param`.
5. **Report target.** In `report.smk`, under `if config["<name>"] == "T":`, append the final
   per-GWAS output to `myoutput` (iterate the right GWAS subset).
6. **Packaging.** Add a `process_<name>(config, gwas)` to `package_results_functions.R` (gated on
   the flag, returns `NULL` when absent) and call it in `package_results.R`. Nesting the result
   under an existing block (e.g. `gwas_qc$<name>_dat`) needs **no** change to `.gd_block_ids` in
   `scripts/reader.R`; a new top-level block does.
7. **Viewer (optional).** Add the flag to `parse_config_flags()` in
   `shiny_apps/output_explorer_embed/R/utils_config.R`, widen the relevant `toggle(...)` gate in
   `app.R`, and render from the packaged block in the appropriate `R/mod_*.R` module.
8. **Tests.** Set `<name>: T` in `tests/config_test.yaml` and `prepare_resources.yaml`.

## Worked example — GCTB SBayesRC (SNP-h² + polygenicity)

- **Tool:** official **gctb 2.5.5** static Linux binary from the GCTB site
  (`https://cnsgenomics.com/software/gctb/download/gctb_2.5.5_Linux.zip`; the host redirects to
  `gctbhub.cloud.edu.au`). Staged by `rule download_gctb` → `{resdir}/software/gctb/gctb_2.5.5_Linux/gctb`.
- **Reference data (shared with GenoPred):** per-population eigen-decomposed LD at
  `…/GenoPred/pipeline_resources/data/sbayesrc_ref/{EUR,EAS,AFR}/` (`ldm.info`, `snp.info`,
  `block*.eigen.bin`) and the baseline annotation `…/sbayesrc_annot/annot_baseline2.2.txt`
  (header `SNP Intercept …`). Defaulted in `rules/sbayesrc.smk` when `sbayesrc_ldref`/
  `sbayesrc_annot` are `NA`; override in config to point elsewhere.
- **Recipe** (`scripts/sbayesrc.R`): build a COJO `.ma` (`SNP A1 A2 freq b se p N`) from
  `{gwas}.cleaned.gz`, then
  1. `gctb --ldm-eigen <ldref> --gwas-summary <gwas>.ma --impute-summary --out <prefix> --thread N`
     → `<prefix>.imputed.ma`
  2. `gctb --ldm-eigen <ldref> --gwas-summary <prefix>.imputed.ma --sbayes RC --annot <annot> --out <prefix> --thread N`
     → `<prefix>.parRes`, `<prefix>.snpRes`
- **Outputs kept:** parse `<prefix>.parRes` (header `Parameter Mean SD`) — `hsq` row = SNP-h²
  (Mean, SD); polygenicity = non-zero-effect variants = `sum(NumSnp2..NumSnp5)` (NumSnp1 is the
  null component); `pi = that / sum(NumSnp1..5)`. There is **no** `Pi`/`NnzSnp` row. Written to
  `{gwas}.GW.sbayesrc.clean.csv` (`gwas,h2_est,h2_se,n_nonzero,n_snp_total,pi_est`); the per-SNP
  `.snpRes` is discarded. Packaged into `gwas_qc$sbayesrc_dat` and shown in the viewer's
  "SNP-h² & rG" tab as a "SBayesRC (h² & polygenicity)" sub-tab.
- **Note:** SBayesRC is genome-wide; the chr22 test config is a mechanics smoke test only
  (h²≈0). Use a whole-genome GWAS for meaningful estimates.

## Running tests (must run on the HPC)

The pipeline runs inside the `genodisc` conda env on the CREATE HPC (the dev VM can't execute it —
the env binaries resolve `/cephfs` paths and the Snakefile guards that the env is active). Submit
via SSH exactly as the web app does (`modules/hpc_interface.py`): SSH to `hpc.create.kcl.ac.uk` and
`sbatch` a wrapper that `source`s `bin/activate_genodisc.sh` and runs, from `pipeline/`:

```bash
snakemake <target> --configfile tests/config_test.yaml --use-conda -c<N>
# e.g. the SBayesRC target (isolated config):
snakemake tests/output/results/COAD01chr22/sbayesrc/COAD01chr22.GW.sbayesrc.clean.csv \
  --configfile tests/config_test_sbayesrc.yaml --use-conda -c4
```

To exercise a single leaf rule against already-built upstream outputs, add
`--allowed-rules <rule> <its download rules> --rerun-triggers mtime`.
