# Contributing to GenoDisc

This guide shows how to add a new optional **analysis ("method")** or **reference dataset** to the
GenoDisc Snakemake pipeline. GCTB **SBayesRC** (SNP-h² + polygenicity) is the worked example; its
completed contribution package is `pipeline/docs/methods/sbayesrc.md`.

A contributor implements the analysis/dataset in the pipeline and fills in a structured
**contribution form** (`pipeline/docs/method_contribution_template.md` or
`dataset_contribution_template.md`). A maintainer then integrates it into the GenoDisc website
(the job-submission panel) and the user documentation — **contributors don't need to do that part**
(see "Maintainer integration" at the end).

AI assistants: see [`AGENTS.md`](AGENTS.md) for the condensed checklist.

## Layout
- `pipeline/Snakefile` — includes the per-topic rule files.
- `pipeline/rules/*.smk` — rules. `dependencies.smk` stages software/reference data into `resdir`;
  `sumstat_qc.smk` defines the GWAS-list DataFrames; `report.smk` assembles the `package_results`
  target.
- `pipeline/scripts/*.R` — worker scripts (one per rule), plus
  `scripts/functions/{utils_functions.R,package_results_functions.R}` and `scripts/reader.R`.
- `pipeline/config.yaml` — the base config (every flag defaults here; user/test configs extend it).
- `pipeline/shiny_apps/output_explorer_embed/` — the results viewer (reads the packaged bundle).

## Adding an analysis — conventions (learned from clump/cojo/finemap)
1. **Config flag.** Add `<name>: F` to `config.yaml` with a comment. Read it in rules via
   `config["<name>"] == "T"` and in R via `read_param(config, '<name>', return_obj = F) == "T"`.
2. **Per-population reference.** `sumstat_qc.smk` builds `gwas_list_df` / `gwas_list_df_eur`. Pick the
   reference by the GWAS `population` column via a `params` lambda interpolated into the shell path
   (copy `rule clump`/`rule cojo`). Subset the GWAS list to the populations the method supports (e.g.
   `gwas_list_df_sbrc` = EUR/EAS/AFR) and iterate that subset in `report.smk`.
3. **Reference staging.** Add a `download_*`/`install_*` rule in `dependencies.smk` (see "Adding a
   reference dataset" below). Append its output to `resource_inputs` under `if config["<name>"]=="T":`.
4. **Rule + worker script.** New `rules/<name>.smk`, `include`d in `Snakefile`; invokes
   `scripts/<name>.R --pipeline_dir {workflow.basedir} --gwas {wildcards.gwas} --config_file
   {params.config_file} ...` with `conda:`, `benchmark:`, `log:`, realistic `resources:` (declare
   `threads` from a config key when the tool is multi-threaded). The worker sources
   `scripts/functions/utils_functions.R` for `read_param`.
5. **Report target.** In `report.smk`, under `if config["<name>"] == "T":`, append the final per-GWAS
   output to `myoutput` (iterate the right GWAS subset).
6. **Packaging.** Add `process_<name>(config, gwas)` to `package_results_functions.R` (gated on the
   flag, returns `NULL` when absent) and call it in `package_results.R`. Nesting under an existing
   block (e.g. `gwas_qc$<name>_dat`) needs **no** `.gd_block_ids` change in `scripts/reader.R`; a new
   top-level block does.
7. **Viewer (optional).** Add the flag to `parse_config_flags()` in
   `shiny_apps/output_explorer_embed/R/utils_config.R`, widen the relevant `toggle(...)` gate in
   `app.R`, and render from the packaged block in the appropriate `R/mod_*.R` module.
8. **Tests.** Add a dedicated test config (e.g. `tests/config_<name>.yaml`) and enable the flag in
   `prepare_resources.yaml`.

Don't rely on training data for external-tool versions/flags — check the official docs and confirm
against the real binary (e.g. GCTB's `.parRes`/MCMC output differs from older docs).

## Adding a reference dataset
Reference data is staged into the pipeline's resource directory (`resdir`) by a download rule, so a
run is self-contained. Follow the pattern used by `download_ldsc_scores`/`download_ldsc_hm3`,
`download_eqtlgen`, and `download_biomart` in `rules/dependencies.smk`:
1. **Download/prep rule** `download_<name>`: `params: resdir=resdir`, `conda: "../envs/main.yaml"`,
   `benchmark:`/`log:`. Shell either `wget -O {params.resdir}/data/<name>/<file> <official URL>` then
   `unzip`/`tar -xf`/`gzip -d` and `rm` the archive, or `Rscript .../scripts/<name>.R --resdir
   {params.resdir}`. Prefer a stable **official** URL over mirrors/drives.
2. **Output under `{resdir}/data/<name>/…`** — a concrete sentinel file (e.g. `.../ldm.info`), or
   `touch(f"{resdir}/data/<name>.done")` when the payload is many files.
3. **`resource_inputs` entry** appended in the matching gated block (or the always-on core `extend`).
   `rule prepare_resources` consumes `resource_inputs`, so this is what makes it stage.
4. **Config key** (optional): a T/F feature toggle, or an explicit path key with an `NA` fallback to
   the staged copy in `resdir`.
5. **Encode version/build/ancestry** in the filename/subdir (e.g. `grch37`, `<POP>/`), the convention
   used throughout `dependencies.smk`.

Fill in `pipeline/docs/dataset_contribution_template.md` (source/URL, version, build, ancestry,
licence/redistribution terms, citation, download rule, config key, size).

## Fill in the contribution form
Copy `pipeline/docs/method_contribution_template.md` (or `dataset_contribution_template.md`) to
`pipeline/docs/methods/<flag>.md` and complete every field — including licence, citation, maintainer
contact, and known failure modes. The viewer-presentation prose and the technical write-up are
drop-in text the maintainer pastes into the docs. `pipeline/docs/methods/sbayesrc.md` is a completed
example. (`pipeline/docs/*` is intentionally separate from the repo-root `docs/` RMarkdown site.)

## Worked example — GCTB SBayesRC
See `pipeline/docs/methods/sbayesrc.md` for the full specifics. In brief: official **gctb 2.5.5**
binary (`rule download_gctb`); eigen-LD + annotation references downloaded from the official GCTB
site by `rule download_sbayesrc_ref`/`download_sbayesrc_annot`; the worker (`scripts/sbayesrc.R`)
builds a COJO `.ma`, runs `gctb --impute-summary` then `gctb --sbayes RC --out-freq 1`, parses the
post-burn-in `hsq` MCMC trace (posterior mode + 95% HPD + mean/SD) and `.parRes` (polygenicity),
writes `{gwas}.GW.sbayesrc.clean.csv` (+ `h2_samples.txt`), packaged at `gwas_qc$sbayesrc_dat`.
SBayesRC is genome-wide — the chr22 test config is a mechanics smoke test (h²≈0); use a whole-genome
GWAS for meaningful estimates.

## Worked example — adding a FUSION TWAS weight panel
Adding a new FUSION expression-weight panel is **not** the generic "Adding a reference dataset" flow
above — TWAS panels have their own machinery. See `pipeline/docs/methods/fetal_brain_twas.md` and
`herv_rtwas.md` for two completed examples (fetal-brain gene-level; HERV retrotranscriptome). Steps:

1. **Selection is config-driven.** A panel is a string the run assembles into the `weights` list in
   `rules/dependencies.smk` (`gtex_weights` + `non_gtex_weights`, plus the `twas_panel_*` appends).
   A GTEx/Broad panel is just an id appended to `gtex_weights`/`non_gtex_weights`; a bespoke panel
   gets its own `twas_panel_<name>: F` flag that appends its id.
2. **Acquisition.** `download_non_gtex_weights` hardcodes the Broad mirror and **pins no checksum** —
   do not copy that for a contributed panel. Add a dedicated rule that pins the official URL **and a
   checksum** (`curl … | md5sum -c -`), extracts, and runs a prep script (see `prep_fetal_brain.R` /
   `prep_herv.R`). Weights land in `{resdir}/data/fusion_snp_weights/<PANEL>/` with a `<PANEL>.pos`
   (`PANEL WGT ID CHR P0 P1 N`), WGT relative to the panel dir (`--weights_dir`).
3. **Sample size N.** `scripts/insert_n_nongtex.R` holds a **hardcoded per-panel N table**; a panel
   missing from it silently gets `N=NA`. A dedicated prep script should set N explicitly.
4. **Feature IDs.** `combine_twas.R` merges the FUSION `ID` to a gene symbol via biomaRt assuming an
   Ensembl id. Gene panels must therefore use Ensembl gene IDs (unversioned) or HGNC symbols (which
   `insert_n_nongtex.R` remaps). Non-gene features (e.g. HERVs) must keep their own id and be routed
   around the gene-based steps (see point 6).
5. **LD reference.** `_twas_ref_prefix` (`rules/twas_pwas.smk`) sends every panel to the
   HapMap3-restricted reference **except** those listed as dense (`psychencode`, `CMC.BRAIN.HERV`).
   If a panel's weights were trained on a full (non-HapMap3) reference, add it to that list or ~80%
   of its weight SNPs are dropped. Verify SNP overlap against the reference before deciding.
6. **Downstream + multiple testing.** TWAS-GSEA/DrugTargetor/CMAP key on gene symbol and are driven
   by `weights_nosplice`; the gene FDR in `combine_twas.R` / `package_results_functions.R` is pooled
   across panels. A non-gene panel must be removed from `weights_nosplice`, excluded from the gene
   FDR pool (given a separate FDR), and guarded at parse time (error if a gene-based analysis is
   requested with no gene-producing panel).
7. **Display + gating.** Add the panel id → label to `tidy_panel_names()` and enumerate it in
   `read_fusion_exp` (`$panels`); extend the "TWAS active" gates (`report.smk`, `dependencies.smk`
   `resource_inputs`, the viewer `utils_config.R`, `create_report.Rmd`) so the panel's output and
   viewer tab appear when it is the only TWAS panel enabled.
8. **Bump `VERSION`** (minor for a new panel; the bundle `schema_version` is unchanged — TWAS nests
   under the existing `mol_assoc/exp/fusion` block).

(Gaps this example surfaced, now partly addressed here: the checksum-vs-practice inconsistency, the
hardcoded N table, the implicit HapMap3 LD routing, the cross-panel FDR, and the `weights_nosplice`
splicing special-case were all previously undocumented.)

## Testing
The pipeline runs inside its conda environment (`pipeline/envs/main.yaml`, plus per-rule envs
materialised by `--use-conda`). From `pipeline/`:
```bash
# 1. Stage the references a config needs (downloads into resdir):
snakemake --configfile prepare_resources.yaml prepare_resources --use-conda -c<N>
# 2. Run a test target (locally or on your own cluster):
snakemake tests/output/results/COAD01chr22/sbayesrc/COAD01chr22.GW.sbayesrc.clean.csv \
  --configfile tests/config_test_sbayesrc.yaml --use-conda -c<N>
```
Use enough cores for multi-threaded rules — Snakemake caps a rule's `threads` at the `-c` value, so
e.g. a 10-thread rule needs `-c10` or higher. To exercise one leaf rule against already-built
upstream outputs: `--allowed-rules <rule> <its download rules> --rerun-triggers mtime`.

## Maintainer integration
Surfacing a contributed method/dataset on the GenoDisc website (the job-submission analysis panel)
and in the user documentation is done by the **maintainer**, using the contributor's filled-in form.
Contributors do not need to touch the web app or the documentation site.
