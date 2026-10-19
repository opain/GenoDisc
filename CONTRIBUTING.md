# Contributing to GenoDisc

This is the **end-to-end** guide for adding a new optional analysis ("method") to GenoDisc — from
the Snakemake pipeline through to the method appearing in the website's job-submission panel and in
the user documentation. GCTB **SBayesRC** (SNP-h² + polygenicity) is the worked example; its
completed contribution package is `pipeline/docs/methods/sbayesrc.md`.

There are two roles:

- **Contributor** — implements the analysis in the pipeline (Part A) and fills in a structured
  **method contribution form** (Part B). A contributor does **not** edit the Django web app or the
  docs site.
- **Maintainer** — integrates a contributed method into the website submission panel and the docs
  pages using the contributor's form (Part C).

AI assistants: see [`AGENTS.md`](AGENTS.md) for the condensed checklist.

---

# Part A — Contributor: implement the pipeline analysis

## Layout
- `pipeline/Snakefile` — includes the per-topic rule files.
- `pipeline/rules/*.smk` — rules. `dependencies.smk` stages software/reference data into `resdir`;
  `sumstat_qc.smk` defines the GWAS-list DataFrames; `report.smk` assembles the `package_results`
  target.
- `pipeline/scripts/*.R` — worker scripts (one per rule), plus
  `scripts/functions/{utils_functions.R,package_results_functions.R}` and `scripts/reader.R`.
- `pipeline/config.yaml` — the base config (every flag defaults here; user/test configs extend it).
- `pipeline/shiny_apps/output_explorer_embed/` — the results viewer (reads the packaged bundle).

## Conventions (learned from clump/cojo/finemap)
1. **Config flag.** Add `<name>: F` to `config.yaml` with a comment. Read it in rules via
   `config["<name>"] == "T"` and in R via `read_param(config, '<name>', return_obj = F) == "T"`.
2. **Per-population reference.** `sumstat_qc.smk` builds `gwas_list_df` / `gwas_list_df_eur`. Pick
   the reference by the GWAS `population` column via a `params` lambda interpolated into the shell
   path (copy `rule clump`/`rule cojo`). Subset the GWAS list to the populations the method supports
   (e.g. `gwas_list_df_sbrc` = EUR/EAS/AFR) and iterate that subset in `report.smk`.
3. **Reference staging.** Add a `download_*`/`install_*` rule in `dependencies.smk` (binaries →
   model on `rule download_gcta`; tarballs → `rule download_ldsc_scores`; GitHub R packages →
   SHA-pinned `install_genoutils`). Append its output to `resource_inputs` under
   `if config["<name>"] == "T":`. **Share** large references with sibling pipelines (e.g. GenoPred)
   rather than re-downloading.
4. **Rule + worker script.** New `rules/<name>.smk`, `include`d in `Snakefile`; invokes
   `scripts/<name>.R --pipeline_dir {workflow.basedir} --gwas {wildcards.gwas} --config_file
   {params.config_file} ...` with `conda:`, `benchmark:`, `log:`, realistic `resources:`. The worker
   sources `scripts/functions/utils_functions.R` for `read_param`.
5. **Report target.** In `report.smk`, under `if config["<name>"] == "T":`, append the final
   per-GWAS output to `myoutput` (iterate the right GWAS subset).
6. **Packaging.** Add `process_<name>(config, gwas)` to `package_results_functions.R` (gated on the
   flag, returns `NULL` when absent) and call it in `package_results.R`. Nesting under an existing
   block (e.g. `gwas_qc$<name>_dat`) needs **no** `.gd_block_ids` change in `scripts/reader.R`; a new
   top-level block does.
7. **Viewer (optional).** Add the flag to `parse_config_flags()` in
   `shiny_apps/output_explorer_embed/R/utils_config.R`, widen the relevant `toggle(...)` gate in
   `app.R`, and render from the packaged block in the appropriate `R/mod_*.R` module.
8. **Tests.** Set `<name>: T` in `tests/config_test.yaml` and `prepare_resources.yaml`.

Don't rely on training data for external-tool versions/flags — check the official docs and confirm
against the real binary (e.g. GCTB's `.parRes`/MCMC output differs from older docs).

## Worked example — GCTB SBayesRC
See `pipeline/docs/methods/sbayesrc.md` for the full specifics (tool version, two-step gctb recipe,
MCMC-trace parsing for the SNP-h² posterior mode/HPD, outputs, packaging, viewer sub-tab). In brief:
official **gctb 2.5.5** binary (`rule download_gctb`); references shared with GenoPred; the worker
(`scripts/sbayesrc.R`) builds a COJO `.ma`, runs `gctb --impute-summary` then `gctb --sbayes RC
--out-freq 1`, parses the post-burn-in `hsq` MCMC trace (mode + 95% HPD + mean/SD) and `.parRes`
(polygenicity), writes `{gwas}.GW.sbayesrc.clean.csv` (+ `h2_samples.txt`), packaged at
`gwas_qc$sbayesrc_dat`. SBayesRC is genome-wide — the chr22 test config is a mechanics smoke test
(h²≈0); use a whole-genome GWAS for meaningful estimates.

## Running tests (must run on the HPC)
The pipeline runs inside the `genodisc` conda env on the CREATE HPC (the dev VM can't execute it —
env binaries resolve `/cephfs` paths and the Snakefile guards that the env is active). Submit via
SSH as the web app does (`modules/hpc_interface.py`): SSH to `hpc.create.kcl.ac.uk` and `sbatch` a
wrapper that `source`s `bin/activate_genodisc.sh` and runs, from `pipeline/`:
```bash
snakemake <target> --configfile tests/config_test.yaml --use-conda -c<N>
# isolated single-method config, e.g.:
snakemake tests/output/results/COAD01chr22/sbayesrc/COAD01chr22.GW.sbayesrc.clean.csv \
  --configfile tests/config_test_sbayesrc.yaml --use-conda -c4
```
To exercise one leaf rule against already-built upstream outputs:
`--allowed-rules <rule> <its download rules> --rerun-triggers mtime`.

---

# Part B — Contributor: fill the method contribution form

The pipeline implementation is not user-visible until the method is added to the website submission
panel and documented. The contributor provides that information as a **form** (not by editing the
web app); the maintainer integrates it.

1. Copy `pipeline/docs/method_contribution_template.md` to `pipeline/docs/methods/<flag>.md`.
2. Complete every field (identity, config params, submission-panel metadata, compute cost, outputs,
   analyses-catalogue entry, viewer presentation prose, technical write-up, verification). No TODOs.
3. The viewer-presentation prose and the technical write-up are drop-in text the maintainer pastes
   into the docs pages — write them user-facing.

`pipeline/docs/methods/sbayesrc.md` is the completed example. (`pipeline/docs/*` is intentionally
separate from the repo-root `docs/` RMarkdown site, so these files are not knitted into it.)

---

# Part C — Maintainer: integrate a contributed method

Driven entirely by the contributor's `pipeline/docs/methods/<flag>.md`. Three areas:

## C1. Website job-submission analysis panel (Django app: `/home/ollie/genodisc-devserver`)
The panel shows a toggle for a pipeline version only if that version's `config.yaml` contains the
toggle's config key (`modules/hpc_interface.py` `_version_capabilities`; `dev` → `repo/current/
pipeline/config.yaml`). Using the form's **submission-panel metadata** + **compute cost**:
- `hpc_app/forms.py` — add `<flag> = forms.BooleanField(required=False)` to
  `AnalysisSelectionFieldsMixin` (no `clean()` change; it is data-driven).
- `modules/analysis_dependencies.py` — add `ANALYSIS_LABELS["<flag>"]`, `ANALYSIS_CONFIG_KEYS["<flag>"]
  = "<config key>"`, and `ANALYSIS_DEPENDENCIES` only if the method has a parent toggle.
- `templates/_job_form.html` — add an analysis card in the form's stated section, carrying
  `data-version-key="<config key>"` so it is version-gated (clump/cojo/finemap cards have none and
  are always shown — a new method MUST set it). Add `'<flag>'` to the `BOOLEAN_FIELDS` JS array.
- `modules/job_service.py` — `build_analysis_selection` (form→JSON), `_build_config_yaml`
  (write `<flag>: T/F`), `flatten_analysis_selection` (required for modify-prefill + the estimate),
  and optionally `summarize_analysis_selection`.
- `modules/resource_estimates.py` — add a `RULE_COST['<rule>']` entry (from the form's runtime/mem/
  per_chr) and a `_rule_load` dispatch (`if flat.get('<flag>'): add('<rule>')`).

## C2. Website docs (Django hardcoded-HTML templates)
- `templates/docs_webapp.html` — add a `<tr>` to the analyses-catalogue table from the form's
  **analyses-catalogue entry** (and a dependency bullet if any).
- `templates/docs_results.html` — if the method has viewer output, add an `<h3 id="<flag>">` section
  inside the relevant `<h2>` viewer-tab block from the form's **viewer presentation** (What it shows
  / How to read it / Caveats), plus a sidebar `<a class="sub">` and a glossary `<tr>`. The figure is
  a screenshot produced by `scripts/docs/capture_viewer_screenshots.py` against the released viewer
  (add a `screenshots.json` entry); defer the `<img>` until the viewer change is released, or
  `scripts/docs/verify_results_page.py` will fail on the missing image.

## C3. Pipeline docs site (RMarkdown; `repo/current/docs/`)
- `docs/pipeline_technical.Rmd` — add the form's **technical write-up** as a Part-1/2/3 section, and
  add the tool to the software/dependencies table.
- `docs/pipeline_guide.Rmd` — add the toggle to the analysis-toggle table and any path params to the
  path-parameter table (from **config parameters**); add the output-target row and the per-GWAS
  output-dir to the output tree (from **outputs**).
- `pipeline/config.yaml` is already self-documenting via the contributor's comments.
- Re-knit and commit the HTML: `rmarkdown::render_site("docs")` (regenerates `docs/*.html` +
  `site_libs/`; GitHub Pages serves the committed HTML — there is no CI build).

## C4. Release
New submissions pin to a cut pipeline version; ship the method to prod by cutting a release
(`bin/cut_release.sh`) once the pipeline + viewer changes are released. The toggle then also appears
for that version (its `config.yaml` carries the key).
