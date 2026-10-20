# AGENTS.md

Guidance for AI coding tools working in the GenoDisc pipeline. The full recipe is in
[`CONTRIBUTING.md`](CONTRIBUTING.md) — read it before adding or changing a pipeline analysis or
dataset. (`CLAUDE.md` holds the generic coding guidelines and points here.)

## Adding an optional analysis (checklist — see CONTRIBUTING.md for detail)
1. `config.yaml`: add `<name>: F` (+ any ref-path keys).
2. `rules/dependencies.smk`: a `download_*`/`install_*` rule staging into `{resdir}/...` + a
   `resource_inputs` entry gated on `config["<name>"]=="T"`. Prefer official download URLs.
3. `rules/<name>.smk` (+ `include` in `Snakefile`) and `scripts/<name>.R` — mirror clump/cojo/finemap;
   pick the reference from the `population` column via a `params` lambda; declare `threads` from a
   config key for multi-threaded tools.
4. `rules/report.smk`: append the final output under `if config["<name>"]=="T":`.
5. `scripts/functions/package_results_functions.R` + `scripts/package_results.R`: `process_<name>()`
   into the bundle (gated; returns `NULL` when absent).
6. Viewer (if surfaced): flag in `R/utils_config.R`, gate in `app.R`, render in a `R/mod_*.R`.
7. Tests: dedicated `tests/config_<name>.yaml`; enable in `prepare_resources.yaml`.
8. Fill the contribution form: copy `pipeline/docs/method_contribution_template.md` (or
   `dataset_contribution_template.md`) → `pipeline/docs/methods/<name>.md` and complete it
   (submission-panel metadata, licence/citation/contact/failure modes, and user-facing docs prose).
   The maintainer uses it to wire the website + docs. Worked example: `pipeline/docs/methods/sbayesrc.md`.

## Invariants
- Flags are the strings `"T"`/`"F"` (not booleans): `config["x"]=="T"`,
  `read_param(cfg,'x',return_obj=F)=="T"`.
- Nesting a packaged result under an existing block id (e.g. `gwas_qc$<name>_dat`) needs no
  `.gd_block_ids` change; a new top-level block does (edit both copies of `reader.R`).
- Don't rely on training data for external-tool versions/flags — check the official docs and confirm
  against the real binary (e.g. GCTB's `.parRes`/MCMC output differs from older docs).
- Stage reference data into `resdir` from a stable official URL; encode version/build/ancestry in the
  path. Reference-data paths/credentials and site-specific run details do not belong in this repo.
