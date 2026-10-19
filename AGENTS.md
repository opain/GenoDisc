# AGENTS.md

Guidance for AI coding tools working in the GenoDisc pipeline. The full recipe is in
[`CONTRIBUTING.md`](CONTRIBUTING.md) — read it before adding or changing a pipeline analysis.
(`CLAUDE.md` holds the generic coding guidelines and points here.)

## Key facts

- **The pipeline only runs on the CREATE HPC**, inside the `genodisc` conda env. The dev VM cannot
  execute it (env binaries resolve `/cephfs`; the Snakefile aborts unless the env is active). Run
  tests by SSHing to `hpc.create.kcl.ac.uk` (the mechanism in `modules/hpc_interface.py`) and
  `sbatch`-ing a wrapper that `source bin/activate_genodisc.sh`.
- `pipeline/repo` note: `…/repo/current` is a **symlink to the dev worktree** (`GenoDisc_v1`);
  edits under `current/` land in the dev checkout. Released versions are frozen worktrees under
  `…/repo/versions/<vX.Y.Z>`; `cut_release.sh` cuts them. Do not edit a versions/ worktree.

## Adding an optional analysis (checklist — see CONTRIBUTING.md for detail)

1. `config.yaml`: add `<name>: F` (+ ref-path keys if any).
2. `rules/dependencies.smk`: a `download_*`/`install_*` rule + a `resource_inputs` entry gated on
   `config["<name>"]=="T"`. Share big references with GenoPred where possible.
3. `rules/<name>.smk` (+ `include` in `Snakefile`) and `scripts/<name>.R` — mirror clump/cojo/finemap;
   pick the LD reference from the `population` column via a `params` lambda.
4. `rules/report.smk`: append the final output under `if config["<name>"]=="T":`.
5. `scripts/functions/package_results_functions.R` + `scripts/package_results.R`: `process_<name>()`
   into the bundle (gated; returns `NULL` when absent).
6. Viewer (if surfaced): flag in `R/utils_config.R`, gate in `app.R`, render in a `R/mod_*.R`.
7. `tests/config_test.yaml` + `prepare_resources.yaml`: `<name>: T`.
8. Fill the contribution form: copy `pipeline/docs/method_contribution_template.md` →
   `pipeline/docs/methods/<name>.md` and complete it (submission-panel metadata + user-facing docs
   prose). The maintainer uses it to wire the website form + docs (CONTRIBUTING.md Part C). Worked
   example: `pipeline/docs/methods/sbayesrc.md`.

## Invariants

- Flags are the strings `"T"`/`"F"` (not booleans): `config["x"]=="T"`,
  `read_param(cfg,'x',return_obj=F)=="T"`.
- Nesting a packaged result under an existing block id (e.g. `gwas_qc$<name>_dat`) needs no
  `.gd_block_ids` change; a new top-level block does (edit both copies of `reader.R`).
- Don't rely on training data for external-tool versions/flags — check the official docs and
  confirm against the real binary (e.g. GCTB's `.parRes` format differs from older docs).
- Keep a running log when doing a large integration (see `pipeline/docs/sbayesrc_integration_log.md`
  for the SBayesRC example) and discuss pipeline changes with the maintainer before committing.
