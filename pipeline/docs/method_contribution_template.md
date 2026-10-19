# GenoDisc method contribution — template

Fill this in when contributing a new **optional analysis ("method")** to the GenoDisc pipeline.
Copy this file to `pipeline/docs/methods/<flag>.md` (where `<flag>` is your config flag, e.g.
`sbayesrc.md`) and complete every field. Leave no `TODO`s.

This form is the single source the **maintainer** uses to (a) add your method to the website's
job-submission *analysis panel* and (b) write the website + pipeline-repo documentation. You do
**not** edit the Django web app or the docs site yourself — the maintainer integrates this form.
First implement the analysis in the pipeline (see `CONTRIBUTING.md` → "For contributors"), then
complete the sections below. See `pipeline/docs/methods/sbayesrc.md` for a worked example.

---

## 1. Method identity
- **Display name:** <!-- short human name, e.g. "GCTB SBayesRC" -->
- **One-line summary:** <!-- one sentence: what it computes / why a user wants it -->
- **Config flag:** `<flag>` <!-- the config.yaml T/F key, e.g. sbayesrc -->
- **Snakemake rule(s):** <!-- rule name(s), e.g. sbayesrc, process_sbayesrc -->
- **Rule module:** `pipeline/rules/<flag>.smk` <!-- new include added to Snakefile -->
- **External tool + version:** <!-- e.g. gctb 2.5.5 (official binary); how it is obtained -->
- **References / citations:** <!-- papers + URLs -->

## 2. Config parameters
List every NEW `config.yaml` key the method adds (the T/F toggle AND any path/option params).
Feeds the `config.yaml` comments and the `pipeline_guide.Rmd` parameter tables.

| Key | Type | Default | Description | Example |
|---|---|---|---|---|
| `<flag>` | T/F | `F` | <!-- what enabling it does --> | `T` |
| <!-- e.g. <flag>_ref --> | path | `NA` | <!-- what it points at; NA behaviour --> | `/path/...` |

## 3. Submission-panel metadata (job-submission analysis toggle)
Everything the maintainer needs to add the toggle to the website form.
- **Form field name:** `<flag>` <!-- usually == config flag -->
- **Display label** (`ANALYSIS_LABELS`): <!-- e.g. "SBayesRC (SNP-h² & polygenicity)" -->
- **Card title:** <!-- the bold name on the toggle card -->
- **Card description** (one line, shown under the title): <!-- accurate, user-facing -->
- **Form section / placement:** <!-- which section of the panel, e.g. "GWAS QC" next to LDSC -->
- **Version-gating config key** (`ANALYSIS_CONFIG_KEYS`): `<flag>` <!-- the config.yaml key that
  proves a pipeline version supports this; the toggle shows only for versions whose config.yaml has it -->
- **Dependencies** (`ANALYSIS_DEPENDENCIES`): <!-- parent toggle(s) that must be on, or "none" -->
- **Supported populations / constraints:** <!-- e.g. EUR/EAS/AFR only; others skipped -->

## 4. Compute cost (for the submission-form resource estimate)
Feeds `resource_estimates.py` (`RULE_COST` + `_rule_load`).
- **Snakemake rule name (for RULE_COST key):** <!-- e.g. sbayesrc -->
- **Typical wall-time (per GWAS):** <!-- e.g. ~1.5 h -->
- **Peak memory:** <!-- e.g. ~32 GB -->
- **Per-chromosome or genome-wide:** <!-- per_chr True/False -->
- **Threads:** <!-- e.g. 4 -->

## 5. Outputs
- **Final per-GWAS target file(s):** <!-- e.g. {outdir}/results/{gwas}/<flag>/{gwas}.GW.<flag>.clean.csv -->
- **Packaged bundle location:** <!-- e.g. gwas_qc$<flag>_dat (nested in an existing block = no reader.R
  change) or a new top-level block id -->
- **Kept vs discarded:** <!-- which raw outputs are retained for provenance; which are dropped -->

## 6. Analyses-catalogue entry (website "choosing analyses" table)
One row for `docs_webapp.html`.
- **Category:** <!-- e.g. GWAS QC / Genetic Correlation / SNP Associations / ... -->
- **"What it does" (one line):** <!-- user-facing description for the catalogue table -->

## 7. Viewer presentation (website "Interpreting your results" page)
Drop-in prose for `docs_results.html`. Omit if the method produces no viewer output.
- **Viewer tab:** <!-- e.g. "SNP-h² & rG" -->
- **Sub-tab name:** <!-- e.g. "SBayesRC (h² & polygenicity)" -->
- **What it shows:** <!-- the columns/plots and their meaning -->
- **How to read it:** <!-- interpretation guidance -->
- **Caveats:** <!-- limitations the user must know -->
- **Figure:** <!-- describe the screenshot the maintainer should capture (which tab/state) -->

## 8. Technical method write-up (pipeline docs site)
A short methods paragraph for `docs/pipeline_technical.Rmd` (what it is, how it works, inputs,
outputs, references, results location). 3–8 sentences; user/method-facing, not implementation detail.

<!-- write the paragraph here -->

## 9. Testing / verification
- **Test config + GWAS used:** <!-- e.g. tests/config_<flag>.yaml on COAD01 (whole genome) -->
- **How it was run end to end:** <!-- command / HPC submission -->
- **Expected / observed result:** <!-- concrete numbers proving it works -->
