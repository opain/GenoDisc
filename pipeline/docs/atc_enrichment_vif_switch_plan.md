# Plan: VIF-OLS ATC drug-class enrichment — interim "all three" + production switch

## Context

The validation study (`docs/atc_enrichment_estimator_validation.Rmd`) concluded that the ATC
drug-class enrichment should move from the floor-0.1 GLS to **CAMERA-style variance-inflated OLS
(VIF-OLS)**, for both engines (MAGMA, TWAS-GSEA). Justification: the estimand (in-class vs rest,
correctly variance-inflated) and a demonstrated GLS pathology (out-of-class leakage inflates small
mixed classes, dilutes large coherent ones). VIF-OLS recovered the CAD lipid control on MAGMA
(rank 2, FDR 0.008) where GLS diluted it; GLS found 0 FDR-significant classes across 18 runs.

**Interim decision (this change): keep all three methods** (VIF-OLS, GLS, Wilcoxon) selectable in the
output, so users — and we — can compare them directly before committing to a hard switch.

## Estimator (VIF-OLS)

Fit `y ~ [1, class_indicator(s), size, log(size)]` by OLS; inflate the class-coefficient SE for
within-class correlation:
`SE = SE_OLS * sqrt(1 + (n1 - 1) * rho_bar)`, `rho_bar = max(0, mean of in-class off-diagonal Sigma)`.
Tails match production: **MAGMA one-sided** (enrichment), **TWAS two-sided**. Direction labels reuse
the engine conventions (MAGMA Enriched/Depleted; TWAS Matches/Opposes, `Reversal_Z = -T`).

**Alignment (the validated, correct path — differs from the current GLS):** align the per-drug
statistics to `drugcorr.rds` by **CID/full drug string**, and explode multi-ATC drugs into **all**
their L2/L3/L4 classes. (The current TWAS GLS keys by a single 7-char `res$ATC` and drops multi-ATC
drugs — bug #2 below. VIF uses the correct join from the outset.)

## Part A — interim implementation (all three methods)

### Pipeline (`repo/current/pipeline`, branch `dev`)

New format scripts (mirror the GLS scripts; emit L2/L3/L4):
- `scripts/format_twas_gsea_drugtargetor_vif_results.R` → `.../twas/drugtargetor/twas_gsea_drugtargetor_vif_{l2,l3,l4}_{panel}_res.csv`
  columns `Code, Name, N_Drugs, Estimate, SE, T, P, FDR, Direction, Reversal_Z`.
- `scripts/format_magma_drugtargetor_vif_results.R` → `.../magma/magma_drug_targetor_vif_{l2,l3,l4}_res.csv`
  columns `Code, Name, N_Drugs, Estimate, SE, T, P, FDR, Direction`.

New rules (mirror the GLS rules, same inputs — drugcorr.rds + per-drug stats + ATC labels):
- `rules/twas_pwas.smk`: `format_twas_gsea_drug_targetor_vif_results` + `..._all_panel`.
- `rules/magma.smk`: `format_magma_drug_targetor_vif_results` (prereq `compute_magma_drugcorr`).

Gating: **reuse the existing GLS flags** (`drug_targetor_atc_gls` for TWAS, `magma_drugtargetor_gls`
for MAGMA) in `rules/report.smk` — no new config flags, so VIF appears wherever GLS already runs and
**no Django change is required**.

### Packaging (`scripts/functions/package_results_functions.R`, `package_results.R`, `reader.R`)

Mirror the GLS readers:
- `read_twas_gsea_atc_vif` → block `tx/atc_vif` (gate `drug_targetor_atc_gls`).
- `read_magma_atc_vif` → block `tx/atc_vif_magma` (gate `magma_drugtargetor_gls`; `Panel = "MAGMA"`).
- Register `tx/atc_vif`, `tx/atc_vif_magma` in `.gd_block_ids` (`reader.R`); assign in `package_results.R`.

### Shiny (`shiny_apps/output_explorer_embed`)

- Add source `"vif"` to `atc_src_choices` (`mod_enrichment.R:734`), labelled **"VIF-OLS (recommended)"**,
  default-selected; keep `"gls"` and `"legacy"`.
- Add a `"vif"` branch in `build_atc_summary_data` (`utils_data.R`) with labels `TWAS-GSEA (VIF-OLS)` /
  `MAGMA (VIF-OLS)` (these match the existing facet-rank + directionality regexes).
- Add `"vif"` branches in `atc_magma_tbl` and `atc_twas_tbl` (`mod_enrichment.R`).
- Add `has_atc_vif` / `has_atc_vif_magma` helpers + level availability; show the level selector for
  `vif` as well as `gls` (`conditionalPanel` at `mod_enrichment.R:1158`).
- Add the two block ids to the Shiny copy of `reader.R`.
- Evidence drill-down is source-agnostic (routes by engine via `tx/drug`) — no change.

### Run + view
Re-run SCZ (SCHI04) end-to-end on the new code (`dev`), which builds the VIF CSVs (reusing cached
drugcorr/stats) and re-packages the bundle with the new blocks; load the bundle in the dev Shiny
Server (`:3838`) for review.

## Part B — production switch (later, separate change)

1. Promote VIF-OLS to the headline method (default + "recommended"); keep GLS + Wilcoxon selectable.
2. **Fix the two alignment bugs** in the shared correlation path (affects GLS too):
   (i) `drugcorr.rds` vs `competitive.clean.csv` drug-set mismatch; (ii) primary-ATC-only keying in
   the TWAS GLS (`res$ATC`, drops multi-ATC drugs). VIF already uses the correct join; back-port to GLS.
3. Docs (pipeline guide/README): describe VIF-OLS as recommended, cite CAMERA (Wu & Smyth 2012);
   reference the validation report.
4. Version release + prod / shinyapps deploy (coordinated). No Django toggle change (VIF rides the
   GLS flags).

## Verification
- `snakemake -n` DAG dry-run builds the new VIF targets for SCHI04.
- Each new R script parses/runs standalone via system `/usr/bin/Rscript` on SCHI04 outputs.
- `package_results` produces `tx/atc_vif` + `tx/atc_vif_magma`; readers return expected columns.
- Shiny: bundle loads; the 3-way toggle switches VIF/GLS/Wilcoxon; summary facets MAGMA→TWAS; tables +
  evidence render for VIF rows.
- Spot-check SCZ: VIF `N05A` ranks well above the GLS ranking; MAGMA VIF `C08C`/`N03A` FDR-significant.
