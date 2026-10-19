# GenoDisc method contribution — GCTB SBayesRC

Completed from `pipeline/docs/method_contribution_template.md`. This is the contribution package
the maintainer integrates into the website submission panel and docs pages. The pipeline analysis
itself is already implemented and committed (pipeline `dev`, commit `0572bee`).

---

## 1. Method identity
- **Display name:** GCTB SBayesRC
- **One-line summary:** Estimates genome-wide SNP-based heritability and polygenicity from GWAS
  summary statistics using a Bayesian mixed model with functional annotations.
- **Config flag:** `sbayesrc`
- **Snakemake rule(s):** `sbayesrc`
- **Rule module:** `pipeline/rules/sbayesrc.smk` (`include:` added to `pipeline/Snakefile`)
- **External tool + version:** **gctb 2.5.5** — the official static Linux binary from the GCTB site
  (`https://cnsgenomics.com/software/gctb/download/gctb_2.5.5_Linux.zip`; host redirects to
  gctbhub.cloud.edu.au). Staged by `rule download_gctb` →
  `{resdir}/software/gctb/gctb_2.5.5_Linux/gctb`. (Deliberately the official binary, not GenoPred's
  gctb 2.5.2 and not the GitHub SBayesRC R package.)
- **References / citations:** Zheng, Z. et al. (2024) "Leveraging functional genomic annotations and
  genome coverage to improve polygenic prediction of complex traits within and between ancestries."
  *Nature Genetics* 56:767–777. doi:10.1038/s41588-024-01704-y. Tool: GCTB
  (https://cnsgenomics.com/software/gctb/), SBayesRC (https://github.com/zhilizheng/SBayesRC).

## 2. Config parameters
| Key | Type | Default | Description | Example |
|---|---|---|---|---|
| `sbayesrc` | T/F | `F` | Run GCTB SBayesRC to estimate SNP-h² and polygenicity (genome-wide). | `T` |
| `sbayesrc_ldref` | path | explicit (NA→GenoPred) | Directory of per-population eigen-decomposed LD references: `<dir>/<POP>/{ldm.info,snp.info,block*.eigen.bin}`. Shared with GenoPred. `NA` falls back to GenoPred's staged copy. | `/scratch/prj/neurohackpain/GenoPred/pipeline_resources/data/sbayesrc_ref` |
| `sbayesrc_annot` | path | explicit (NA→GenoPred) | GCTB functional annotation file (header `SNP Intercept ...`), i.e. baseline 2.2. `NA` falls back to GenoPred's staged copy. | `/scratch/.../sbayesrc_annot/annot_baseline2.2.txt` |

## 3. Submission-panel metadata
- **Form field name:** `sbayesrc`
- **Display label** (`ANALYSIS_LABELS`): `SBayesRC (SNP-h² & polygenicity)`
- **Card title:** SBayesRC
- **Card description:** Estimate SNP-based heritability and polygenicity (number of non-zero-effect
  variants) with GCTB SBayesRC. Genome-wide; EUR/EAS/AFR only.
- **Form section / placement:** **GWAS QC** section, directly after the *LDSC heritability* toggle
  (it is an alternative SNP-h² estimator; its output appears beside LDSC in the viewer).
- **Version-gating config key** (`ANALYSIS_CONFIG_KEYS`): `sbayesrc` → config key `sbayesrc`.
  (The toggle shows only for pipeline versions whose `config.yaml` contains `sbayesrc:` — i.e. `dev`
  now, and any release cut after this. Add `data-version-key="sbayesrc"` to the card so it is
  correctly hidden for older releases.)
- **Dependencies** (`ANALYSIS_DEPENDENCIES`): none (standalone, like `clump`/`finemap`).
- **Supported populations / constraints:** EUR, EAS, AFR only (published eigen-LD references).
  AMR/SAS GWAS are skipped (the report iterates a `gwas_list_df_sbrc` population subset).

## 4. Compute cost
- **Snakemake rule name (RULE_COST key):** `sbayesrc`
- **Typical wall-time (per GWAS):** ~1.5 h (whole-genome; chain length 3000)
- **Peak memory:** ~32 GB (≈32768 MB)
- **Per-chromosome or genome-wide:** genome-wide → `per_chr = False`
- **Threads:** 4
  (For `resource_estimates.py`: `RULE_COST['sbayesrc'] = {'low_s': 5400.0, 'high_s': 5400.0,
  'max_rss_mb': 32768, 'per_chr': False}`, and in `_rule_load`: `if flat.get('sbayesrc'): add('sbayesrc')`.)

## 5. Outputs
- **Final per-GWAS target file:** `{outdir}/results/{gwas}/sbayesrc/{gwas}.GW.sbayesrc.clean.csv`
  (one row: `gwas,h2_mode,h2_mean,h2_sd,h2_hpd_low,h2_hpd_high,n_post_samples,n_nonzero,n_snp_total,pi_est`).
  Also written alongside: `{gwas}.GW.sbayesrc.h2_samples.txt` (post-burn-in hsq draws for the density
  plot), plus `{gwas}.parRes` and `{gwas}.sbrc.stdout.log` for provenance.
- **Packaged bundle location:** `gwas_qc$sbayesrc_dat` — `$val` (the scalars above) + `$h2_samples`
  (the posterior draws). Nested in the existing `gwas_qc` block, so **no `reader.R` `.gd_block_ids`
  change** is needed.
- **Kept vs discarded:** keep `.parRes`, the stdout log, the clean.csv and h2_samples. Discard the
  large `.ma`/`.imputed.ma` and the per-SNP `.snpRes`.

## 6. Analyses-catalogue entry (docs_webapp.html table)
- **Category:** GWAS QC
- **"What it does" (one line):** GCTB SBayesRC — annotation-informed Bayesian estimate of SNP-based
  heritability and polygenicity (number of non-zero-effect variants). EUR/EAS/AFR.

## 7. Viewer presentation (docs_results.html)
- **Viewer tab:** SNP-h² & rG
- **Sub-tab name:** SBayesRC (h² & polygenicity)
- **What it shows:** For each GWAS: SNP-h² reported as the **posterior mode** with a **95% HPD**
  credible interval and the **posterior mean (SD)**; plus polygenicity as the **number of
  non-zero-effect SNPs** and the proportion **π**. A **posterior-density plot** of SNP-h² is shown
  with the mode (solid line) and 95% HPD bounds (dashed).
- **How to read it:** SNP-h² is bounded at 0 and its posterior is skewed when the true value is near
  zero, so the **mode** is the point estimate and the **HPD** is the uncertainty measure (prefer
  these to mean±SE near zero). The reported **SD is the standard error** — it is the standard
  deviation of the MCMC samples and is **not** divided by √N. Use the density plot to judge how
  symmetric the posterior is (symmetric → mean ≈ mode). SBayesRC's annotation-informed estimate
  complements the LDSC estimate on the adjacent sub-tab.
- **Caveats:** Genome-wide estimate on the HapMap3 SNP set; **EUR/EAS/AFR only**. This is a
  heritability/polygenicity summary, **not a fine-mapping output** (no per-variant credible sets).
- **Figure:** screenshot of the "SBayesRC (h² & polygenicity)" sub-tab for the ALS (or COAD01)
  example showing the table + posterior-density plot (light theme). Capture with
  `scripts/docs/capture_viewer_screenshots.py` once the viewer sub-tab is released to the public
  Shiny server; add a `screenshots.json` entry with `"tab": "SNP-h² & rG"`,
  `"sub_tab": "SBayesRC (h² & polygenicity)"`.

## 8. Technical method write-up (pipeline_technical.Rmd — Part 1, after LDSC)
GCTB **SBayesRC** (Zheng et al. 2024) estimates SNP-based heritability and polygenicity from GWAS
summary statistics using a Bayesian multi-component mixed model that incorporates functional genomic
annotations and a genome-wide eigen-decomposed LD reference. GenoDisc runs the official **gctb 2.5.5**
binary in two steps per GWAS: (1) `--impute-summary` imputes/QCs the summary statistics to the LD
reference SNP set, then (2) `--sbayes RC` fits the model with the baseline annotation. SNP-h² is
summarised from the post-burn-in MCMC trace of `hsq` as the posterior **mode** (robust for the
skewed near-zero posterior) with a **95% HPD** interval and posterior **mean/SD** (the SD is the
posterior SD, i.e. the standard error). Polygenicity is the posterior-mean number of
non-zero-effect variants (sum of the non-null mixture components) and its proportion π. The LD
reference (UKB HapMap3 eigen decomposition) is chosen per-GWAS from the population column
(EUR/EAS/AFR); other populations are skipped. The eigen-LD and annotation references are shared with
the GenoPred pipeline. Results are stored in `{outdir}/results/{gwas}/sbayesrc/`. Posterior-summary
presentation follows advice from J. Zeng (GCTB author).

## 9. Testing / verification
- **Test config + GWAS:** `pipeline/tests/config_sbayesrc_wg.yaml` (whole-genome COAD01, SBayesRC +
  LDSC); also `pipeline/tests/config_test_sbayesrc.yaml` (chr22 smoke test).
- **How it was run:** on the CREATE HPC via SSH + `sbatch` (the env can't run on the dev VM); target
  `snakemake tests/output_sbrc_wg/results/bundle.tar.gz --configfile tests/config_sbayesrc_wg.yaml
  --use-conda -c8`.
- **Observed result (COAD01, whole genome):** SNP-h² mode **0.035**, 95% HPD **[0.033, 0.039]**,
  posterior mean (SD) 0.0356 (0.0016), **~10.6k** non-zero SNPs of 1,154,502 (π ≈ 0.009); LDSC h²
  0.066 for comparison. End-to-end bundle packaged with `gwas_qc$sbayesrc_dat`. (The chr22 smoke test
  yields h²≈0 as expected — a mechanics check only.)
