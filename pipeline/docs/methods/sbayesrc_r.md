# GenoDisc method contribution — SBayesRC (R package)

A completed method-contribution form for the **SBayesRC R-package** method (`sbayesrc_r`). This is an
**alternative implementation** of SBayesRC to the gctb-binary method (`sbayesrc`, see
`methods/sbayesrc.md`): same model and same SNP-h²/polygenicity output, but run through the **SBayesRC
R package** (zhilizheng, v0.2.6) instead of the gctb 2.5.5 binary. The R package uses the same
eigen-decomposed LD reference, and is a faster engine (measured ~1.7× at matched threads, and it
multithreads the MCMC where the binary does not).

## 1. Method identity
- **Flag:** `sbayesrc_r`
- **Rule module:** `pipeline/rules/sbayesrc_r.smk`
- **Worker:** `pipeline/scripts/sbayesrc_r.R`
- **External tool + version:** the **SBayesRC R package** v0.2.6 (`github.com/zhilizheng/SBayesRC`),
  installed from the release tarball by `install_sbayesrc_rpkg` (`rules/dependencies.smk`). (A gctb
  binary is *not* used by this route.)
- **References / citations:** Zheng et al. 2024, *Nature Genetics* (`10.1038/s41588-024-01704-y`).
- **Tool licence:** MIT (SBayesRC R package).
- **Reference-data licence / redistribution:** eigen-LD reference is UK-Biobank-derived — download from
  source, do not re-host (same panel as the gctb method; on CREATE it reuses GenoPred's staged copy).
- **Citation users should give:** Zheng et al. 2024 (as above).
- **Maintainer / contact:** Oliver Pain (oliver.pain@kcl.ac.uk).

## 2. Config parameters
| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `sbayesrc_r` | T/F | `F` | Run the SBayesRC R-package method. |
| `sbayesrc_r_ldref` | path | `NA` | Per-population eigen-LD dir (`<dir>/<POP>/`). `NA` = the copy staged in `{resdir}/data/sbayesrc_r_ref` by `download_sbayesrc_r_ref`. On CREATE, set to GenoPred's staged `sbayesrc_ref`. |
| `sbayesrc_r_annot` | path | `NA` | Baseline2.2 annotation. `NA` = the shared `{resdir}/data/sbayesrc_annot/annot_baseline2.2.txt`. |
| `sbayesrc_r_threads` | int | `10` | OpenMP threads (`OMP_NUM_THREADS`); the R package multithreads the MCMC. |
| `sbayesrc_r_seed` | int | `1` | R-level seed (`set.seed`). |

## 3. Submission-panel metadata
- **Form field (`forms.py`):** `sbayesrc_r = forms.BooleanField(required=False)`.
- **`ANALYSIS_LABELS`:** `"sbayesrc_r": "SBayesRC (R package)"`.
- **Card (`_job_form.html`):** GWAS QC section, after the SBayesRC (GCTB) card; title "SBayesRC (R package)",
  desc "Faster R-package implementation of SBayesRC: genome-wide SNP-h² + polygenicity. EUR/EAS/AFR only."
  Version-gated: `data-version-key="sbayesrc_r"`.
- **`ANALYSIS_CONFIG_KEYS`:** `"sbayesrc_r": "sbayesrc_r"`.
- **`ANALYSIS_DEPENDENCIES`:** none (standalone).
- **Populations:** EUR/EAS/AFR (reuses `gwas_list_df_sbrc`).

## 4. Compute cost (`RULE_COST` + `_rule_load`)
- **Rule name:** `sbayesrc_r`. Genome-wide, one run per GWAS, `per_chr: False`.
- **Measured (ADHD, whole-genome, EUR):** ~27 min at 10 threads (default), ~50 min at 4 threads; peak
  RSS ~8 GB. (vs the gctb method's 1h26m at 4 threads on the same data + identical LD reference.)
- **Suggested `RULE_COST['sbayesrc_r']`:** `{'low_s': 1650.0, 'high_s': 1650.0, 'max_rss_mb': 8192, 'per_chr': False}`
  (10-thread default). `_rule_load`: `if flat.get('sbayesrc_r'): add('sbayesrc_r')`.

## 5. Outputs
- **Final target:** `{outdir}/results/{gwas}/sbayesrc_r/{gwas}.GW.sbayesrc_r.clean.csv`
  (columns: `gwas, h2_mean, h2_sd, n_nonzero, n_snp_total, pi_est`). The SBayesRC `.par` + annotation
  enrichment + tidy/impute logs are retained in that dir for provenance; per-SNP effects (PGS weights)
  are NOT kept.
- **Packaged bundle:** `gwas_qc$sbayesrc_r_dat` (`process_sbayesrc_r`).

## 6. Analyses-catalogue entry (docs_webapp.html)
GWAS QC · "SBayesRC (R package)" · "Faster R-package implementation of SBayesRC — genome-wide SNP-based
heritability and polygenicity. EUR/EAS/AFR only."

## 7. Viewer presentation (docs_results.html)
An `<h3 id="sbayesrc-r">` under the "SNP-h² & rG" tab: a table of SNP-h² (posterior **mean ± SD**) and
polygenicity (non-zero SNPs, π). Unlike the GCTB SBayesRC tab there is **no posterior mode / HPD /
density** — the R package reports `.par` summary statistics, not a retained MCMC trace.

## 8. Technical method write-up (pipeline_technical.Rmd — after the SBayesRC section)
The SBayesRC R package (Zheng et al. 2024) runs the same annotation-informed Bayesian mixed model as the
gctb-binary method, via `SBayesRC::tidy()` → `impute()` → `sbayesrc()` over the same UK-Biobank HapMap3
eigen-LD reference and baseline2.2 annotation. It is multithreaded (OpenMP) and faster than the binary.
GenoDisc retains the `.par` parameter file (which the package writes but GenoPred discards) to report
SNP-h² (posterior mean/SD) and polygenicity. Results in `{outdir}/results/{gwas}/sbayesrc_r/`.

## 9. Testing / verification
- **Test config:** `pipeline/tests/config_sbayesrc_r.yaml` (chr22 smoke; reuses GenoPred's refs).
- **Observed (ADHD, whole genome, 10 threads):** SNP-h² mean 0.270 (SD 0.022), 14,284 non-zero SNPs
  (π ≈ 0.012); ~27 min. Consistent with the gctb method (0.228) to within the expected engine/tuning
  difference (both plausible for ADHD; same LD reference).

## 10. Known failure modes
- **Unsupported population:** EUR/EAS/AFR only (AMR/SAS skipped).
- **SBayesRC tuning failures:** the worker wraps `sbayesrc()` in GenoPred's retry logic (`bTune = FALSE`,
  `tuneStep = c(0.995, 0.9, 0.8, 0.7, 0.6)`) for the "all correlations negative" / "best parameter at
  minimum threshold" / "invalid tune outputs" errors.
- **First-run env build:** the dedicated `sbayesrc_r` conda env + SBayesRC compile is a one-off cost;
  reused thereafter from the shared `--conda-prefix`.
- **SBayesRC vs LDSC / vs the gctb method:** different modelling assumptions/engines can give
  substantially different SNP-h² estimates; an LDSC intercept far from 1 warrants caution interpreting any.
