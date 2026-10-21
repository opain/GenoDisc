# GenoDisc method contribution — LDAK SumHer

A completed method-contribution form for **LDAK SumHer** SNP-heritability (`sumher`). This is a **third,
independent** SNP-h² estimator alongside LDSC (`ldsc`) and the two SBayesRC routes (`sbayesrc`,
`sbayesrc_r`). It runs LDAK's `--sum-hers` with the **BLD-LDAK** heritability model over a **pre-computed
HapMap3 tagging file**, on the observed scale. SumHer assumes a *different* heritability model from LDSC
and SBayesRC, so its estimate is **not expected to agree** with theirs — the three are complementary, not
cross-checks of one another.

## 1. Method identity
- **Flag:** `sumher`
- **Rule module:** `pipeline/rules/sumher.smk`
- **Worker:** `pipeline/scripts/sumher.R`
- **External tool + version:** the **LDAK** binary `ldak7.01.linux`, pinned to commit
  `8934ac25b2264e926219da52c05beddd29a6546a` of `github.com/dougspeed/LDAK` and staged by `download_ldak`
  (`rules/dependencies.smk`). dougspeed re-versions the Linux executable on every change, so the pin is a
  commit SHA (not `main`). Verified: the binary prints **"Version 6"**, sha256
  `c40268ad5e24d4de6ba64876e22411181f894e45637aad2f7ed89babcef7cbec` (pinned + checked in the rule), and
  accepts `--sum-hers --tagfile --summary --check-sums --cutoff --max-threads`.
- **References / citations:** SumHer — Speed & Balding 2019, *Nature Genetics* 51:277–284
  (`10.1038/s41588-018-0279-5`). BLD-LDAK model — Speed, Holmes & Balding 2020, *Nature Genetics*
  52:458–462 (`10.1038/s41588-020-0600-y`, "Evaluating and improving heritability models using summary
  statistics").
- **Tool licence:** GPLv3 (LDAK; `LICENCE.txt` in the repo). Redistribution of the binary must preserve
  GPLv3 terms and offer corresponding source (the repo's `Source_Code/` satisfies this).
- **Reference-data licence / redistribution:** pre-computed BLD-LDAK HapMap3 tagging files are published by
  the LDAK author at `genetics.ghpc.au.dk/doug/TaggingFiles/` — downloaded from source, not re-hosted.
- **Citation users should give:** Speed & Balding 2019 + Speed, Holmes & Balding 2020 (as above).
- **Maintainer / contact:** Oliver Pain (oliver.pain@kcl.ac.uk).

## 2. Config parameters
| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `sumher` | T/F | `F` | Run LDAK SumHer SNP-heritability. |
| `sumher_ldak` | path | `NA` | LDAK executable. `NA` = the copy staged at `{resdir}/software/ldak/ldak7.01.linux` by `download_ldak`. |
| `sumher_ref` | path | `NA` | Base dir of tagging files (`<base>/<ldakpop>/bld.ldak.hapmap.<ldakpop>.tagging`). `NA` = `{resdir}/data/sumher_tagging/bld_ldak/hapmap` staged by `download_sumher_tagging`. |
| `sumher_cutoff` | float | `0.01` | LDAK `--cutoff`: excludes large-effect predictors (they bias heritability models). LDAK's default. |
| `sumher_max_missing` | float | `0.2` | Warn threshold on the fraction of tagging predictors absent after GenoDisc QC (LDAK's author advises <0.2 for pre-computed taggings). |
| `sumher_threads` | int | `1` | LDAK `--max-threads`. SumHer from a pre-computed tagfile is light, so 1 suffices. |

## 3. Submission-panel metadata (Django front-end — follow-on, NOT wired in this task)
- **Form field (`forms.py`):** `sumher = forms.BooleanField(required=False)`.
- **`ANALYSIS_LABELS`:** `"sumher": "LDAK SumHer (SNP-h²)"`.
- **Card (`_job_form.html`):** GWAS QC section, after the SBayesRC cards; title "LDAK SumHer", desc
  "Independent SNP-based heritability via LDAK's SumHer (BLD-LDAK model). EUR/EAS/SAS/AFR."
  Version-gated: `data-version-key="sumher"`.
- **`ANALYSIS_CONFIG_KEYS`:** `"sumher": "sumher"`.
- **`ANALYSIS_DEPENDENCIES`:** none (standalone).
- **Populations:** EUR/EAS/SAS/AFR (uses `gwas_list_df_sumher`; AMR skipped — no tagging file).

## 4. Compute cost (`RULE_COST` + `_rule_load`)
- **Rule name:** `sumher`. Genome-wide, one run per GWAS, `per_chr: False`, single-threaded by default.
- **Measured (COAD01 whole-genome, EUR, 1 thread):** 65 s wall / 46 s CPU; peak RSS ~2.1 GB (Snakemake
  benchmark). Reading the ~1.2M-row tagging file for the missingness count + the LDAK regression dominate.
- **Suggested `RULE_COST['sumher']`:** `{'low_s': 65.0, 'high_s': 65.0, 'max_rss_mb': 2560, 'per_chr': False}`.
  `_rule_load`: `if flat.get('sumher'): add('sumher')`.

## 5. Outputs
- **Final target:** `{outdir}/results/{gwas}/sumher/{gwas}.sumher.clean.csv` (columns: `gwas, h2, h2_sd,
  model, n_tagging, n_matched, prop_missing, missing_warning, n_excluded_cutoff, logl_null, logl_alt,
  cutoff, ldak_banner`). The LDAK `.hers`/`.extra`/`.enrich`/`.cats`/`.share` + a `.ldak.log` are retained
  in that dir for provenance.
- **Packaged bundle:** `gwas_qc$sumher_dat` (`process_sumher`) = `list(val, enrich, meta)`.

## 6. Analyses-catalogue entry (docs_webapp.html)
GWAS QC · "LDAK SumHer" · "Independent SNP-based heritability via LDAK's SumHer (BLD-LDAK model) over
pre-computed HapMap3 tagging. A different heritability model from LDSC/SBayesRC. EUR/EAS/SAS/AFR."

## 7. Viewer presentation (docs_results.html)
An `<h3 id="sumher">` under the "SNP-h² & rG" tab: a table of SNP-h² (point estimate ± SD), the model
(BLD-LDAK), and reference-SNP coverage (matched / total). A neutral line notes the three heritability
methods are not expected to agree; a visible warning appears when reference-SNP missingness exceeds the
configured threshold. Enrichments are retained but **not surfaced** yet.

## 8. Technical method write-up (pipeline_technical.Rmd — after the SBayesRC section)
LDAK SumHer (Speed & Balding 2019) estimates SNP-heritability from GWAS summary statistics under the
BLD-LDAK heritability model (Speed, Holmes & Balding 2020), which combines Baseline-LD annotations with
LDAK's SNP-weighting. The worker builds an LDAK summary file (`Predictor A1 A2 n Z A1Freq`, rsID-keyed,
`Z = BETA/SE`) from `{gwas}.cleaned.gz` and runs `ldak --sum-hers … --tagfile <BLD-LDAK HapMap3 tagging>
--check-sums NO --cutoff 0.01` on the observed scale (no intercept/genomic-control; no
prevalence/ascertainment). SNP-h² and its SD are read from the final (all-SNPs) row of `.hers`. Results in
`{outdir}/results/{gwas}/sumher/`.

## 9. Testing / verification
- **Test config:** `pipeline/tests/config_test_sumher.yaml` (COAD01chr22, EUR→gbr tagging). A chr22-only
  summary against a genome-wide tagging file intentionally trips both the missingness-warning
  (prop_missing 0.985) and non-convergence (`converged=FALSE`) paths; `h2_sd` correctly NA (−9999 sentinel).
- **Whole-genome (COAD01, EUR):** SumHer h² **−0.01467 (SD 0.00484)**, `converged=TRUE`, coverage
  1,163,623/1,168,975 (prop_missing 0.46%, no warning), 0 excluded by cutoff. For comparison, in the same
  run LDSC h² = 0.0664 (SE 0.0044), intercept 0.8868. (Numbers only — a cross-method difference is expected
  under different heritability models, not characterised here as explained or resolved.)

## 10. Known failure modes
- **High reference-SNP missingness:** if a large fraction of the tagging file's predictors are absent after
  GenoDisc QC (e.g. a sparse or single-chromosome GWAS), the estimate is unreliable — surfaced via
  `prop_missing` + `missing_warning` (default threshold 0.2, per LDAK's author).
- **Large-effect loci:** heritability models are biased by very-large-effect predictors; `--cutoff` (0.01)
  down-weights/excludes them. The excluded count is recorded (`n_excluded_cutoff`).
- **Unsupported population:** EUR/EAS/SAS/AFR only (AMR skipped — no BLD-LDAK HapMap3 tagging).
- **Ambiguous sample size:** observed-scale only; GenoDisc uploads may carry total N or Neff, so
  liability-scale conversion (prevalence/ascertainment) is deliberately not attempted.
- **SumHer vs LDSC vs SBayesRC:** different heritability models — divergent estimates are expected and do
  not by themselves indicate an error.
