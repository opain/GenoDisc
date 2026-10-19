# SBayesRC integration — running log

Append-only log of every file touched and every decision/finding while adding GCTB SBayesRC
(SNP-h² + polygenicity) to the GenoDisc pipeline. Source material for CONTRIBUTING.md / AGENTS.md
/ CLAUDE.md once it works end to end on COAD01.

## Locked decisions (from user, 2026-10-18)
- Run standard genome-wide SBayesRC only for **SNP-h² (hsq) + polygenicity** (Pi / N non-zero).
  NOT fine-mapping; do not keep per-SNP PIP/effects; do not use GWFM mode.
- Present as a **new sub-tab in the viewer "SNP-h² & rG" tab**, beside LDSC.
- Tool: **official gctb binary v2.5.5** (host gctbhub.cloud.edu.au, `gctb_2.5.5_Linux.zip`).
  NOT GenoPred's gctb 2.5.2, NOT the GitHub SBayesRC R package.
- Populations EUR/EAS/AFR via the per-GWAS `population` column; skip AMR/SAS with a logged note.
- Reuse GenoPred's staged refs (no re-download):
  LD `/scratch/prj/neurohackpain/GenoPred/pipeline_resources/data/sbayesrc_ref/{EUR,EAS,AFR}/`,
  annot `/scratch/.../GenoPred/pipeline_resources/data/sbayesrc_annot/annot_baseline2.2.txt`.
- Do NOT commit until the user reviews the pipeline changes.

## Reference facts
- Official recipe (binary): `gctb --ldm-eigen <LDdir> --gwas-summary X.ma --impute-summary --out X --thread N`
  then `gctb --ldm-eigen <LDdir> --gwas-summary X.imputed.ma --sbayes RC --annot <annot> --out X --thread N`.
- `.ma` columns: `SNP A1 A2 freq b se p N`.
- Parse `<out>.parRes`: `hsq` row (Mean,SD) = SNP-h²; `Pi` / `NnzSnp` rows = polygenicity.
- GenoPred reference files: pgs_methods.smk:686-732, Scripts/pgs_methods/sbayesrc.R,
  dependencies.smk:1090-1198, envs/sbayesrc.yaml.

## Timeline
- 2026-10-18: Plan approved. Created this log. Beginning with read-only verification of GenoPred
  refs + template files before writing code.
- 2026-10-18: Verified GenoPred refs readable (ollie/smbmountaccess). Each of EUR/EAS/AFR has 591
  `blockN.eigen.bin` + `ldm.info` + `snp.info`. `ldm.info` header:
  `Block Chrom StartSnpIdx StartSnpID EndSnpIdx EndSnpID NumSnps` (tab-sep).
  Annotation `annot_baseline2.2.txt` (~2.0 GB) header is `SNP Intercept Coding_UCSC ...` — i.e.
  the gctb-binary `--annot` layout (SNP, Intercept, then annotations). **Risk #2 (annot format)
  largely resolved: GenoPred's annot is already binary-compatible.** Risk #1 (binary reads the
  eigen blocks) still to confirm at first run.
- 2026-10-18: Read templates. Cleaned sumstats (`{gwas}.cleaned.gz`) columns:
  `CHR A1 A2 FREQ P BETA SE N SNP BP REF.FREQ` (per-SNP N present) → `.ma` = SNP A1 A2 FREQ BETA SE P N
  renamed to `SNP A1 A2 freq b se p N`. `read_param` only auto-defaults NA for resdir/refdir, so
  SBayesRC ldref/annot NA→GenoPred defaults handled in rules/sbayesrc.smk + sbayesrc.R. report.smk
  and resource_inputs gate each analysis with `if config[x]=="T"`; report iterates gwas_list_df_eur
  — SBayesRC will iterate a new EUR/EAS/AFR subset `gwas_list_df_sbrc`. download_gcta (deps.smk:362)
  is the binary-download template. finemap.R is the worker-script template.
- 2026-10-18: Downloaded official `gctb_2.5.5_Linux.zip` (host redirects cnsgenomics→gctbhub);
  unzips to `gctb_2.5.5_Linux/gctb` (static ELF, "GCTB 2.5.5, last updated 12 Dec 2025", MIT).
  **PROBE CONFIRMED on the real binary**: `--ldm-eigen <dir> --gwas-summary X.ma --impute-summary
  --out X --thread N` accepted and **it reads GenoPred's EUR eigen dir directly** ("591 LD Blocks …
  1154522 SNPs"), matched 16696 of 90446 chr22 SNPs. **Risk #1 resolved.** Standard-SBayesRC flags
  (`--sbayes RC`, `--impute-summary`) exist in 2.5.5 (README's headline example is the newer `--gwfm
  RC` mode, which we are NOT using). Dev VM has ~13 GB RAM; full SBayesRC (2 GB annot + MCMC over
  1.15M SNPs; GenoPred budgets 20 GB) needs a SLURM compute node → end-to-end test via snakemake on
  a node, rule declares mem_mb≈32000.
- 2026-10-18/19: Ran the full 2-step recipe manually on the COAD01chr22 `.ma` against GenoPred's
  EUR eigen ref + annot. Impute → `<out>.imputed.ma` (~1 min). `--sbayes RC` with the 2 GB annot
  completed in ~2 min and fit in the VM's ~13 GB (so the chr22 smoke test can even run locally;
  a real genome-wide GWAS still wants a node — rule keeps mem_mb=32000).
  **Real `.parRes` format (gctb 2.5.5) differs from the docs**: header `Parameter Mean SD`; rows
  `NumSnp1..5, Vg1..5, SigmaSq, hsq, ResVar` — **NO `Pi`/`NnzSnp` row**. So:
    - SNP-h2 = `hsq` row (Mean, SD).
    - polygenicity = non-zero-effect variants = sum(NumSnp2..NumSnp5) (NumSnp1 = null component);
      pi = that / sum(NumSnp1..5).
  Rewrote the sbayesrc.R parser to this; validated on the real file → h2_est=3.4e-5, h2_se=1.45e-4,
  n_nonzero≈33.4, n_snp_total=16774, pi≈0.002 (h2≈0 expected: chr22-only data). clean.csv columns:
  gwas,h2_est,h2_se,n_nonzero,n_snp_total,pi_est. `process_sbayesrc()` reads these.

## Files touched (pipeline wiring)
- `config.yaml` — `sbayesrc: F`, `sbayesrc_ldref: NA`, `sbayesrc_annot: NA` (+ comments).
- `Snakefile` — `include: "rules/sbayesrc.smk"`.
- `rules/dependencies.smk` — new `rule download_gctb` (gctb 2.5.5, GCTA-pattern) + resource_inputs
  entry under `if config.get("sbayesrc")=="T"`.
- `rules/sumstat_qc.smk` — `gwas_list_df_sbrc` = population ∈ {EUR,EAS,AFR}.
- `rules/report.smk` — append sbayesrc clean.csv target (iterates gwas_list_df_sbrc).
- `rules/sbayesrc.smk` (new) — `rule sbayesrc` (genome-wide; mem_mb=32000, threads=4; ldref base +
  annot default to GenoPred staged paths when config NA; ldref = <base>/<population>).
- `scripts/sbayesrc.R` (new) — build `.ma`, run gctb impute + sbayes RC, parse `.parRes` → clean.csv.
- `scripts/functions/package_results_functions.R` — new `process_sbayesrc()` (gated).
- `scripts/package_results.R` — `gwas_qc$sbayesrc_dat <- process_sbayesrc(...)`.
- `tests/config_test.yaml`, `prepare_resources.yaml` — `sbayesrc: T`.

## Files touched (viewer)
- `R/utils_config.R` — `sbayesrc <- flag(config, "sbayesrc")` + added to returned flag list.
- `app.R` — `toggle("SNP-h² & rG", any(cf$ldsc, cf$ldsc_gencor, cf$sbayesrc))` (tab shows when
  sbayesrc ran even if ldsc did not).
- `R/mod_h2_gencor.R` — new "SBayesRC (h² & polygenicity)" sub-tab in the h2_gencor_tabs
  tabsetPanel + `output$sbayesrc_table` / `output$sbayesrc_legend` (reads gwas_qc$sbayesrc_dat$val,
  gated on cf$sbayesrc; shows SNP-h² (SE), non-zero SNPs, polygenicity π).
All edited R files parse clean.

## End-to-end test (HPC via SSH, as the website does)
- Dev VM can't execute the pipeline: the genodisc conda env binaries' shebangs point at /cephfs
  (absent here); snakemake's env-active guard (dependencies.smk:23) also blocks a local run. So the
  run goes to the HPC — SSH er_slurm_group1967@hpc.create.kcl.ac.uk (key /home/ollie/.ssh/slurm_svc_acct),
  partition neurohack_cpu (same mechanism the Django app uses in modules/hpc_interface.py).
- 2026-10-19: sbatch'd the sbayesrc test target on a compute node (4 cpu, 32 GB):
  `snakemake tests/output/results/COAD01chr22/sbayesrc/COAD01chr22.GW.sbayesrc.clean.csv
   --configfile tests/config_test.yaml --use-conda -c4` → SLURM job 38091704. COMPLETED but
  snakemake exited 1 at config-parse on an UNRELATED pre-existing validation:
  "magma_interaction.enabled is true but pathway_gmt_dir=resources/data/msigdb contains no *.gmt
  files" (this checkout's resources/ lacks the MSigDB gmts; nothing to do with SBayesRC). It aborts
  all targets. → Added `tests/config_test_sbayesrc.yaml`: a minimal overlay enabling ONLY sbayesrc
  (base config.yaml defaults everything else off + magma_interaction disabled). Re-ran via
  sbatch → SLURM job 38091737.
- 2026-10-19: Iterated on test isolation (jobs 38091737/775/825/854/938). Findings:
  * `repo/current` is a SYMLINK to `GenoDisc_v1` (that IS the dev worktree — edits land correctly).
  * The bare checkout's `resources/` is unstaged and `tests/output` is partial, so snakemake wants
    to rebuild unrelated upstream (prep_1kg → sumstat_prep_i), which fails (no 1kg data here) and —
    critically — snakemake DELETED the pre-staged `cleaned.gz` when sumstat_prep_i failed. None of
    this touches SBayesRC.
  * Resolution: restored `cleaned.gz` (real COAD01chr22 data, reconstructed from the earlier .ma)
    and ran ONLY the leaf rule: `snakemake <clean.csv target> --configfile tests/config_test_sbayesrc.yaml
    --use-conda -c4 --allowed-rules sbayesrc download_gctb --rerun-triggers mtime`.
- 2026-10-19: **END-TO-END SUCCESS (SLURM job 38091938, compute node, SNAKEMAKE_EXIT=0, 1/1 steps).**
  `download_gctb` staged gctb 2.5.5; `rule sbayesrc` ran the two-step gctb recipe + parser and wrote
  `tests/output/results/COAD01chr22/sbayesrc/COAD01chr22.GW.sbayesrc.clean.csv`:
  `gwas,h2_est,h2_se,n_nonzero,n_snp_total,pi_est = COAD01chr22,7.6e-05,0.00024,112.99,16774,0.00674`
  (h2≈0 / ~113 non-zero SNPs — expected for chr22-only smoke data; `.parRes` retained alongside).
  `process_sbayesrc()` verified to read this into `gwas_qc$sbayesrc_dat$val` (h2_est/h2_se/n_nonzero/
  n_snp_total/pi_est). Viewer sub-tab reads the same path. All edited R parses clean.

## Verified
- download_gctb rule (gctb 2.5.5 official binary) — runs on HPC, reads GenoPred eigen refs + annot.
- sbayesrc rule — runs end-to-end via snakemake on COAD01chr22, produces clean.csv.
- process_sbayesrc() — packages clean.csv into gwas_qc$sbayesrc_dat$val.
- Viewer (utils_config flag, app.R gate, mod_h2_gencor sub-tab) — parses; mirrors working LDSC path.

## Caveats / notes for reviewer
- chr22-only test is a MECHANICS smoke test (h2≈0 is expected); meaningful h²/polygenicity needs a
  whole-genome GWAS (e.g. config_test_full) — SBayesRC is genome-wide.
- The real pipeline produces cleaned.gz via sumstat_prep_i normally; the isolated test restored it
  because this bare checkout's upstream resources aren't staged (pre-existing env limitation).
- A full-bundle browser render of the new viewer sub-tab is the remaining manual check (same as
  prior viewer changes); the data path is validated.
- Reference data shared with GenoPred (eigen LD + annot) via hardcoded GenoPred paths when config
  sbayesrc_ldref/annot are NA — confirm that cross-pipeline path coupling is acceptable, or set
  explicit paths in prod config.
- NOT COMMITTED — pending user review of pipeline changes (per instruction).

## Docs drafted (from this log)
- `CONTRIBUTING.md` (repo root, new) — recipe for adding an optional analysis, with the SBayesRC
  worked example.
- `AGENTS.md` (repo root, new) — condensed AI checklist + invariants; points to CONTRIBUTING.md.
- `CLAUDE.md` (repo root, existing generic guidelines) — appended a "Project recipes" pointer to
  CONTRIBUTING.md / AGENTS.md and the HPC/SSH test note.
- Added `tests/config_test_sbayesrc.yaml` (isolated SBayesRC-only test config).

## Post-review changes (user, 2026-10-19)
- Explicit ref paths now live in `config.yaml` (`sbayesrc_ldref`/`sbayesrc_annot` set to GenoPred's
  staged dirs, no longer NA); `rules/sbayesrc.smk` reads them directly (removed the hardcoded
  GenoPred fallback + NA logic).
- Added `tests/config_sbayesrc_wg.yaml` (whole-genome COAD01, SBayesRC + LDSC, resdir=pipeline_resources)
  and launched a REAL whole-genome end-to-end run to a packaged bundle for the Shiny app:
  SLURM job 38092139 (8 cpu, 48 GB) → target tests/output_sbrc_wg/results/bundle.tar.gz.
- 2026-10-19: **WHOLE-GENOME END-TO-END SUCCESS** (job 38092139, ~1h36m, 7/7 steps, SNAKEMAKE_EXIT=0).
  Bundle: `tests/output_sbrc_wg/results/bundle.tar.gz` (+ manifest.json), `gwas_qc.rds` contains
  sbayesrc_dat. Real genome-wide values (NOT the chr22 ~0 smoke):
  SBayesRC h2_est=0.0322 (SE 0.0014), n_nonzero=10208.6 of n_snp_total=1,154,502, pi_est=0.00884;
  LDSC h2=0.0664 (0.0044) for comparison (SBayesRC lower — expected; HM3 + annotation model).
- 2026-10-19: CONFIG RECONCILIATION — a parallel "lost-signal reporting" change landed in
  repo/current and **reverted my config.yaml explicit-paths edit back to `sbayesrc_ldref/annot: NA`**
  (my report.smk/sumstat_qc.smk/package_results.R sbayesrc edits survived). Per the "intentional,
  don't revert" note I left config.yaml as NA and instead **restored the NA→GenoPred fallback in
  rules/sbayesrc.smk** so the pipeline works with NA and still honours explicit paths if set.
  NOTE FOR USER: your requested explicit-paths config edit was overwritten by that parallel change —
  re-apply in config.yaml if you still want explicit paths rather than the NA fallback.
- NOTE: repo/current is a SHARED working tree now carrying BOTH this SBayesRC work and an unrelated
  uncommitted lost-signal feature (config.yaml removed_p_thresh; package_results.R lost_signals/
  manhattan_raw; GenoUtils SHA d734c605). A commit must stage ONLY the SBayesRC files.

## Presentation of SNP-h2 (J. Zeng advice, 2026-10-19)
- Posterior SE = SD of MCMC samples (do NOT divide by sqrt(N)) — our h2_sd is already that.
- SNP-h2 is bounded at 0 and its posterior is skewed near zero → report posterior MODE + 95% HPD,
  plus a posterior-density plot. **Discard burn-in before computing mode/HPD.**
- gctb `--write-mcmc-txt` was flaky and dumps ~4 GB per-SNP files, so instead run
  `--sbayes RC ... --out-freq 1`, capture gctb stdout, and parse the per-iteration `Iter...hsq`
  table (iter=$1, hsq=$13; chain 3000, burn-in 1000 from the log) → 2000 post-burn-in hsq draws.
  Compute mode=peak of density(.,from=0), 95% HPD (narrowest interval), mean, SD; persist the draws
  to `{gwas}.GW.sbayesrc.h2_samples.txt`.
- clean.csv now: gwas,h2_mode,h2_mean,h2_sd,h2_hpd_low,h2_hpd_high,n_post_samples,n_nonzero,
  n_snp_total,pi_est. process_sbayesrc() carries val + h2_samples. Viewer sub-tab shows mode / 95%
  HPD / posterior mean (SD) / polygenicity + a ggplot density (mode solid, HPD dashed).
- config.yaml explicit paths re-applied (user confirmed); fallback kept in rules/sbayesrc.smk.
- Re-ran WG (job 38094363) to regenerate the bundle with the new outputs (only sbayesrc + package
  rerun; resdir=pipeline_resources fully staged).
