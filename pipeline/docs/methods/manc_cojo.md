# GenoDisc method contribution — Manc-COJO

> Status: Manc-COJO is **the COJO method used** — the website's "COJO" toggle runs it,
> and its results populate the existing COJO surface in the results package and output
> viewer. It was validated to reproduce GCTA-COJO's single-ancestry selection exactly
> (identical SNPs and joint effects on test data) before the swap. GCTA-COJO (config
> `cojo: T`) is retained for internal comparison but is no longer requested by the website.

---

## 1. Method identity
- **Display name:** Manc-COJO
- **One-line summary:** Identifies independently associated SNPs by stepwise conditional & joint analysis of GWAS summary statistics — reproduces GCTA-COJO's single-ancestry selection but much faster.
- **Config flag:** `manc_cojo`
- **Snakemake rule(s):** `manc_cojo`, `manc_cojo_all_chr`, `process_manc_cojo` (the last reuses `scripts/process_cojo.R --method manc`)
- **Rule module:** rules live in `pipeline/rules/sumstat_qc.smk` (alongside the GCTA `cojo` rules), not a separate `.smk`.
- **External tool + version:** `manc_cojo` v1.1.0 — official statically-linked Linux release binary, staged by `rule download_manc_cojo` to `{resdir}/software/manc_cojo/manc_cojo_v1.1.0` (sha256-pinned, verified on download).
- **References / citations:** Wang et al. 2026, *Nature Genetics*, doi:[10.1038/s41588-026-02799-1](https://doi.org/10.1038/s41588-026-02799-1). Code: https://github.com/light156/multi-ancestry-COJO ; docs: https://light156.github.io/multi-ancestry-COJO-docs/
- **Tool licence:** GPL-3.0
- **Reference-data licence / redistribution:** none new — reuses the same 1KG+HGDP per-population LD reference already staged for GCTA-COJO/clumping.
- **Citation users should give:** Wang et al. 2026 (Manc-COJO), in addition to GCTA-COJO (Yang et al. 2012) while both are reported.
- **Maintainer / contact:** Oliver Pain (oliver.pain@kcl.ac.uk)

## 2. Config parameters
| Key | Type | Default | Description | Example |
|---|---|---|---|---|
| `manc_cojo` | T/F | `F` | Run Manc-COJO stepwise selection alongside GCTA-COJO (independent of `cojo`). | `T` |

No other new keys: the reference, selection p-value (5e-8), MAF filter (0.01), window (10000 kb) and collinearity (0.9) are the same as the GCTA-COJO rule and are not separately configurable. The config parser now also rejects any `cojo`/`manc_cojo` value that is not strictly `T`/`F`.

## 3. Submission-panel metadata (job-submission analysis toggle)
**Not yet surfaced.** Deliberately out of scope for this change — Manc-COJO runs only when `manc_cojo: T` is set directly in the pipeline config. No `ANALYSIS_LABELS` / form card / `ANALYSIS_CONFIG_KEYS` entry is requested until the comparison concludes.
- **Version-gating config key** (when later added): `manc_cojo`
- **Dependencies:** none (independent of `cojo`)
- **Supported populations / constraints:** same as GCTA-COJO — any population with a staged 1KG+HGDP reference (EUR, …). The reference is picked by the per-GWAS `population` column.

## 4. Compute cost (for the submission-form resource estimate)
- **Snakemake rule name (for RULE_COST key):** `manc_cojo`
- **Typical wall-time (per GWAS):** seconds per chromosome (chr22 test: ~2.8 s). Dramatically faster than GCTA-COJO (upstream report: >100×).
- **Peak memory:** <0.1 GB on the chr22 test; scales with per-chromosome SNP count.
- **Per-chromosome or genome-wide:** per_chr True (one job per chromosome, same as GCTA-COJO).
- **Threads:** 1 (`--thread-num {threads}`). Per-chromosome cost is I/O-bound (reading the reference `.bed`/`.bim`), so extra threads give no speed-up (measured 1 thread ~82 s vs 4 threads ~99 s over 22 chr).

## 5. Outputs
- **Final per-GWAS target file(s):** `{outdir}/results/{gwas}/manc_cojo/{gwas}.GW.manc_cojo.clean.csv` and `{gwas}.GW.manc_cojo.status.csv` (same schema as the GCTA `cojo/` equivalents). Per-chromosome `{gwas}_chr{chr}.jma.cojo`, `.log`, and a `{gwas}_chr{chr}.manc_cojo.status` side-car are also written.
- **Packaged bundle location:** `snp_assoc$cojo` + `snp_assoc$cojo_status` (the existing COJO slot). `package_results.R` sources these from `manc_cojo/` when `manc_cojo: T`, else `cojo/` — identical schema, so no reader/viewer schema change. (Nested in the existing `snp_assoc` block → no `reader.R` change.)
- **Kept vs discarded:** the `.jma.cojo` (selected SNPs + joint effects) and `.log` (incl. peak-memory line) are kept for the comparison; `--output-all` extras (`.cma.cojo`/`.ldr.cojo`/`.badsnps`) are not requested.

## 6. Analyses-catalogue entry (website "choosing analyses" table)
- **Category:** SNP Associations
- **"What it does":** Conditional & joint analysis (Manc-COJO) — identifies conditionally independent association signals within a locus. (Surfaced via the existing "COJO conditional analysis" toggle, whose internal field name stays `cojo`; the generated config emits `manc_cojo: T`.)

## 7. Viewer presentation (website "Interpreting your results" page)
Shown under the existing **COJO** surface in the output viewer (Lead-variants tab / COJO method option) — the viewer's `cf$cojo` is true when either engine ran, and Manc populates the same `snp_assoc$cojo` object, so the presentation is unchanged. No new viewer tab or figure.

## 8. Technical method write-up (pipeline docs site)
Manc-COJO v1.1.0 (Wang et al. 2026) reproduces GCTA-COJO's single-ancestry stepwise conditional selection of independently associated variants, but with a much faster implementation. When `manc_cojo: T`, GenoDisc runs it alongside GCTA-COJO on the same cleaned summary statistics (`{gwas}.cleaned.cojo`, header `SNP A1 A2 freq b se p N`), the same per-population 1KG+HGDP LD reference, and the same parameters (selection p < 5e-8, MAF ≥ 0.01, 10 Mb window, collinearity 0.9). Results are written to `{outdir}/results/{gwas}/manc_cojo/` and packaged into the existing COJO slot. Small, legitimate differences from GCTA-COJO are possible in principle (e.g. exclusion of multi-allelic/monomorphic SNPs), though on the test GWAS the two selected an identical set of SNPs with identical joint effects. Unlike GCTA-COJO, Manc-COJO does not abort a chromosome when the number of selected signals reaches the LD reference sample size — it continues and resolves collinearity, so such chromosomes are flagged (`exceeds_reference_n`) but their results are retained. Manc-COJO is now the COJO method the website requests; GCTA-COJO (`cojo: T`) is retained for internal comparison.

## 9. Testing / verification
- **Test config + GWAS used:** `tests/config_test_full.yaml` on `COAD01` (whole genome, EUR); `tests/config_test.yaml` on `COAD01chr22` for dry-run wiring.
- **How it was run end to end:** `snakemake <...>/manc_cojo/COAD01.GW.manc_cojo.clean.csv --configfile tests/config_test_full.yaml --use-conda -c1` on CREATE; gating `./manc_cojo --help` and a chr22 selection confirmed on the HPC node. `snakemake -n` confirms GCTA-only, Manc-only and both-enabled configs each plan correctly, and an invalid `manc_cojo` value fails at parse time.
- **Expected / observed result:** whole-genome COAD01 comparison (`misc/cojo_comparison/compare_cojo.R`), both methods freshly run on identical input at `--cojo-p 5e-8`: GCTA-COJO and Manc-COJO each select **41 SNPs, all 41 shared (Jaccard 1.00)**, with allele-aligned shared-SNP joint-effect (bJ) correlation = **1.00** — the two methods agree exactly on this dataset. Timing: Manc ~5.3× faster (95.7 s vs 504.3 s at 4 threads vs GCTA 1), and ~6.3× at 1 thread each (82 s vs 516 s); Manc's per-chromosome cost is I/O-bound (reading the `.bed`), not compute. Peak memory ~1/7th of GCTA (173 MB vs ~1.2 GB). GCTA output is byte-identical before/after the `process_cojo.R --method` refactor (verified by diff against the pre-change script).
  - *Note:* an initial comparison showed a spurious GCTA 58 vs Manc 41 with 17 GCTA-only SNPs on 6 chromosomes. Root cause was **stale `.jma.cojo` from a prior run** (GCTA's `.jma.cojo` was not a declared Snakemake output, so when a re-run selected nothing on a chromosome the old file persisted and `process_cojo.R` read it). Fixed: both the `cojo` and `manc_cojo` rules now delete their per-chromosome `.jma.cojo`/`.cma.cojo`/`.ldr.cojo`/`.badsnps` at the start of the shell, so a shrinking re-run cannot leave orphaned selections.

## 10. Known failure modes
- **Output folder must pre-exist:** Manc-COJO errors if the `--out` directory is absent; the rule `mkdir -p`s it first.
- **`--maf` default differs:** Manc-COJO's `--maf` defaults to 0 (GCTA applies 0.01); the rule passes `--maf 0.01` explicitly to match.
- **No "too many SNPs" abort:** GCTA-COJO's abort (→ `reference_too_small`, chromosome dropped) does not occur; Manc-COJO keeps going. The rule instead counts selected SNPs vs the reference `.fam` size and flags `exceeds_reference_n` while keeping results.
- **SNPs absent from the LD reference:** multi-allelic (shared-ID) and monomorphic reference SNPs are excluded, so a few variants GCTA keeps may be dropped (expected small divergence).
- **Binary portability:** the pinned release binary is statically linked x86-64; the sha256 is verified on download so a silent upstream re-release fails loudly.
