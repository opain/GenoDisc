> **AUTOMATICALLY GENERATED DRAFT — CHECK AND EDIT BEFORE USE.**
> This methods text was produced by the GenoDisc write-up generator (v0.1.0) from the recorded run configuration, software provenance and result objects (results schema v2, pipeline 1.4.0). It describes what the pipeline actually ran. It is a starting point, not a finished section: verify every number, version and citation, fill the `[CITE: ...]` placeholders, and edit for style before use. Citation keys (`@key`) refer to `references.bib` in this folder.

# Methods

Summary statistics were quality-controlled and harmonised with the GenoDisc cleaning pipeline against the 1000 Genomes + HGDP reference panel (1000 Genomes Phase 3 + HGDP) [@1000genomes2015]. Variants were removed if the imputation INFO score was below 0.9, if the reported or reference minor-allele frequency (MAF) was below 0.01, or if the reported and reference MAF differed by more than 0.2; alleles were harmonised to the reference and strand-ambiguous and duplicate variants resolved. The genome build was identified as GRCh37. Of 9,000,000 variants in the input, 7,000,000 remained after QC.

SNP-based heritability on the observed scale and the LD-score regression intercept were estimated with LD Score Regression (LDSC) (git commit aa33296) [@bulik2015ldsc], using the European-ancestry HapMap3 LD scores supplied with the pipeline.

Independent genome-wide significant loci were identified by LD-based clumping in PLINK (v1.90b6.21) [@chang2015plink] using the 1000 Genomes + HGDP reference panel (1000 Genomes Phase 3 + HGDP) [@1000genomes2015], matched to the GWAS population (index-variant P < 1.0e-05, r2 < 0.1, 500 kb window). The nearest gene to each index variant was annotated.

Gene-level association was tested with MAGMA (v1.10) [@deleeuw2015magma], mapping variants to genes with a 35 kb upstream / 10 kb downstream window and aggregating variant P-values accounting for LD. Gene P-values were FDR-corrected.

Transcriptome-wide association analysis (TWAS) was performed with FUSION (git commit 4635dd1) [@gusev2016fusion] across 2 expression SNP-weight panels (Blood (GTEx) and Brain - Cortex (GTEx); 9,000 gene–tissue models) using GTEx v8 FUSION expression SNP-weights (GTEx v8) [@gtex2020]. Gene-level associations were corrected for multiple testing across gene panels using the Benjamini–Hochberg procedure. Colocalisation of each significant association was assessed with coloc (v3.1) [@giambartolomei2014coloc].

All analyses were orchestrated by the GenoDisc post-GWAS pipeline (version 1.4.0) and are available through the GenoDisc web server [CITE: GenoDisc paper]. Results were explored with the GenoDisc Output Viewer [CITE: GenoDisc paper].

## Reconciliation — check before use

The following analyses were requested in the configuration but produced no results for any submitted GWAS, so they are not described above. This is expected when an analysis is restricted to certain ancestries (e.g. European-only steps run on a non-European GWAS); otherwise check the pipeline logs for failures:
- `cojo`
