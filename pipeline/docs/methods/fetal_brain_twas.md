# Fetal-brain gene-level TWAS panel (`FETAL.BRAIN.RNASEQ`)

FUSION expression SNP-weights for **human fetal brain**, added as a first-class
non-GTEx TWAS panel. Enable with `twas_panel_fetal_brain: T`.

| | |
|---|---|
| Tissue / age | Fetal brain, 12–19 post-conception weeks |
| Sample size (N) | 120 fetal brains |
| Ancestry | EUR-trained |
| Feature type | Genes (Ensembl gene IDs) |
| Features in panel | 1351 cis-heritable (P<0.01) genes |
| Genome build | hg19 / GRCh37, rsID-keyed |
| LD reference used | GenoDisc HapMap3-restricted EUR (`1KG.Phase3.EUR.MAF_001.hm3`) — the weights are HapMap3-based, so ~100% of weight SNPs are retained (verified) |
| Multiple testing | Pooled with the other gene panels (standard cross-panel FDR) |
| Downstream | Flows through TWAS-GSEA / DrugTargetor like any gene panel |

**Source.** figshare [11637036](https://doi.org/10.6084/m9.figshare.11637036),
file `op-fusion-fetal_brain-full_sample-gene_level.tar.gz`
(md5 `e16e9de8f10be07eadd5c3f20479868a`), downloaded + verified by rule
`prep_fetal_brain_weights` and formatted by `scripts/prep_fetal_brain.R`. The
transcript-level file is deliberately not used.

**Citation.** Bray N, Pain O (2020). *Human fetal brain expression SNP-weights
for use in TWAS.* figshare. doi:10.6084/m9.figshare.11637036. Underlying fetal
brain eQTL data: O'Brien/Hall et al.

**Licence.** CC BY 4.0 (redistribution and derived use permitted with attribution).

**Notes.** A handful of Ensembl IDs may not map to a current gene symbol in
`biomart_genes_grch37.tsv`; those features are still tested and shown by Ensembl
ID (position-based), they just carry no symbol. Feature counts actually *tested*
for a given GWAS depend on GWAS×reference SNP overlap (reported at verification).
