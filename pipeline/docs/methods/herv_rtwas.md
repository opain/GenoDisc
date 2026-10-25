# HERV retrotranscriptome TWAS panel (`CMC.BRAIN.HERV`)

FUSION expression SNP-weights for **human endogenous retroviruses (HERVs)** in
adult brain, added as a first-class non-GTEx TWAS panel (retrotranscriptome-wide
association study, rTWAS). Enable with `twas_panel_herv: T`.

| | |
|---|---|
| Tissue | Adult dorsolateral prefrontal cortex (CommonMind Consortium) |
| Sample size (N) | 563 (EUR) |
| Ancestry | EUR (AFR weights not included) |
| Feature type | **HERVs only** (1238) — `FAMILY_cytoband` ids, e.g. `ERVLE_1p36.32a` |
| Features in archive | 8194 = 6956 genes + 1238 HERVs; the gene features are dropped (they duplicate the existing `CMC.BRAIN.RNASEQ` DLPFC panel) |
| Genome build | hg19 / GRCh37, rsID-keyed |
| LD reference used | GenoDisc **dense** EUR (`1KG.Phase3.EUR.MAF_001.fusion`) — the weights were trained on a full 1000G reference, so the HapMap3 reference would drop ~80% of weight SNPs; the dense reference retains ~99% (verified; ~1% attrition vs the paper's bundled LD reference, which is therefore not shipped) |
| Multiple testing | **Separate FDR** from the gene panels (HERVs are not genes) |
| Downstream | **Excluded** from gene-based steps (TWAS-GSEA, DrugTargetor, CMAP, pathway) — HERVs have no gene-set membership. A parse-time error is raised if a gene-based analysis is requested with the HERV panel as the only TWAS panel. |

**Source.** figshare [22179655](https://doi.org/10.18742/22179655), file
`FUSION_weights_EUR_all_1.tar.gz` (md5 `e5d1c3135a8a56d9ceb3a1ad2558810c`),
downloaded + verified by rule `prep_herv_weights` and formatted by
`scripts/prep_herv.R` (which flattens the nested archive path, subsets to HERV
features, and drops the gene features). AFR weights and the FOCUS weights are not
used.

**Citation.** Duarte RRR, Pain O, et al. (2024). *Integrating human endogenous
retroviruses into transcriptome-wide association studies highlights novel risk
factors for major psychiatric conditions.* Nature Communications.
doi:10.1038/s41467-024-48153-z. Tutorial: https://rodrigoduarte88.github.io/neuro_rTWAS/

**Licence.** CC BY 4.0 (redistribution and derived use permitted with attribution).

**Display.** HERV features are surfaced in the viewer by their HERV id (the
`Gene Symbol`/`Ensembl ID` columns are populated with the HERV id during
packaging, so tables, the summary heatmap and cross-GWAS views show them without
viewer-code changes). The panel appears as "Brain - DLPFC HERVs (Duarte)".
