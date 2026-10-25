#!/usr/bin/Rscript

# Methods-paragraph builders. Each takes (ctx, g, meta) where `g` is the collected
# output[[gwas]] and `meta` is the gwas_list row; returns a Markdown paragraph for a
# block that RAN, or NULL for a block that did not (never "not performed"). Builders
# record the block id in ctx$ran_blocks and cite via writeup_tool()/writeup_ref().
#
# No quantitative value is literal here: numbers come from `g`, versions/citations
# from the registry, and fixed algorithm constants from writeup/params.yaml.
#
# gstr() (defined in format.R) is glue()->character; used instead of the R>=4.1
# native pipe because the pipeline runs R 4.0.2.

library(glue)

# mark a block as having run (so reconciliation can detect requested-but-absent)
.ran <- function(ctx, id) ctx$ran_blocks <- union(ctx$ran_blocks, id)

.is_cross_ancestry <- function(g) isTRUE(pluck(g, "gwas_qc", "cross_ancestry_weights"))

# ---- preprocessing (always runs) --------------------------------------------

methods_preprocessing <- function(ctx, g, meta) {
  v <- pluck(g, "gwas_qc", "cleaner_dat", "val")
  .ran(ctx, "preprocessing")
  opt <- ctx$opt
  # val$build is the list returned by extract_build(): list(build=, overlap=).
  bq <- if (!is.null(v)) v$build else NULL
  build <- if (is.list(bq)) bq$build else bq
  if (is.null(build) || length(build) != 1 || is.na(build)) build <- NA
  ref <- writeup_tool(ctx, "ref_1kg_hgdp")
  n_orig  <- if (!is.null(v)) v$n_var_orig else NA
  n_final <- if (!is.null(v)) v$n_snp_final else NA
  info <- cfg_val(ctx, "info", 0.9); maf <- cfg_val(ctx, "maf", 0.01); maf_diff <- cfg_val(ctx, "maf_diff", 0.2)
  build_clause <- if (!is.na(build)) gstr(" The genome build was identified as {build}.") else ""
  count_clause <- if (!is.na(n_orig) && !is.na(n_final)) {
    gstr(" Of {fmt_int(n_orig)} variants in the input, {fmt_int(n_final)} remained after QC.")
  } else if (!is.na(n_final)) {
    gstr(" {fmt_int(n_final)} variants remained after QC.")
  } else ""
  gstr(
    "Summary statistics were quality-controlled and harmonised with the GenoDisc ",
    "cleaning pipeline against the {ref}. Variants were removed if the ",
    "imputation INFO score was below {fmt_num(info, opt)}, if the reported or reference ",
    "minor-allele frequency (MAF) was below {fmt_num(maf, opt)}, or if the reported and ",
    "reference MAF differed by more than {fmt_num(maf_diff, opt)}; alleles were harmonised to ",
    "the reference and strand-ambiguous and duplicate variants resolved.",
    "{build_clause}{count_clause}"
  )
}

# ---- heritability: LDSC -----------------------------------------------------

methods_ldsc_h2 <- function(ctx, g, meta) {
  v <- pluck(g, "gwas_qc", "ldsc_dat", "val")
  if (is.null(v)) return(NULL)
  .ran(ctx, "ldsc_h2")
  tool <- writeup_tool(ctx, "ldsc")
  gstr(
    "SNP-based heritability on the observed scale and the LD-score regression ",
    "intercept were estimated with {tool}, using the European-ancestry HapMap3 ",
    "LD scores supplied with the pipeline."
  )
}

methods_sbayesrc <- function(ctx, g, meta) {
  v <- pluck(g, "gwas_qc", "sbayesrc_dat", "val")
  if (is.null(v)) return(NULL)
  .ran(ctx, "sbayesrc")
  tool <- writeup_tool(ctx, "sbayesrc")
  nchain <- cfg_val(ctx, "sbayesrc_chain_length", NA); burn <- cfg_val(ctx, "sbayesrc_burn_in", NA)
  mcmc <- if (!is.na(nchain) && !is.na(burn)) gstr(" (MCMC chain length {fmt_int(nchain)}, burn-in {fmt_int(burn)})") else ""
  gstr(
    "Genome-wide SNP-based heritability and polygenicity (the number of variants ",
    "with non-zero effect) were estimated with {tool}, a Bayesian mixed model ",
    "incorporating functional genomic annotations{mcmc}."
  )
}

methods_sbayesrc_r <- function(ctx, g, meta) {
  v <- pluck(g, "gwas_qc", "sbayesrc_r_dat", "val")
  if (is.null(v)) return(NULL)
  .ran(ctx, "sbayesrc_r")
  tool <- writeup_tool(ctx, "sbayesrc_r")
  gstr(
    "Genome-wide SNP-based heritability, polygenicity and per-annotation ",
    "heritability enrichment were additionally estimated with the {tool}."
  )
}

methods_sumher <- function(ctx, g, meta) {
  v <- pluck(g, "gwas_qc", "sumher_dat", "val")
  if (is.null(v)) return(NULL)
  .ran(ctx, "sumher")
  tool <- writeup_tool(ctx, "sumher")
  cutoff <- cfg_val(ctx, "sumher_cutoff", NA)
  cut_clause <- if (!is.na(cutoff)) gstr(" with large-effect predictors excluded at a per-SNP heritability cutoff of {fmt_num(cutoff, ctx$opt)}") else ""
  gstr(
    "SNP-based heritability was also estimated with {tool} under the BLD-LDAK ",
    "heritability model from pre-computed HapMap3 tagging files{cut_clause}. The three ",
    "heritability methods assume different models and are not expected to agree."
  )
}

methods_gencor <- function(ctx, g, meta) {
  sec    <- pluck(g, "gwas_qc", "ldsc_gencor_dat", "table")
  within <- pluck(g, "gwas_qc", "ldsc_gencor_within_dat", "table")
  if (is.null(sec) && is.null(within)) return(NULL)
  .ran(ctx, "gencor")
  tool <- writeup_tool(ctx, "ldsc_rg")
  parts <- character(0)
  if (!is.null(sec))    parts <- c(parts, "a curated set of reference GWAS")
  if (!is.null(within)) parts <- c(parts, "each pair of submitted European-ancestry GWAS")
  gstr(
    "Genetic correlations were estimated using bivariate (cross-trait) ",
    "{tool} between the submitted summary statistics and {oxford(parts)}. ",
    "Correlations were corrected for multiple testing across all tested pairs using the ",
    "Benjamini–Hochberg procedure."
  )
}

# ---- SNP associations -------------------------------------------------------

methods_clump <- function(ctx, g, meta) {
  d <- pluck(g, "snp_assoc", "clump")
  if (is.null(d)) return(NULL)
  .ran(ctx, "clump")
  tool <- writeup_tool(ctx, "plink")
  ref  <- writeup_tool(ctx, "ref_1kg_hgdp")
  p1 <- param(ctx, "clump", "p1"); r2 <- param(ctx, "clump", "r2"); kb <- param(ctx, "clump", "kb")
  gstr(
    "Independent genome-wide significant loci were identified by LD-based clumping ",
    "in {tool} using the {ref}, matched to the GWAS population ",
    "(index-variant P < {fmt_p(p1, ctx$opt)}, r2 < {fmt_num(r2, ctx$opt)}, ",
    "{fmt_int(kb)} kb window). The nearest gene to each index variant was annotated."
  )
}

methods_cojo <- function(ctx, g, meta) {
  d <- pluck(g, "snp_assoc", "cojo")
  if (is.null(d)) return(NULL)
  .ran(ctx, "cojo")
  key <- if (cfg_on(ctx, "manc_cojo")) "manc_cojo" else "gcta_cojo"
  tool <- writeup_tool(ctx, key)
  ref  <- writeup_tool(ctx, "ref_1kg_hgdp")
  p <- param(ctx, "cojo", "p")
  gstr(
    "Conditionally independent association signals were identified by stepwise ",
    "conditional and joint (COJO) analysis with {tool}, using the {ref}, ",
    "matched to the GWAS population (selection threshold P < {fmt_p(p, ctx$opt)})."
  )
}

methods_finemap <- function(ctx, g, meta) {
  d <- pluck(g, "snp_assoc", "susie")
  if (is.null(d)) return(NULL)
  .ran(ctx, "finemap")
  tool <- writeup_tool(ctx, "susie")
  Lp <- param(ctx, "susie", "L_primary"); Ls <- param(ctx, "susie", "L_secondary")
  mhc <- param(ctx, "susie", "mhc_note")
  gstr(
    "Statistical fine-mapping of genome-wide significant loci was performed with ",
    "{tool}, allowing up to {fmt_int(Lp)} causal variants per locus (with a ",
    "single-effect pass allowing {fmt_int(Ls)}); {mhc}. Credible sets were reported ",
    "at the 95% level."
  )
}

# ---- molecular associations -------------------------------------------------

.twas_panel_datasets <- function(ctx) {
  ds <- character(0)
  if (cfg_on(ctx, "twas_panel_psychencode")) ds <- c(ds, "psychencode_twas")
  if (length(cfg_list(ctx, "gtex_weights")) > 0) ds <- c(ds, "gtex_v8")
  if (cfg_on(ctx, "twas_panel_fetal_brain")) ds <- c(ds, "fetal_brain_panel")
  if (cfg_on(ctx, "twas_panel_herv")) ds <- c(ds, "herv_panel")
  ds
}

methods_twas <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "exp", "fusion")
  if (is.null(d) || is.null(d$res)) return(NULL)
  .ran(ctx, "twas")
  fusion <- writeup_tool(ctx, "fusion")
  panels <- d$panels
  n_panel <- if (!is.null(panels)) nrow(panels) else NA
  panel_names <- if (!is.null(panels)) oxford(panels$Panel) else "the selected expression panels"
  n_gene <- if (!is.null(panels)) sum(suppressWarnings(as.numeric(panels$N_gene)), na.rm = TRUE) else NA
  ds <- .twas_panel_datasets(ctx)
  ds_cite <- if (length(ds) > 0) paste0(" using ", oxford(vapply(ds, function(k) writeup_tool(ctx, k), character(1)))) else ""
  herv <- if (cfg_on(ctx, "twas_panel_herv")) "; HERV (retrotranscriptome) features were FDR-corrected separately from genes" else ""
  cross <- if (.is_cross_ancestry(g)) " As the GWAS is of non-European ancestry, colocalisation was not assessed (the weights are European-trained)." else gstr(" Colocalisation of each significant association was assessed with {writeup_tool(ctx, 'coloc')}.")
  gstr(
    "Transcriptome-wide association analysis (TWAS) was performed with {fusion} across ",
    "{fmt_int(n_panel)} expression SNP-weight {pluralise('panel', n_panel)} ({panel_names}; ",
    "{fmt_int(n_gene)} gene–tissue models){ds_cite}. Gene-level associations were corrected ",
    "for multiple testing across gene panels using the Benjamini–Hochberg procedure{herv}.{cross}"
  )
}

methods_pwas <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "protein", "fusion")
  if (is.null(d) || is.null(d$results)) return(NULL)
  .ran(ctx, "pwas")
  fusion <- writeup_tool(ctx, "fusion")
  ref <- writeup_tool(ctx, "rosmap_banner")
  panels <- d$panels
  panel_names <- if (!is.null(panels)) oxford(panels$Panel) else "the selected protein panels"
  cross <- if (.is_cross_ancestry(g)) " Colocalisation was not assessed for this non-European GWAS." else gstr(" Colocalisation was assessed with {writeup_tool(ctx, 'coloc')}.")
  gstr(
    "Proteome-wide association analysis (PWAS) was performed with {fusion} using ",
    "{ref} protein SNP-weight panels ({panel_names}). Associations were FDR-corrected ",
    "across panels.{cross}"
  )
}

.smr_exp_datasets <- function(ctx) {
  m <- c(
    smr_expression_panel_psychencode = "psychencode_smr",
    smr_expression_panel_metabrain_basalganglia = "metabrain",
    smr_expression_panel_metabrain_cerebellum = "metabrain",
    smr_expression_panel_metabrain_cortex = "metabrain",
    smr_expression_panel_metabrain_hippocampus = "metabrain",
    smr_expression_panel_metabrain_spinalcord = "metabrain",
    smr_expression_panel_eqtlgen = "eqtlgen"
  )
  unique(unname(m[vapply(names(m), function(f) cfg_on(ctx, f), logical(1))]))
}

methods_smr_exp <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "exp", "smr")
  if (is.null(d) || is.null(d$results)) return(NULL)
  .ran(ctx, "smr_exp")
  tool <- writeup_tool(ctx, "smr")
  ds <- .smr_exp_datasets(ctx)
  ds_cite <- if (length(ds) > 0) paste0(" using ", oxford(vapply(ds, function(k) writeup_tool(ctx, k), character(1))), " cis-eQTL data") else ""
  heidi <- if (.is_cross_ancestry(g)) " The HEIDI heterogeneity test was not applied for this non-European GWAS." else " The HEIDI test was used to distinguish pleiotropy/causality from linkage."
  gstr(
    "Summary-data-based Mendelian randomisation (SMR) of gene expression was performed ",
    "with {tool}{ds_cite}. SMR P-values were FDR-corrected across panels.{heidi}"
  )
}

methods_smr_protein <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "protein", "smr")
  if (is.null(d) || is.null(d$results)) return(NULL)
  .ran(ctx, "smr_protein")
  tool <- writeup_tool(ctx, "smr")
  ref <- writeup_tool(ctx, "rosmap_banner")
  gstr(
    "SMR of protein levels was performed with {tool} using {ref} pQTL data, with ",
    "FDR correction of SMR P-values and the HEIDI heterogeneity test applied as above."
  )
}

methods_magma_gene <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "magma")
  if (is.null(d)) return(NULL)
  .ran(ctx, "magma_gene")
  tool <- writeup_tool(ctx, "magma")
  up <- param(ctx, "magma", "window_upstream_kb"); dn <- param(ctx, "magma", "window_downstream_kb")
  gstr(
    "Gene-level association was tested with {tool}, mapping variants to genes with a ",
    "{fmt_int(up)} kb upstream / {fmt_int(dn)} kb downstream window and aggregating variant ",
    "P-values accounting for LD. Gene P-values were FDR-corrected."
  )
}

# ---- drug / gene-set / pathway enrichment -----------------------------------

methods_enrichment <- function(ctx, g, meta) {
  magma_drug <- pluck(g, "tx", "drug", "magma")
  gsea_drug  <- pluck(g, "tx", "drug", "twas_gsea")
  gcsc_drug  <- pluck(g, "tx", "drug", "gcsc")
  cmap       <- pluck(g, "tx", "cmap", "drug")
  if (is.null(magma_drug) && is.null(gsea_drug) && is.null(gcsc_drug) && is.null(cmap)) return(NULL)
  .ran(ctx, "enrichment")
  sentences <- character(0)
  # Cite the Drug Targetor database only within the sentences that actually use it
  # (MAGMA/TWAS-GSEA drug enrichment), not GCSC/CMAP.
  if (!is.null(magma_drug)) {
    sentences <- c(sentences, gstr(
      "Competitive drug-set enrichment was tested with {writeup_tool(ctx, 'magma')} ",
      "over the {writeup_tool(ctx, 'drugtargetor')}."))
  }
  if (!is.null(gsea_drug)) {
    sentences <- c(sentences, gstr(
      "Directional drug-set enrichment of TWAS associations was tested with ",
      "{writeup_tool(ctx, 'twas_gsea')} over the {writeup_tool(ctx, 'drugtargetor')}."))
  }
  if (!is.null(gcsc_drug)) {
    sentences <- c(sentences, gstr(
      "Gene co-regulation score regression drug enrichment was computed with ",
      "{writeup_tool(ctx, 'gcsc')}."))
  }
  if (!is.null(cmap)) {
    sentences <- c(sentences, gstr(
      "TWAS associations were also tested against {writeup_tool(ctx, 'cmap')} drug signatures."))
  }
  # ATC-class model-based enrichment (GLS / VIF-OLS), if run.
  if (!is.null(pluck(g, "tx", "atc_gls")) || !is.null(pluck(g, "tx", "atc_gls_magma")) ||
      !is.null(pluck(g, "tx", "atc_vif")) || !is.null(pluck(g, "tx", "atc_vif_magma"))) {
    sentences <- c(sentences, gstr(
      "Model-based ATC-class enrichment was additionally tested using correlation-aware ",
      "VIF-OLS ({writeup_tool(ctx, 'atc_enrichment_vif')}) and GLS ",
      "({writeup_tool(ctx, 'atc_enrichment_gls')}) contrasts, with FDR correction applied ",
      "within each ATC level."))
  }
  paste0("Drug-repurposing enrichment analyses were performed as follows. ",
         paste(sentences, collapse = " "),
         " All enrichment P-values were FDR-corrected.")
}

methods_pathway <- function(ctx, g, meta) {
  mg <- pluck(g, "tx", "pathway", "magma")
  gs <- pluck(g, "tx", "pathway", "twas_gsea")
  if (is.null(mg) && is.null(gs)) return(NULL)
  .ran(ctx, "pathway")
  tools <- character(0)
  if (!is.null(mg)) tools <- c(tools, writeup_tool(ctx, "magma"))
  if (!is.null(gs)) tools <- c(tools, writeup_tool(ctx, "twas_gsea"))
  ds <- writeup_tool(ctx, "msigdb")
  gstr(
    "Pathway enrichment against user-supplied gene-set collections ({ds}) was tested with ",
    "{oxford(tools)}. P-values were pooled across all gene sets and FDR-corrected."
  )
}

methods_tissue <- function(ctx, g, meta) {
  d <- pluck(g, "tissue", "specific")
  if (is.null(d) || is.null(d$res)) return(NULL)
  .ran(ctx, "tissue")
  tool <- writeup_tool(ctx, "magma")
  gtex <- writeup_tool(ctx, "gtex_v8")
  gstr(
    "Tissue-specific enrichment of gene-level association was tested with {tool} using ",
    "{gtex} tissue expression, with FDR correction across tissues followed by conditional ",
    "analysis of the significant tissues."
  )
}

# ---- closing paragraph ------------------------------------------------------

methods_closing <- function(ctx, manifest) {
  ver <- manifest$pipeline_tag %||% (manifest$pipeline_version %||% "unknown")
  commit <- manifest$pipeline_commit
  commit_clause <- if (!is.null(commit) && !is.na(commit)) gstr(" (commit {substr(as.character(commit), 1, 10)})") else ""
  gstr(
    "All analyses were orchestrated by the GenoDisc post-GWAS pipeline ",
    "(version {ver}{commit_clause}) and are available through the GenoDisc web server ",
    "[CITE: GenoDisc paper]. Results were explored with the GenoDisc Output Viewer ",
    "[CITE: GenoDisc paper]."
  )
}
