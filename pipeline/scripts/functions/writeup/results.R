#!/usr/bin/Rscript

# Results-paragraph builders — strictly descriptive. Each takes (ctx, g, meta) and
# returns a Markdown paragraph for a block that RAN (with explicit null-result
# sentences), or NULL for a block that did not run. "Number significant" is always
# the per-block FDR count (stored P.FDR < alpha), matching the pipeline tables and
# the Output Viewer. No interpretation, causal or novelty language.

library(glue)

ALPHA <- function(ctx) param(ctx, "fdr_alpha") %||% 0.05

# first existing column among `names`, else NULL
.col <- function(df, names) {
  for (n in names) if (!is.null(df) && n %in% names(df)) return(df[[n]])
  NULL
}

# top-N rows of df ordered by `order_col` (ascending P by default)
.top <- function(df, order_col, n = 3L, decreasing = FALSE) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || !order_col %in% names(df)) return(df[0, , drop = FALSE])
  v <- suppressWarnings(as.numeric(df[[order_col]]))
  ord <- order(v, decreasing = decreasing, na.last = NA)
  df[head(ord, n), , drop = FALSE]
}

# stable table reference string
.tref <- function(label, block_id, tab) {
  gstr("(Table W-{label}; see the {tab} of the Output Viewer, packaged block `{block_id}`)")
}

# ---- heritability -----------------------------------------------------------

results_ldsc_h2 <- function(ctx, g, meta) {
  v <- pluck(g, "gwas_qc", "ldsc_dat", "val"); if (is.null(v)) return(NULL)
  opt <- ctx$opt
  gstr(
    "LD-score regression estimated the observed-scale SNP-based heritability at ",
    "{fmt_num(v$obs_h2_est, opt)} (SE {fmt_num(v$obs_h2_se, opt)}) with an LD-score ",
    "regression intercept of {fmt_num(v$int_est, opt)} (SE {fmt_num(v$int_se, opt)}) ",
    "{.tref('ldsc-h2', 'gwas_qc', 'Heritability tab')}."
  )}

results_sbayesrc <- function(ctx, g, meta) {
  v <- pluck(g, "gwas_qc", "sbayesrc_dat", "val"); if (is.null(v)) return(NULL)
  opt <- ctx$opt
  gstr(
    "SBayesRC estimated the SNP-based heritability at {fmt_num(v$h2_mode, opt)} ",
    "(95% HPD {fmt_ci(v$h2_hpd_low, v$h2_hpd_high, opt)}), with {fmt_int(v$n_nonzero)} of ",
    "{fmt_int(v$n_snp_total)} variants estimated to have non-zero effect ",
    "{.tref('sbayesrc', 'gwas_qc', 'Heritability tab')}."
  )}

results_sumher <- function(ctx, g, meta) {
  v <- pluck(g, "gwas_qc", "sumher_dat", "val"); if (is.null(v)) return(NULL)
  opt <- ctx$opt
  warn <- if (isTRUE(v$missing_warning)) gstr(" (note: {fmt_pct(v$prop_missing, opt)} of reference tagging predictors were absent after QC)") else ""
  gstr(
    "SumHer estimated the SNP-based heritability at {fmt_num(v$h2, opt)} ",
    "(SE {fmt_num(v$h2_sd, opt)}) under the {v$model} model{warn} ",
    "{.tref('sumher', 'gwas_qc', 'Heritability tab')}."
  )}

results_gencor <- function(ctx, g, meta) {
  sec    <- pluck(g, "gwas_qc", "ldsc_gencor_dat", "table")
  within <- pluck(g, "gwas_qc", "ldsc_gencor_within_dat", "table")
  tab <- if (!is.null(sec)) sec else within
  if (is.null(tab)) return(NULL)
  opt <- ctx$opt; n <- nrow(tab)
  nsig <- count_below(tab, "rg_p_fdr", ALPHA(ctx))
  if (n == 0) return(gstr("No genetic correlations were estimated {.tref('gencor', 'gwas_qc', 'Genetic correlation tab')}."))
  top <- .top(tab, "rg_p_fdr", 1L)
  top_clause <- if (nrow(top) == 1) {
    lab <- .col(top, c("label", "name"))
    gstr(" The strongest correlation was with {lab[1]} (rg = {fmt_num(top$rg[1], opt)}, SE {fmt_num(top$rg_se[1], opt)}, P.FDR = {fmt_p(top$rg_p_fdr[1], opt)}).")
  } else ""
  null_clause <- if (nsig == 0) " No genetic correlation was significant after FDR correction." else ""
  gstr(
    "Genetic correlation was tested for {fmt_int(n)} trait {pluralise('pair', n)}; ",
    "{fmt_int(nsig)} {was_were(nsig)} significant after FDR correction (P.FDR < {fmt_num(ALPHA(ctx), opt)}).",
    "{top_clause}{null_clause} {.tref('gencor', 'gwas_qc', 'Genetic correlation tab')}"
  )}

# ---- SNP associations -------------------------------------------------------

results_clump <- function(ctx, g, meta) {
  d <- pluck(g, "snp_assoc", "clump"); if (is.null(d)) return(NULL)
  n <- nrow(d); opt <- ctx$opt
  if (n == 0) return(gstr("No independent genome-wide significant loci were identified {.tref('clump', 'snp_assoc', 'SNP associations tab')}."))
  top <- .top(d, "P", 1L)
  gene <- .col(top, c("NearestGene")); snp <- .col(top, c("SNP", "rsid"))
  lead <- if (!is.null(snp)) snp[1] else if (!is.null(gene)) gene[1] else "the lead variant"
  gstr(
    "Clumping identified {fmt_int(n)} independent genome-wide significant ",
    "{pluralise('locus', n, 'loci')}. The most significant was {lead}",
    "{if (!is.null(gene)) paste0(' (nearest gene ', gene[1], ')') else ''} ",
    "(P = {fmt_p(top$P[1], opt)}) {.tref('clump', 'snp_assoc', 'SNP associations tab')}."
  )}

results_cojo <- function(ctx, g, meta) {
  d <- pluck(g, "snp_assoc", "cojo"); if (is.null(d)) return(NULL)
  n <- nrow(d); opt <- ctx$opt
  st <- pluck(g, "snp_assoc", "cojo_status")
  fail_clause <- if (!is.null(st) && isTRUE(st$any_failed)) gstr(" ({fmt_int(st$n_failed)} chromosome(s) were skipped because the number of signals exceeded the LD reference size.)") else ""
  if (n == 0) return(gstr("No conditionally independent signals were identified by COJO{fail_clause} {.tref('cojo', 'snp_assoc', 'SNP associations tab')}."))
  pj <- if ("pJ" %in% names(d)) "pJ" else "P"
  top <- .top(d, pj, 1L)
  gene <- .col(top, c("NearestGene")); snp <- .col(top, c("SNP"))
  gstr(
    "COJO identified {fmt_int(n)} conditionally independent association ",
    "{pluralise('signal', n)}. The strongest joint signal was {if (!is.null(snp)) snp[1] else 'the top variant'}",
    "{if (!is.null(gene)) paste0(' (nearest gene ', gene[1], ')') else ''} ",
    "(joint P = {fmt_p(top[[pj]][1], opt)}){fail_clause} {.tref('cojo', 'snp_assoc', 'SNP associations tab')}."
  )}

results_finemap <- function(ctx, g, meta) {
  d <- pluck(g, "snp_assoc", "susie", "L1"); if (is.null(d)) return(NULL)
  opt <- ctx$opt
  # a real credible set has a non-NA cs (the no-signal case is a single all-NA row)
  cs <- .col(d, c("cs"))
  has_cs <- if (!is.null(cs)) !is.na(cs) else rep(FALSE, nrow(d))
  n <- sum(has_cs)
  if (n == 0) return(gstr("No 95% credible sets were resolved by fine-mapping {.tref('finemap', 'snp_assoc', 'SNP associations tab')}."))
  dd <- d[has_cs, , drop = FALSE]
  top <- .top(dd, "TopPIP", 1L, decreasing = TRUE)
  gene <- .col(top, c("Gene"))
  gstr(
    "Fine-mapping resolved {fmt_int(n)} credible {pluralise('set', n)} at the 95% level. ",
    "The credible set with the highest posterior inclusion probability ",
    "(TopPIP = {fmt_num(top$TopPIP[1], opt)}) mapped to {if (!is.null(gene)) gene[1] else 'an unannotated region'} ",
    "{.tref('finemap', 'snp_assoc', 'SNP associations tab')}."
  )}

# ---- molecular associations -------------------------------------------------

# shared descriptive sentence for a FDR + (optional) support column molecular block
.mol_result <- function(ctx, df, fdr_col, p_col, label_col, block_id, tab, analysis,
                        support_col = NULL, support_desc = NULL, support_rule = NULL) {
  opt <- ctx$opt; n <- nrow(df)
  if (n == 0) return(gstr("No {analysis} associations were tested {.tref(block_id, block_id, tab)}."))
  nsig <- count_below(df, fdr_col, ALPHA(ctx))
  support_clause <- ""
  if (!is.null(support_col) && support_col %in% names(df) && nsig > 0) {
    sig <- df[suppressWarnings(as.numeric(df[[fdr_col]])) < ALPHA(ctx) & !is.na(df[[fdr_col]]), , drop = FALSE]
    nsupp <- if (is.null(support_rule)) count_below(sig, support_col, ALPHA(ctx)) else sum(support_rule(suppressWarnings(as.numeric(sig[[support_col]]))), na.rm = TRUE)
    support_clause <- gstr(", of which {fmt_int(nsupp)} also {support_desc}")
  }
  if (nsig == 0) {
    return(gstr(
      "{analysis} tested {fmt_int(n)} associations; none were significant after FDR ",
      "correction (P.FDR < {fmt_num(ALPHA(ctx), opt)}) {.tref(block_id, block_id, tab)}."
    ))
  }
  top <- .top(df, p_col, 1L)
  lab <- .col(top, label_col)
  panel <- .col(top, c("PANEL", "Panel"))
  panel_clause <- if (!is.null(panel)) gstr(" in {panel[1]}") else ""
  gstr(
    "{analysis} tested {fmt_int(n)} associations; {fmt_int(nsig)} {was_were(nsig)} significant ",
    "after FDR correction (P.FDR < {fmt_num(ALPHA(ctx), opt)}){support_clause}. The strongest ",
    "association was {if (!is.null(lab)) lab[1] else 'the top feature'}{panel_clause} ",
    "(P = {fmt_p(top[[p_col]][1], opt)}, P.FDR = {fmt_p(top[[fdr_col]][1], opt)}) ",
    "{.tref(block_id, block_id, tab)}."
  )}

results_twas <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "exp", "fusion", "res"); if (is.null(d)) return(NULL)
  .mol_result(ctx, d, "TWAS.P.FDR", "TWAS.P", c("Gene Symbol", "Ensembl ID"),
              "mol_assoc/exp/fusion", "TWAS tab", "Transcriptome-wide association analysis",
              support_col = "COLOC_logical", support_desc = "showed colocalisation support",
              support_rule = function(v) v == 1 | v == TRUE)
}

results_pwas <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "protein", "fusion", "results"); if (is.null(d)) return(NULL)
  .mol_result(ctx, d, "pwas_all.P.FDR", "pwas_all.P", c("Gene Symbol", "Ensembl ID"),
              "mol_assoc/protein/fusion", "PWAS tab", "Proteome-wide association analysis",
              support_col = "COLOC_logical", support_desc = "showed colocalisation support",
              support_rule = function(v) v == 1 | v == TRUE)
}

results_smr_exp <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "exp", "smr", "results"); if (is.null(d)) return(NULL)
  .mol_result(ctx, d, "p_SMR.FDR", "p_SMR", c("Gene Symbol", "Ensembl ID"),
              "mol_assoc/exp/smr", "SMR tab", "SMR of gene expression",
              support_col = "p_HEIDI", support_desc = "passed the HEIDI test (p_HEIDI > 0.05)",
              support_rule = function(v) v > 0.05)
}

results_smr_protein <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "protein", "smr", "results"); if (is.null(d)) return(NULL)
  .mol_result(ctx, d, "p_SMR.FDR", "p_SMR", c("Gene Symbol", "Ensembl ID"),
              "mol_assoc/protein/smr", "SMR tab", "SMR of protein levels",
              support_col = "p_HEIDI", support_desc = "passed the HEIDI test (p_HEIDI > 0.05)",
              support_rule = function(v) v > 0.05)
}

results_magma_gene <- function(ctx, g, meta) {
  d <- pluck(g, "mol_assoc", "magma"); if (is.null(d)) return(NULL)
  .mol_result(ctx, d, "P.FDR", "P", c("ID"),
              "mol_assoc/magma", "MAGMA tab", "MAGMA gene-level analysis")
}

# ---- enrichment / pathway / tissue ------------------------------------------

.enrich_line <- function(ctx, df, name_col, analysis, block_id) {
  if (is.null(df) || nrow(df) == 0) return(gstr("{analysis} returned no results."))
  n <- nrow(df); nsig <- count_below(df, "P.FDR", ALPHA(ctx))
  if (nsig == 0) return(gstr("{analysis} tested {fmt_int(n)} sets; none were FDR-significant."))
  top <- .top(df, "P", 1L); lab <- .col(top, name_col)
  gstr("{analysis} tested {fmt_int(n)} sets, {fmt_int(nsig)} FDR-significant (top: {if (!is.null(lab)) lab[1] else 'NA'}, P.FDR = {fmt_p(top$P.FDR[1], ctx$opt)}).")
}

results_enrichment <- function(ctx, g, meta) {
  lines <- character(0)
  md <- pluck(g, "tx", "drug", "magma"); if (!is.null(md)) lines <- c(lines, .enrich_line(ctx, md, c("Name"), "MAGMA drug enrichment", "tx/drug/magma"))
  gd <- pluck(g, "tx", "drug", "twas_gsea"); if (!is.null(gd)) lines <- c(lines, .enrich_line(ctx, gd, c("Name"), "TWAS-GSEA drug enrichment", "tx/drug/twas_gsea"))
  cd <- pluck(g, "tx", "drug", "gcsc"); if (!is.null(cd)) lines <- c(lines, .enrich_line(ctx, cd, c("Name"), "GCSC drug enrichment", "tx/drug/gcsc"))
  cm <- pluck(g, "tx", "cmap", "drug"); if (!is.null(cm)) lines <- c(lines, .enrich_line(ctx, cm, c("Name"), "CMAP signature enrichment", "tx/cmap/drug"))
  if (length(lines) == 0) return(NULL)
  paste0(paste(lines, collapse = " "),
         gstr(" (see the Drug enrichment tab of the Output Viewer.)"))}

results_pathway <- function(ctx, g, meta) {
  mg <- pluck(g, "tx", "pathway", "magma"); gs <- pluck(g, "tx", "pathway", "twas_gsea")
  if (is.null(mg) && is.null(gs)) return(NULL)
  lines <- character(0)
  if (!is.null(mg)) lines <- c(lines, .enrich_line(ctx, mg, c("Name"), "MAGMA pathway enrichment", "tx/pathway/magma"))
  if (!is.null(gs)) lines <- c(lines, .enrich_line(ctx, gs, c("Name"), "TWAS-GSEA pathway enrichment", "tx/pathway/twas_gsea"))
  paste0(paste(lines, collapse = " "), " (see the Pathway tab of the Output Viewer.)")}

results_tissue <- function(ctx, g, meta) {
  d <- pluck(g, "tissue", "specific"); if (is.null(d) || is.null(d$res)) return(NULL)
  res <- d$res; n <- nrow(res); opt <- ctx$opt
  nsig <- count_below(res, "P.FDR", ALPHA(ctx))
  keep <- d$keep
  keep_clause <- if (!is.null(keep) && length(keep) > 0) gstr(" After conditional analysis, {fmt_int(length(keep))} tissue(s) remained: {oxford(as.character(keep))}.") else ""
  if (n == 0) return(gstr("No tissues were tested for enrichment {.tref('tissue', 'tissue', 'Tissue tab')}."))
  if (nsig == 0) return(gstr("MAGMA tissue enrichment tested {fmt_int(n)} tissues; none were FDR-significant {.tref('tissue', 'tissue', 'Tissue tab')}."))
  top <- .top(res, "P", 1L); lab <- .col(top, c("Tissue"))
  gstr(
    "MAGMA tissue enrichment tested {fmt_int(n)} tissues; {fmt_int(nsig)} {was_were(nsig)} ",
    "FDR-significant. The most enriched was {if (!is.null(lab)) lab[1] else 'NA'} ",
    "(P = {fmt_p(top$P[1], opt)}, P.FDR = {fmt_p(top$P.FDR[1], opt)}).{keep_clause} ",
    "{.tref('tissue', 'tissue', 'Tissue tab')}"
  )}
