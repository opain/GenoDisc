#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(scales)
})

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if(length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})
source(file.path(.script_dir, "_paths.R"))
tdir <- tables_dir
fdir <- figures_dir

summ <- fread(file.path(tdir, "upstream_summary_per_trait.tsv"))
trait_order <- c("ADHD","ASD","BIP","MDD","SCZ","ALZ","PRK","ALS","MIG")
summ[, Trait := factor(Trait, levels = trait_order)]

############################
# Fig U1: per-trait counts by module (log scale)
############################
melt_summ <- melt(summ, id.vars = "Trait")
# Drop the credible-set count columns (we plot genes instead) and keep useful counts
keep_metrics <- c(
  "Clump_GWsig_loci"            = "GW-sig clumped loci",
  "COJO_independent_signals"    = "COJO independent signals",
  "SuSiE_L10_genes"             = "SuSiE L=10 mapped genes",
  "TWAS_FUSION_hc_genes"        = "TWAS-FUSION HC genes (FDR & COLOC)",
  "SMR_expression_hc_genes"     = "SMR expression HC genes (FDR & HEIDI)",
  "PWAS_FUSION_hc_genes"        = "PWAS-FUSION HC genes (FDR & COLOC)",
  "SMR_protein_hc_genes"        = "SMR protein HC genes (FDR & HEIDI)",
  "MAGMA_FDR_genes"             = "MAGMA gene-level FDR<0.05",
  "Tissue_FDR_enriched"         = "Enriched tissues (FDR<0.05)",
  "Tissue_retained_conditional" = "Tissues retained after conditional analysis"
)
melt_keep <- melt_summ[as.character(variable) %in% names(keep_metrics)]
melt_keep[, label := factor(keep_metrics[as.character(variable)],
                            levels = unname(keep_metrics))]

p_u1 <- ggplot(melt_keep, aes(x = Trait, y = pmax(value, 0.5), fill = Trait)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = value), vjust = -0.3, size = 2.6) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  scale_y_log10(labels = comma) +
  labs(title = "Per-trait upstream signal counts",
       subtitle = "y-axis log scale (values of 0 plotted at 0.5 for visibility).\nHC = high-confidence: TWAS/PWAS FDR<0.05 AND COLOC_logical=TRUE; SMR FDR<0.05 AND HEIDI>0.01.",
       x = NULL, y = "Count (log scale)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fdir, "figU1_upstream_counts_per_trait.png"),
       p_u1, width = 13, height = 9, dpi = 200)

############################
# Fig U2: cross-trait overlap heatmap, per method
############################
ov <- fread(file.path(tdir, "cross_trait_overlap_genes.tsv"))
ov[, TraitA := factor(TraitA, levels = trait_order)]
ov[, TraitB := factor(TraitB, levels = trait_order)]
# Keep upper triangle for readability
ov <- ov[as.integer(TraitA) <= as.integer(TraitB)]
plot_method <- function(method_name){
  d <- ov[Method == method_name]
  if(nrow(d) == 0) return(NULL)
  ggplot(d, aes(x = TraitA, y = TraitB, fill = N_shared)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = N_shared), size = 2.8) +
    scale_fill_gradient(low = "white", high = "#b2182b", na.value = "grey90") +
    labs(title = method_name, x = NULL, y = NULL, fill = "Shared\ngenes") +
    theme_bw(base_size = 10) +
    theme(legend.position = "right", axis.text.x = element_text(angle = 30, hjust = 1))
}
plist <- lapply(c("SuSiE_L10","TWAS_FUSION","SMR_expression","PWAS_FUSION","SMR_protein","MAGMA"),
                plot_method)
plist <- plist[!sapply(plist, is.null)]
if(requireNamespace("gridExtra", quietly = TRUE)){
  png(file.path(fdir, "figU2_cross_trait_overlap.png"), width = 2400, height = 1600, res = 200)
  do.call(gridExtra::grid.arrange, c(plist, ncol = 3))
  dev.off()
}

############################
# Fig U3: recurring high-confidence genes across traits (top 30 by N_traits)
############################
recur <- fread(file.path(tdir, "cross_trait_recurring_genes.tsv"))
plot_recurring <- function(method_name, max_genes = 30){
  d <- recur[Method == method_name][order(-N_traits, Gene)][seq_len(min(.N, max_genes))]
  if(nrow(d) == 0) return(NULL)
  # Pull per-trait presence
  long_g <- fread(file.path(tdir, "upstream_high_confidence_genes_long.tsv"))
  m <- long_g[Method == method_name][Gene %in% d$Gene]
  m[, Present := 1L]
  w <- dcast(m, Gene ~ Trait, value.var = "Present", fill = 0L)
  long <- melt(w, id.vars = "Gene", variable.name = "Trait", value.name = "Present")
  long[, Trait := factor(Trait, levels = trait_order)]
  long[, Gene := factor(Gene, levels = rev(d$Gene))]
  ggplot(long, aes(x = Trait, y = Gene, fill = factor(Present))) +
    geom_tile(colour = "white") +
    scale_fill_manual(values = c(`0` = "grey95", `1` = "#b2182b"), guide = "none") +
    labs(title = method_name, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7))
}
plist <- lapply(c("SuSiE_L10","TWAS_FUSION","SMR_expression","PWAS_FUSION","SMR_protein","MAGMA"),
                plot_recurring)
plist <- plist[!sapply(plist, is.null)]
if(requireNamespace("gridExtra", quietly = TRUE)){
  png(file.path(fdir, "figU3_recurring_genes_across_traits.png"), width = 2400, height = 2400, res = 200)
  do.call(gridExtra::grid.arrange, c(plist, ncol = 3))
  dev.off()
}

############################
# Fig U4: tissue enrichment across traits
############################
tw <- fread(file.path(tdir, "tissue_enrichment_all_traits.tsv"))
tw[, Trait := factor(Trait, levels = trait_order)]
# Restrict to tissues nominally significant somewhere
sig_tissues <- unique(tw[P < 0.05, Tissue])
tws <- tw[Tissue %in% sig_tissues]
# Order tissues by number of traits where they're FDR-sig
tissue_rank <- tws[, .(N_fdr = sum(FDR_sig, na.rm=TRUE)), by = Tissue][order(-N_fdr)]
tws[, Tissue := factor(Tissue, levels = rev(tissue_rank$Tissue))]
tws[, neglog10P := -log10(P)]
tws[, neglog10P_cap := pmin(neglog10P, 10)]
p_u4 <- ggplot(tws, aes(x = Trait, y = Tissue, fill = neglog10P_cap)) +
  geom_tile(colour = "grey95") +
  geom_text(data = tws[P < 0.05 & P.FDR >= 0.05], aes(label = "*"), size = 3) +
  geom_text(data = tws[P.FDR < 0.05 & !Retained], aes(label = "**"), size = 3) +
  geom_text(data = tws[Retained == TRUE], aes(label = "++"), size = 3, colour="black") +
  scale_fill_gradient(low = "white", high = "#b2182b", na.value = "grey95",
                      name = "-log10(P)\n(capped at 10)") +
  labs(title = "MAGMA tissue-specific enrichment per trait",
       subtitle = "* nominal P<0.05; ** FDR<0.05; ++ retained after conditional analysis (= robust tissue signal).",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(fdir, "figU4_tissue_enrichment.png"), p_u4,
       width = 11, height = max(7, 0.16 * length(sig_tissues) + 3), dpi = 200, limitsize = FALSE)

cat("Saved upstream figures:\n")
cat(list.files(fdir, pattern = "figU.*\\.png"), sep = "\n")
