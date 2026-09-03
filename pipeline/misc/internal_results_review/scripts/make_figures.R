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

gd  <- fread(file.path(tdir, "twas_gsea_drug_all_traits.tsv"))
ga  <- fread(file.path(tdir, "twas_gsea_atc_all_traits.tsv"))
md  <- fread(file.path(tdir, "magma_drug_all_traits.tsv"))
ma  <- fread(file.path(tdir, "magma_atc_all_traits.tsv"))
gcd <- fread(file.path(tdir, "gcsc_drug_all_traits.tsv"))
gca <- fread(file.path(tdir, "gcsc_atc_all_traits.tsv"))
cm  <- fread(file.path(tdir, "cmap_moa_all_traits.tsv"))
cd  <- fread(file.path(tdir, "cmap_drug_all_traits.tsv"))
sm  <- fread(file.path(tdir, "hits_per_trait_per_method.tsv"))

trait_order <- c("ADHD","ASD","BIP","MDD","SCZ","ALZ","PRK","ALS","MIG")
for(x in list(gd, ga, md, ma, gcd, gca, cm, cd, sm)) x[, Trait := factor(Trait, levels = trait_order)]

############################
# Fig 1: hits per trait, per method
############################
sm[, Method := factor(Method, levels = c(
  "TWAS-GSEA drug (dir)","TWAS-GSEA drug (nondir)","TWAS-GSEA ATC (dir)",
  "MAGMA drug","MAGMA ATC","GCSC drug","GCSC ATC",
  "CMAP drug (per-sig)","CMAP per-MOA"))]
sm_long <- melt(sm, id.vars = c("Trait","Method"),
                variable.name = "Sig", value.name = "N")
p1 <- ggplot(sm_long, aes(x = Trait, y = N, fill = Sig)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = N), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.3) +
  facet_wrap(~ Method, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Nominal = "grey60", FDR_sig = "#b2182b")) +
  labs(title = "Hits per trait by method",
       subtitle = "Grey = nominal P<0.05; red = within-trait FDR<0.05. TWAS-GSEA / CMAP drug & ATC counts are summed across panels.",
       x = NULL, y = "Number of drugs / ATC classes / signatures / MOAs") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fdir, "fig1_hits_per_trait_by_method.png"),
       p1, width = 13, height = 11, dpi = 200)

############################
# Fig 2: TWAS-GSEA ATC heatmap, best panel per (trait, ATC)
############################
ga_best <- ga[order(P), .SD[1], by = .(Trait, ATC_Code, ATC_Description)]
top_atcs <- unique(ga_best[P < 0.05, ATC_Code])
gah <- ga_best[ATC_Code %in% top_atcs]
gah[, ATC_label := str_trunc(paste0(ATC_Code, ": ", ATC_Description), 65)]
ord <- gah[P < 0.05, .N, by = ATC_label][order(-N)]
gah[, ATC_label := factor(ATC_label, levels = rev(ord$ATC_label))]
gah[, Zcap := pmin(pmax(Reversal_Z, -4), 4)]

p2 <- ggplot(gah, aes(x = Trait, y = ATC_label, fill = Zcap)) +
  geom_tile(colour = "grey90") +
  geom_text(data = gah[P < 0.05 & P.FDR >= 0.05], aes(label = "*"), size = 4) +
  geom_text(data = gah[P.FDR < 0.05], aes(label = "**"), size = 4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, na.value = "grey90",
                       name = "Reversal Z\n(+ opposes,\n- matches)",
                       limits = c(-4, 4), oob = squish) +
  labs(title = "TWAS-GSEA ATC findings across traits (best panel per cell)",
       subtitle = "* nominal P<0.05; ** within-trait FDR<0.05. Cells take the panel with smallest P for each (trait, ATC).",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(fdir, "fig2_twas_gsea_atc_heatmap.png"), p2,
       width = 11, height = max(6, 0.18 * length(top_atcs) + 3), dpi = 200, limitsize = FALSE)

############################
# Fig 3: cross-trait drug-level TWAS-GSEA heatmap (drugs nominally sig in >= 2 traits, best panel per cell)
############################
gd_best <- gd[order(P), .SD[1], by = .(Trait, Name_clean)]
drug_n_traits <- gd_best[, .(N_traits_nom = sum(P < 0.05, na.rm=TRUE)), by = Name_clean]
recur <- drug_n_traits[N_traits_nom >= 2, Name_clean]
gdr <- gd_best[Name_clean %in% recur]
gdr[, Drug_label := str_to_title(Name_clean)]
ord_d <- gdr[P < 0.05, .N, by = Drug_label][order(-N)]
gdr[, Drug_label := factor(Drug_label, levels = rev(ord_d$Drug_label))]
gdr[, Zcap := pmin(pmax(Reversal_Z, -4), 4)]
p3 <- ggplot(gdr, aes(x = Trait, y = Drug_label, fill = Zcap)) +
  geom_tile(colour = "grey90") +
  geom_text(data = gdr[P < 0.05 & P.FDR >= 0.05], aes(label = "*"), size = 3) +
  geom_text(data = gdr[P.FDR < 0.05], aes(label = "**"), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, na.value = "grey90",
                       name = "Reversal Z",
                       limits = c(-4, 4), oob = squish) +
  labs(title = "TWAS-GSEA drugs nominally significant in >=2 traits (best panel per cell)",
       subtitle = "* nominal P<0.05; ** within-trait FDR<0.05.",
       x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.5))
ggsave(file.path(fdir, "fig3_drug_cross_trait_heatmap.png"), p3,
       width = 9, height = max(6, 0.13 * length(recur) + 3), dpi = 200, limitsize = FALSE)

############################
# Fig 4: panel concordance for TWAS-GSEA ATC
############################
panels <- unique(ga$Panel)
gar <- ga[, .(Reversal_Z, P), by = .(Trait, ATC_Code, ATC_Description, Panel)]
# Wide form: one column per panel
gaw <- dcast(gar, Trait + ATC_Code + ATC_Description ~ Panel, value.var = "Reversal_Z")
# Compute pairwise Spearman between panels' Reversal_Z within each trait
ppairs <- t(combn(panels, 2))
res <- list()
for(tr in unique(ga$Trait)){
  sub <- gaw[Trait == tr]
  for(i in seq_len(nrow(ppairs))){
    p1n <- ppairs[i,1]; p2n <- ppairs[i,2]
    x1 <- sub[[p1n]]; x2 <- sub[[p2n]]
    keep <- !is.na(x1) & !is.na(x2)
    if(sum(keep) < 5) next
    rho <- suppressWarnings(cor(x1[keep], x2[keep], method = "spearman"))
    res[[length(res)+1]] <- data.table(Trait = tr, Panel1 = p1n, Panel2 = p2n, rho = rho, n = sum(keep))
  }
}
panel_cor <- rbindlist(res)
panel_cor[, Trait := factor(Trait, levels = trait_order)]
panel_cor[, Pair := paste(pmin(Panel1, Panel2), pmax(Panel1, Panel2), sep = " vs\n")]
p4 <- ggplot(panel_cor, aes(x = Trait, y = Pair, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.6) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman\nrho") +
  labs(title = "Cross-panel concordance of TWAS-GSEA ATC Reversal_Z",
       subtitle = "Higher rho = panels agree on which classes are opposes/matches direction.",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(fdir, "fig4_panel_concordance.png"), p4,
       width = 11, height = 5, dpi = 200)

############################
# Fig 5: Method FDR-sig overlap (drug and ATC, per trait)
############################
drug_ov <- fread(file.path(tdir, "method_comparison_drug.tsv"))
atc_ov  <- fread(file.path(tdir, "method_comparison_atc.tsv"))

ov_drug <- drug_ov[, .(
  GSEA_only  = sum(gsea_FDR_sig & !magma_FDR_sig & !gcsc_FDR_sig, na.rm = TRUE),
  MAGMA_only = sum(magma_FDR_sig & !gsea_FDR_sig & !gcsc_FDR_sig, na.rm = TRUE),
  GCSC_only  = sum(gcsc_FDR_sig & !gsea_FDR_sig & !magma_FDR_sig, na.rm = TRUE),
  Any2_or_more = sum((gsea_FDR_sig + magma_FDR_sig + gcsc_FDR_sig) >= 2, na.rm = TRUE)
), by = Trait]
ov_drug_long <- melt(ov_drug, id.vars = "Trait", variable.name = "Overlap", value.name = "N")
ov_drug_long[, Trait := factor(Trait, levels = trait_order)]

ov_atc <- atc_ov[, .(
  GSEA_only  = sum(gsea_FDR_sig & !magma_FDR_sig & !gcsc_FDR_sig, na.rm = TRUE),
  MAGMA_only = sum(magma_FDR_sig & !gsea_FDR_sig & !gcsc_FDR_sig, na.rm = TRUE),
  GCSC_only  = sum(gcsc_FDR_sig & !gsea_FDR_sig & !magma_FDR_sig, na.rm = TRUE),
  Any2_or_more = sum((gsea_FDR_sig + magma_FDR_sig + gcsc_FDR_sig) >= 2, na.rm = TRUE)
), by = Trait]
ov_atc_long <- melt(ov_atc, id.vars = "Trait", variable.name = "Overlap", value.name = "N")
ov_atc_long[, Trait := factor(Trait, levels = trait_order)]

p5a <- ggplot(ov_drug_long, aes(x = Trait, y = N, fill = Overlap)) +
  geom_col(position = position_stack()) +
  scale_fill_manual(values = c(GSEA_only = "#fdae61", MAGMA_only = "#4575b4",
                               GCSC_only = "#74add1", Any2_or_more = "#1a9850")) +
  labs(title = "Drug-level FDR<0.05 hits per trait by method",
       y = "Number of drugs", x = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

p5b <- ggplot(ov_atc_long, aes(x = Trait, y = N, fill = Overlap)) +
  geom_col(position = position_stack()) +
  scale_fill_manual(values = c(GSEA_only = "#fdae61", MAGMA_only = "#4575b4",
                               GCSC_only = "#74add1", Any2_or_more = "#1a9850")) +
  labs(title = "ATC-level FDR<0.05 hits per trait by method",
       y = "Number of ATC classes", x = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

library(grid)
png(file.path(fdir, "fig5_method_overlap.png"), width = 2200, height = 1400, res = 200)
gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
if(gridExtra_ok){
  gridExtra::grid.arrange(p5a, p5b, ncol = 1)
} else {
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(p5a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p5b, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
}
dev.off()

############################
# Fig 6: CMAP MOA cross-trait heatmap (best panel/cell-line per (trait, MOA))
############################
cm_best <- cm[order(P), .SD[1], by = .(Trait, MOA)]
top_moas <- unique(cm_best[P < 0.05, MOA])
cmh <- cm_best[MOA %in% top_moas]
cmh[, MOA_label := str_trunc(MOA, 60)]
ord_m <- cmh[P < 0.05, .N, by = MOA_label][order(-N)]
cmh[, MOA_label := factor(MOA_label, levels = rev(ord_m$MOA_label))]
cmh[, Zcap := pmin(pmax(Reversal_Z, -4), 4)]
# Limit to MOAs sig in >=2 traits for readability
recur_moa <- cmh[P < 0.05, .N, by = MOA][N >= 2, MOA]
cmh_r <- cmh[MOA %in% recur_moa]

p6 <- ggplot(cmh_r, aes(x = Trait, y = MOA_label, fill = Zcap)) +
  geom_tile(colour = "grey90") +
  geom_text(data = cmh_r[P < 0.05 & P.FDR >= 0.05], aes(label = "*"), size = 3) +
  geom_text(data = cmh_r[P.FDR < 0.05], aes(label = "**"), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, na.value = "grey90",
                       name = "Reversal Z",
                       limits = c(-4, 4), oob = squish) +
  labs(title = "CMAP MOAs nominally significant in >=2 traits (best panel/cell per cell)",
       subtitle = "* nominal P<0.05; ** within-trait FDR<0.05.",
       x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.5))
ggsave(file.path(fdir, "fig6_cmap_moa_cross_trait.png"), p6,
       width = 9, height = max(6, 0.18 * length(recur_moa) + 3), dpi = 200, limitsize = FALSE)

############################
# Fig 7: gene-set size vs significance (TWAS-GSEA drug)
############################
p7 <- ggplot(gd[!is.na(P)], aes(x = N_Genes, y = -log10(P))) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  scale_x_log10() +
  facet_grid(Panel ~ Trait, scales = "free_y") +
  labs(title = "TWAS-GSEA drug: -log10(P) vs gene-set size",
       x = "Drug-target gene-set size (log)",
       y = "-log10(P)") +
  theme_bw(base_size = 8) + theme(strip.text = element_text(size = 6))
ggsave(file.path(fdir, "fig7_gene_set_size_vs_P.png"), p7,
       width = 14, height = 8, dpi = 200)

cat("Figures saved to:\n")
cat(list.files(fdir, pattern = "fig.*\\.png"), sep = "\n")
