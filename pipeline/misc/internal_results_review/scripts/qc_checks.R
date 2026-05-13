#!/usr/bin/env Rscript
suppressMessages({ library(data.table) })

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if(length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})
source(file.path(.script_dir, "_paths.R"))
tdir   <- tables_dir
out_md <- file.path(report_dir, "qc_checks.md")

gd  <- fread(file.path(tdir, "twas_gsea_drug_all_traits.tsv"))
ga  <- fread(file.path(tdir, "twas_gsea_atc_all_traits.tsv"))
md  <- fread(file.path(tdir, "magma_drug_all_traits.tsv"))
ma  <- fread(file.path(tdir, "magma_atc_all_traits.tsv"))
gcd <- fread(file.path(tdir, "gcsc_drug_all_traits.tsv"))
cm  <- fread(file.path(tdir, "cmap_moa_all_traits.tsv"))
cd  <- fread(file.path(tdir, "cmap_drug_all_traits.tsv"))

sink(out_md)
cat("# Quality-control summaries (multi-panel)\n\n")

cat("## 1. Direction invariant check (sign(Reversal_Z) agrees with Direction)\n\n")
for(dat_name in list(list(name="TWAS-GSEA drug", d=gd),
                     list(name="TWAS-GSEA ATC",  d=ga),
                     list(name="CMAP MOA",       d=cm))){
  d <- as.data.table(dat_name$d)[!is.na(Direction)]
  bad <- sum((d$Direction == "Opposes disease" & d$Reversal_Z < 0) |
             (d$Direction == "Matches disease" & d$Reversal_Z > 0))
  zero_z <- sum(d$Reversal_Z == 0 & !is.na(d$Direction))
  cat(sprintf("- %s: %d bad / %d rows (zero-Z edge cases: %d, all from P=1 Wilcoxon)\n",
              dat_name$name, bad, nrow(d), zero_z))
}

cat("\n## 2. Gene-set size vs significance, TWAS-GSEA drug\n\n")
cat("Spearman correlation between -log10(P) and N_Genes per trait per panel:\n\n")
qc <- gd[, .(rho = suppressWarnings(cor(-log10(P), N_Genes, method="spearman", use="pairwise.complete.obs")),
             N = .N), by = .(Trait, Panel)]
print(qc)

cat("\n## 3. Oncology drug overrepresentation among top-30 TWAS-GSEA drug hits\n\n")
gd[, Oncology := grepl("antineopl|kinase inhibit|cytotoxic|antimetab|alkylating", ATC_Description, ignore.case=TRUE)]
qc3 <- gd[, .(
  pct_overall_oncology = round(mean(Oncology, na.rm=TRUE)*100, 1),
  pct_top30_oncology   = round(mean(Oncology[order(P)][1:30], na.rm=TRUE)*100, 1),
  pct_top100_oncology  = round(mean(Oncology[order(P)][1:100], na.rm=TRUE)*100, 1)
), by = .(Trait)]
print(qc3)

cat("\n## 4. Drugs with very large target sets (N_Genes > 500)\n\n")
cat("By construction, large-target drugs are often pleiotropic. Check whether any reach FDR<0.05:\n\n")
print(gd[N_Genes > 500 & P.FDR < 0.05, .(Trait, Panel, Name, N_Genes, Estimate, P, P.FDR, Direction)])

cat("\n## 5. Direction consistency for cross-trait recurring ATC classes\n\n")
ga_best <- ga[order(P), .SD[1], by = .(Trait, ATC_Code)]
ga_best[, Nom := P < 0.05]
cross <- ga_best[, .(n_nom = sum(Nom, na.rm=TRUE),
                     n_opposes  = sum(Nom & Direction == "Opposes disease", na.rm=TRUE),
                     n_matches  = sum(Nom & Direction == "Matches disease", na.rm=TRUE),
                     min_FDR = min(P.FDR, na.rm=TRUE)),
                 by = .(ATC_Code, ATC_Description)]
cross <- cross[n_nom >= 3][order(-n_nom, min_FDR)]
cross[, dir_consistency := round(abs(n_opposes - n_matches) / pmax(n_nom, 1), 2)]
print(head(cross, 25))

cat("\n## 6. Panel concordance: Spearman correlation between panels on TWAS-GSEA ATC Reversal_Z\n\n")
panels <- unique(ga$Panel)
ppairs <- t(combn(panels, 2))
res <- list()
for(tr in unique(ga$Trait)){
  sub <- dcast(ga[Trait == tr, .(Reversal_Z, Panel, ATC_Code)],
               ATC_Code ~ Panel, value.var = "Reversal_Z")
  for(i in seq_len(nrow(ppairs))){
    p1 <- ppairs[i,1]; p2 <- ppairs[i,2]
    if(!all(c(p1, p2) %in% names(sub))) next
    x1 <- sub[[p1]]; x2 <- sub[[p2]]
    keep <- !is.na(x1) & !is.na(x2)
    if(sum(keep) < 5) next
    rho <- suppressWarnings(cor(x1[keep], x2[keep], method="spearman"))
    res[[length(res)+1]] <- data.table(Trait=tr, Panel1=p1, Panel2=p2, rho=round(rho,3))
  }
}
print(rbindlist(res))

cat("\n## 7. PRK CMAP outlier check\n\n")
cat("PRK has 4781 FDR-sig CMAP signatures vs 0-16 for other traits. Check distribution.\n\n")
cat("Top-10 P-values of PRK CMAP drug per panel:\n")
print(cd[Trait == "PRK"][order(P)][, head(.SD, 3), by = Panel][, .(Panel, cmap_name, cell_iname, P)])
cat("\nPRK CMAP P distribution by panel (quantiles):\n")
print(cd[Trait == "PRK", .(N=.N, P_1pct = quantile(P,0.01), P_5pct = quantile(P,0.05), median_P = median(P)), by = Panel])

cat("\n## 8. Cross-method drug-level overlap (FDR<0.05)\n\n")
md_cmp <- fread(file.path(tdir, "method_comparison_drug.tsv"))
ovd <- md_cmp[, .(
  GSEA_FDR_sig  = sum(gsea_FDR_sig, na.rm=TRUE),
  MAGMA_FDR_sig = sum(magma_FDR_sig, na.rm=TRUE),
  GCSC_FDR_sig  = sum(gcsc_FDR_sig, na.rm=TRUE),
  Any2          = sum((gsea_FDR_sig + magma_FDR_sig + gcsc_FDR_sig) >= 2, na.rm=TRUE),
  All3          = sum(gsea_FDR_sig & magma_FDR_sig & gcsc_FDR_sig, na.rm=TRUE)
), by = Trait]
print(ovd)

cat("\n## 9. Cross-method ATC-level overlap (FDR<0.05)\n\n")
ma_cmp <- fread(file.path(tdir, "method_comparison_atc.tsv"))
ova <- ma_cmp[, .(
  GSEA_FDR_sig  = sum(gsea_FDR_sig, na.rm=TRUE),
  MAGMA_FDR_sig = sum(magma_FDR_sig, na.rm=TRUE),
  GCSC_FDR_sig  = sum(gcsc_FDR_sig, na.rm=TRUE),
  Any2          = sum((gsea_FDR_sig + magma_FDR_sig + gcsc_FDR_sig) >= 2, na.rm=TRUE),
  All3          = sum(gsea_FDR_sig & magma_FDR_sig & gcsc_FDR_sig, na.rm=TRUE)
), by = Trait]
print(ova)

sink()
cat("Wrote: ", out_md, "\n")
