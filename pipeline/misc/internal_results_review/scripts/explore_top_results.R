#!/usr/bin/env Rscript
suppressMessages({ library(data.table) })

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if(length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})
source(file.path(.script_dir, "_paths.R"))
tdir <- tables_dir
gd <- fread(file.path(tdir, "twas_gsea_drug_all_traits.tsv"))
ga <- fread(file.path(tdir, "twas_gsea_atc_all_traits.tsv"))
md <- fread(file.path(tdir, "magma_drug_all_traits.tsv"))
ma <- fread(file.path(tdir, "magma_atc_all_traits.tsv"))

cat("\n============================================================\n")
cat("TWAS-GSEA DRUG: top 10 per trait (P-sorted)\n")
cat("============================================================\n")
for(tr in sort(unique(gd$Trait))){
  cat("\n--- ", tr, " ---\n", sep="")
  print(head(gd[Trait == tr, .(Panel, Name, N_Genes, Estimate, Direction, Reversal_Z, P, P.FDR, ATC_Description)][order(P)], 10))
}

cat("\n\n============================================================\n")
cat("TWAS-GSEA ATC: top 10 per trait (P-sorted), with direction\n")
cat("============================================================\n")
for(tr in sort(unique(ga$Trait))){
  cat("\n--- ", tr, " ---\n", sep="")
  print(head(ga[Trait == tr, .(ATC_Code, ATC_Description, N_Drugs, Estimate, Direction, P, P.FDR)][order(P)], 10))
}

cat("\n\n============================================================\n")
cat("MAGMA DrugTargetor DRUG: top 10 per trait (P-sorted)\n")
cat("============================================================\n")
for(tr in sort(unique(md$Trait))){
  cat("\n--- ", tr, " ---\n", sep="")
  print(head(md[Trait == tr, .(Name, N_Genes, BETA, Z, P, P.FDR, ATC_Description)][order(P)], 10))
}

cat("\n\n============================================================\n")
cat("MAGMA DrugTargetor ATC: top 10 per trait (P-sorted)\n")
cat("============================================================\n")
for(tr in sort(unique(ma$Trait))){
  cat("\n--- ", tr, " ---\n", sep="")
  print(head(ma[Trait == tr, .(ATC_Code, ATC_Description, N_Drugs, Z, P, P.FDR)][order(P)], 10))
}
