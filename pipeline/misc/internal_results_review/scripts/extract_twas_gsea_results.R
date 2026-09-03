#!/usr/bin/env Rscript
# Extract drug-repurposing results from /work/results/results_package.rds:
#   - TWAS-GSEA drug-level and ATC-level (directional + non-directional), per panel
#   - MAGMA DrugTargetor drug-level and ATC-level
#   - GCSC drug-level and ATC-level
#   - CMAP TWAS-GSEA per-signature and per-MOA, per panel
#
# All tables now carry Direction (text) and Reversal_Z (numeric) where the
# underlying test is directional. Reversal_Z > 0 always = opposes disease
# (candidate therapeutic direction); Direction in {Opposes disease, Matches
# disease, NA}. See README in genodisc_patches/ for sign conventions.

suppressMessages({
  library(data.table)
  library(stringr)
})

# Source shared path resolution; works under both Rscript and interactive use.
.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if(length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})
source(file.path(.script_dir, "_paths.R"))
outdir <- tables_dir

x <- readRDS(rds_path)
traits <- setdiff(names(x), "configuration")
cat("Traits:", paste(traits, collapse = ", "), "\n")

clean_name <- function(s) trimws(tolower(s))

############################
# TWAS-GSEA drug-level (directional + non-directional, multi-panel)
############################
collect_drug_gsea <- function(slot, label){
  out <- list()
  for(tr in traits){
    d <- x[[tr]]$tx$drug[[slot]]
    if(is.null(d) || nrow(d) == 0) next
    d <- as.data.table(copy(d))
    setnames(d, c("N Genes","ATC Code","ATC Description"),
                c("N_Genes","ATC_Code","ATC_Description"))
    d[, Trait := tr]
    d[, Result_Type := label]
    d[, Name_clean := clean_name(Name)]
    d[, ChEMBL := NULL]
    setcolorder(d, c("Trait","Result_Type","Panel","Name","Name_clean","N_Genes",
                     "Estimate","SE","P","P.FDR","Direction","Reversal_Z",
                     "ATC_Code","ATC_Description"))
    out[[tr]] <- d
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}
gsea_drug       <- collect_drug_gsea("twas_gsea",        "TWAS-GSEA drug (directional)")
gsea_drug_nondir <- collect_drug_gsea("twas_gsea_nondir", "TWAS-GSEA drug (non-directional)")

fwrite(gsea_drug,        file.path(outdir, "twas_gsea_drug_all_traits.tsv"),        sep = "\t")
fwrite(gsea_drug_nondir, file.path(outdir, "twas_gsea_drug_nondir_all_traits.tsv"), sep = "\t")
cat("twas_gsea_drug rows:",        nrow(gsea_drug),        "\n")
cat("twas_gsea_drug_nondir rows:", nrow(gsea_drug_nondir), "\n")

############################
# TWAS-GSEA ATC-level (directional + non-directional, multi-panel)
############################
collect_atc_gsea <- function(slot, label){
  out <- list()
  for(tr in traits){
    d <- x[[tr]]$tx$atc[[slot]]
    if(is.null(d) || nrow(d) == 0) next
    d <- as.data.table(copy(d))
    setnames(d,
      c("ATC Code","ATC Description","N Drugs","Class Median T","Non-class Median T"),
      c("ATC_Code","ATC_Description","N_Drugs","Class_Median_T","Non_Class_Median_T"))
    d[, Trait := tr]
    d[, Result_Type := label]
    setcolorder(d, c("Trait","Result_Type","Panel","ATC_Code","ATC_Description",
                     "N_Drugs","Estimate","Class_Median_T","Non_Class_Median_T",
                     "P","P.FDR","Direction","Reversal_Z"))
    out[[tr]] <- d
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}
gsea_atc       <- collect_atc_gsea("twas_gsea",        "TWAS-GSEA ATC (directional)")
gsea_atc_nondir <- collect_atc_gsea("twas_gsea_nondir", "TWAS-GSEA ATC (non-directional)")
fwrite(gsea_atc,        file.path(outdir, "twas_gsea_atc_all_traits.tsv"),        sep = "\t")
fwrite(gsea_atc_nondir, file.path(outdir, "twas_gsea_atc_nondir_all_traits.tsv"), sep = "\t")
cat("twas_gsea_atc rows:",        nrow(gsea_atc),        "\n")
cat("twas_gsea_atc_nondir rows:", nrow(gsea_atc_nondir), "\n")

############################
# MAGMA DrugTargetor (drug + ATC, single tissue-agnostic SNP-based)
############################
magma_drug <- rbindlist(lapply(traits, function(tr){
  d <- x[[tr]]$tx$drug$magma
  if(is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.table(copy(d))
  setnames(d, c("N Genes","ATC Code","ATC Description"),
              c("N_Genes","ATC_Code","ATC_Description"))
  d[, Trait := tr]; d[, Result_Type := "MAGMA DrugTargetor drug"]
  d[, Name_clean := clean_name(Name)]
  d[, ChEMBL := NULL]
  d[, Z := -qnorm(P)]
  setcolorder(d, c("Trait","Result_Type","Name","Name_clean","N_Genes",
                   "BETA","SE","Z","P","P.FDR","ATC_Code","ATC_Description"))
  d
}), use.names = TRUE, fill = TRUE)

magma_atc <- rbindlist(lapply(traits, function(tr){
  d <- x[[tr]]$tx$atc$magma
  if(is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.table(copy(d))
  setnames(d, c("ATC Code","N Drugs","ATC Description"),
              c("ATC_Code","N_Drugs","ATC_Description"))
  d[, Trait := tr]; d[, Result_Type := "MAGMA DrugTargetor ATC"]
  d[, Z := -qnorm(P)]
  setcolorder(d, c("Trait","Result_Type","ATC_Code","ATC_Description","N_Drugs","Z","P","P.FDR"))
  d
}), use.names = TRUE, fill = TRUE)

fwrite(magma_drug, file.path(outdir, "magma_drug_all_traits.tsv"), sep = "\t")
fwrite(magma_atc,  file.path(outdir, "magma_atc_all_traits.tsv"),  sep = "\t")
cat("magma_drug rows:", nrow(magma_drug), "\n")
cat("magma_atc rows:",  nrow(magma_atc),  "\n")

############################
# GCSC (drug + ATC, brain+blood aggregated, no panel)
############################
gcsc_drug <- rbindlist(lapply(traits, function(tr){
  d <- x[[tr]]$tx$drug$gcsc
  if(is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.table(copy(d))
  setnames(d, c("ATC Code","ATC Description"),
              c("ATC_Code","ATC_Description"))
  d[, Trait := tr]; d[, Result_Type := "GCSC drug"]
  d[, Name_clean := clean_name(Name)]
  d[, ChEMBL := NULL]
  setcolorder(d, c("Trait","Result_Type","Name","Name_clean","Enrichment","SE","Z","P","P.FDR","ATC_Code","ATC_Description"))
  d
}), use.names = TRUE, fill = TRUE)

gcsc_atc <- rbindlist(lapply(traits, function(tr){
  d <- x[[tr]]$tx$atc$gcsc
  if(is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.table(copy(d))
  setnames(d, c("ATC Code","N Drugs","ATC Description"),
              c("ATC_Code","N_Drugs","ATC_Description"))
  d[, Trait := tr]; d[, Result_Type := "GCSC ATC"]
  d[, Z := -qnorm(P)]
  setcolorder(d, c("Trait","Result_Type","ATC_Code","ATC_Description","N_Drugs","Z","P","P.FDR"))
  d
}), use.names = TRUE, fill = TRUE)
fwrite(gcsc_drug, file.path(outdir, "gcsc_drug_all_traits.tsv"), sep = "\t")
fwrite(gcsc_atc,  file.path(outdir, "gcsc_atc_all_traits.tsv"),  sep = "\t")
cat("gcsc_drug rows:", nrow(gcsc_drug), "\n")
cat("gcsc_atc rows:",  nrow(gcsc_atc),  "\n")

############################
# CMAP TWAS-GSEA (per-signature + per-MOA, multi-panel + multi-cell-line)
############################
cmap_drug <- rbindlist(lapply(traits, function(tr){
  d <- x[[tr]]$tx$cmap$drug
  if(is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.table(copy(d))
  d[, Trait := tr]; d[, Result_Type := "CMAP TWAS-GSEA per-signature"]
  d[, Name_clean := clean_name(cmap_name)]
  setcolorder(d, c("Trait","Result_Type","Panel","Name","Name_clean","cmap_name","cell_iname",
                   "pert_itime","pert_idose","moa","N_Mem_Avail",
                   "Estimate","SE","Z","P","P.FDR","Direction","Reversal_Z"))
  d
}), use.names = TRUE, fill = TRUE)
fwrite(cmap_drug, file.path(outdir, "cmap_drug_all_traits.tsv"), sep = "\t")
cat("cmap_drug rows:", nrow(cmap_drug), "\n")

cmap_moa <- rbindlist(lapply(traits, function(tr){
  d <- x[[tr]]$tx$cmap$moa
  if(is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.table(copy(d))
  setnames(d, c("N Drugs","Class Median T","Non-class Median T"),
              c("N_Drugs","Class_Median_T","Non_Class_Median_T"))
  d[, Trait := tr]; d[, Result_Type := "CMAP TWAS-GSEA per-MOA"]
  setcolorder(d, c("Trait","Result_Type","Panel","MOA","Cell_Line","N_Drugs",
                   "Estimate","Class_Median_T","Non_Class_Median_T","P","P.FDR",
                   "Direction","Reversal_Z"))
  d
}), use.names = TRUE, fill = TRUE)
fwrite(cmap_moa, file.path(outdir, "cmap_moa_all_traits.tsv"), sep = "\t")
cat("cmap_moa rows:", nrow(cmap_moa), "\n")

############################
# Per-trait top-N tables (sorted by P)
############################
top_n <- 30
fwrite(gsea_drug[order(P), .SD[1:min(top_n, .N)], by = Trait],
       file.path(outdir, "twas_gsea_drug_top_by_trait.tsv"), sep = "\t")
fwrite(gsea_atc[order(P), .SD[1:min(top_n, .N)], by = Trait],
       file.path(outdir, "twas_gsea_atc_top_by_trait.tsv"), sep = "\t")
fwrite(cmap_moa[order(P), .SD[1:min(top_n, .N)], by = Trait],
       file.path(outdir, "cmap_moa_top_by_trait.tsv"), sep = "\t")
fwrite(cmap_drug[order(P), .SD[1:min(top_n, .N)], by = Trait],
       file.path(outdir, "cmap_drug_top_by_trait.tsv"), sep = "\t")
fwrite(gcsc_drug[order(P), .SD[1:min(top_n, .N)], by = Trait],
       file.path(outdir, "gcsc_drug_top_by_trait.tsv"), sep = "\t")
fwrite(gcsc_atc[order(P), .SD[1:min(top_n, .N)], by = Trait],
       file.path(outdir, "gcsc_atc_top_by_trait.tsv"), sep = "\t")

############################
# Cross-trait recurrence (per-drug / per-ATC aggregating across panels)
# For each drug, find its BEST panel within each trait (smallest P) and use
# that as the trait-level summary. Then count traits where best-P is below
# nominal/FDR thresholds.
############################

best_per_drug_trait <- gsea_drug[order(P), .SD[1], by = .(Trait, Name_clean)]
best_per_drug_trait[, FDR_Sig := P.FDR < 0.05]
best_per_drug_trait[, Nom_Sig := P < 0.05]
drug_cross <- best_per_drug_trait[, .(
  N_traits_nominal = sum(Nom_Sig, na.rm = TRUE),
  N_traits_fdr     = sum(FDR_Sig, na.rm = TRUE),
  N_opposes        = sum(Direction == "Opposes disease" & Nom_Sig, na.rm = TRUE),
  N_matches        = sum(Direction == "Matches disease" & Nom_Sig, na.rm = TRUE),
  Min_P            = min(P, na.rm = TRUE),
  Min_FDR          = min(P.FDR, na.rm = TRUE),
  Mean_Reversal_Z  = mean(Reversal_Z, na.rm = TRUE),
  Traits_nominal   = paste(unique(Trait[Nom_Sig]), collapse = ","),
  Traits_fdr       = paste(unique(Trait[FDR_Sig]), collapse = ",")
), by = .(Name_clean, ATC_Code, ATC_Description)]
drug_cross <- drug_cross[order(-N_traits_nominal, Min_P)]
fwrite(drug_cross, file.path(outdir, "twas_gsea_drug_cross_trait.tsv"), sep = "\t")

best_per_atc_trait <- gsea_atc[order(P), .SD[1], by = .(Trait, ATC_Code)]
best_per_atc_trait[, FDR_Sig := P.FDR < 0.05]
best_per_atc_trait[, Nom_Sig := P < 0.05]
atc_cross <- best_per_atc_trait[, .(
  N_traits_nominal = sum(Nom_Sig, na.rm = TRUE),
  N_traits_fdr     = sum(FDR_Sig, na.rm = TRUE),
  N_opposes        = sum(Direction == "Opposes disease" & Nom_Sig, na.rm = TRUE),
  N_matches        = sum(Direction == "Matches disease" & Nom_Sig, na.rm = TRUE),
  Min_P            = min(P, na.rm = TRUE),
  Min_FDR          = min(P.FDR, na.rm = TRUE),
  Mean_Reversal_Z  = mean(Reversal_Z, na.rm = TRUE),
  Traits_nominal   = paste(unique(Trait[Nom_Sig]), collapse = ","),
  Traits_fdr       = paste(unique(Trait[FDR_Sig]), collapse = ",")
), by = .(ATC_Code, ATC_Description)]
atc_cross <- atc_cross[order(-N_traits_nominal, Min_P)]
fwrite(atc_cross, file.path(outdir, "twas_gsea_atc_cross_trait.tsv"), sep = "\t")

# CMAP MOA cross-trait
best_per_moa_trait <- cmap_moa[order(P), .SD[1], by = .(Trait, MOA)]
best_per_moa_trait[, FDR_Sig := P.FDR < 0.05]
best_per_moa_trait[, Nom_Sig := P < 0.05]
moa_cross <- best_per_moa_trait[, .(
  N_traits_nominal = sum(Nom_Sig, na.rm = TRUE),
  N_traits_fdr     = sum(FDR_Sig, na.rm = TRUE),
  N_opposes        = sum(Direction == "Opposes disease" & Nom_Sig, na.rm = TRUE),
  N_matches        = sum(Direction == "Matches disease" & Nom_Sig, na.rm = TRUE),
  Min_P            = min(P, na.rm = TRUE),
  Min_FDR          = min(P.FDR, na.rm = TRUE),
  Mean_Reversal_Z  = mean(Reversal_Z, na.rm = TRUE),
  Traits_nominal   = paste(unique(Trait[Nom_Sig]), collapse = ","),
  Traits_fdr       = paste(unique(Trait[FDR_Sig]), collapse = ",")
), by = .(MOA)]
moa_cross <- moa_cross[order(-N_traits_nominal, Min_P)]
fwrite(moa_cross, file.path(outdir, "cmap_moa_cross_trait.tsv"), sep = "\t")

############################
# Method-comparison tables
############################
gd <- best_per_drug_trait[, .(Trait, Drug = Name_clean,
                              gsea_P = P, gsea_FDR = P.FDR, gsea_Reversal_Z = Reversal_Z,
                              gsea_Direction = Direction, gsea_Panel = Panel)]
md <- magma_drug[, .(Trait, Drug = Name_clean, magma_P = P, magma_FDR = P.FDR, magma_Z = Z)]
gc_d <- gcsc_drug[, .(Trait, Drug = Name_clean, gcsc_P = P, gcsc_FDR = P.FDR, gcsc_Z = Z)]
drug_cmp <- Reduce(function(a,b) merge(a,b, by = c("Trait","Drug"), all = TRUE),
                   list(gd, md, gc_d))
drug_cmp[, gsea_FDR_sig  := !is.na(gsea_FDR)  & gsea_FDR  < 0.05]
drug_cmp[, magma_FDR_sig := !is.na(magma_FDR) & magma_FDR < 0.05]
drug_cmp[, gcsc_FDR_sig  := !is.na(gcsc_FDR)  & gcsc_FDR  < 0.05]
fwrite(drug_cmp, file.path(outdir, "method_comparison_drug.tsv"), sep = "\t")

ga <- best_per_atc_trait[, .(Trait, ATC_Code, ATC_Description,
                             gsea_P = P, gsea_FDR = P.FDR, gsea_Reversal_Z = Reversal_Z,
                             gsea_Direction = Direction, gsea_Panel = Panel)]
ma <- magma_atc[, .(Trait, ATC_Code, magma_P = P, magma_FDR = P.FDR, magma_Z = Z)]
gc_a <- gcsc_atc[, .(Trait, ATC_Code, gcsc_P = P, gcsc_FDR = P.FDR, gcsc_Z = Z)]
atc_cmp <- Reduce(function(a,b) merge(a,b, by = c("Trait","ATC_Code"), all = TRUE),
                  list(ga, ma, gc_a))
atc_cmp[, gsea_FDR_sig  := !is.na(gsea_FDR)  & gsea_FDR  < 0.05]
atc_cmp[, magma_FDR_sig := !is.na(magma_FDR) & magma_FDR < 0.05]
atc_cmp[, gcsc_FDR_sig  := !is.na(gcsc_FDR)  & gcsc_FDR  < 0.05]
fwrite(atc_cmp, file.path(outdir, "method_comparison_atc.tsv"), sep = "\t")

############################
# Coverage summary by method
############################
summary_tbl <- rbindlist(list(
  gsea_drug      [, .(Method = "TWAS-GSEA drug (dir)",   Nominal = sum(P < 0.05, na.rm = TRUE),   FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)],
  gsea_drug_nondir[, .(Method = "TWAS-GSEA drug (nondir)", Nominal = sum(P < 0.05, na.rm = TRUE), FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)],
  gsea_atc       [, .(Method = "TWAS-GSEA ATC (dir)",    Nominal = sum(P < 0.05, na.rm = TRUE),   FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)],
  magma_drug     [, .(Method = "MAGMA drug",              Nominal = sum(P < 0.05, na.rm = TRUE),  FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)],
  magma_atc      [, .(Method = "MAGMA ATC",               Nominal = sum(P < 0.05, na.rm = TRUE),  FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)],
  gcsc_drug      [, .(Method = "GCSC drug",               Nominal = sum(P < 0.05, na.rm = TRUE),  FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)],
  gcsc_atc       [, .(Method = "GCSC ATC",                Nominal = sum(P < 0.05, na.rm = TRUE),  FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)],
  cmap_drug      [, .(Method = "CMAP drug (per-sig)",     Nominal = sum(P < 0.05, na.rm = TRUE),  FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)],
  cmap_moa       [, .(Method = "CMAP per-MOA",            Nominal = sum(P < 0.05, na.rm = TRUE),  FDR_sig = sum(P.FDR < 0.05, na.rm = TRUE)), by = .(Trait)]
))
fwrite(summary_tbl, file.path(outdir, "hits_per_trait_per_method.tsv"), sep = "\t")

cat("\nDONE\n")
cat("FDR-significant hits per trait by method:\n")
print(dcast(summary_tbl, Trait ~ Method, value.var = "FDR_sig"))
