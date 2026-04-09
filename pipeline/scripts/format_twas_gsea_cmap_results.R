#!/usr/bin/Rscript

# Format the per-(gwas, weight) CMAP TWAS-GSEA output:
#  - emit a per-signature CSV with parsed cmap_name / cell_iname / pert_itime /
#    pert_idose / moa columns
#  - emit a per-MOA enrichment CSV using the same Wilcoxon framework as
#    format_twas_gsea_drugtargetor_results.R, with MOA in place of ATC

suppressMessages({
  library(optparse)
  library(data.table)
  library(stringr)
})

option_list <- list(
  make_option("--twas",        action = "store", default = NA, type = 'character', help = "GWAS ID [required]"),
  make_option("--panel",       action = "store", default = NA, type = 'character', help = "Weight panel [required]"),
  make_option("--config_file", action = "store", default = NA, type = 'character', help = "Path to config file [required]")
)
opt <- parse_args(OptionParser(option_list = option_list))

source('scripts/functions/utils_functions.R')

config <- readLines(opt$config_file)
outdir <- gsub('outdir: ', '', config[grepl('outdir: ', config)])
cmap_compoundinfo <- read_param(config = opt$config_file, param = 'cmap_compoundinfo', return_obj = FALSE)
if(is.na(cmap_compoundinfo) || cmap_compoundinfo == 'NA') {
  stop("cmap_compoundinfo path must be set in the config file to format CMAP TWAS-GSEA results.")
}

competitive_path <- paste0(outdir, '/results/', opt$twas, '/twas/cmap/twas_gsea_cmap_', opt$panel, '.competitive.txt')
res <- fread(competitive_path)

# T is signed: positive = drug mimics disease signature, negative = drug
# reverses it (the repurposing direction). Use a two-sided P, matching the
# directional DrugTargetor formatter convention.
res$P      <- 2 * pnorm(-abs(res$T))
res$P.FDR  <- p.adjust(res$P, method = 'fdr')
res$Z      <- res$Estimate / res$SE
res$Panel  <- opt$panel

# Parse GeneSet of the form "{idx}.{cmap_name}_{cell_iname}_{pert_itime}_{pert_idose}".
# cmap_name can itself contain hyphens, but the prop generator stripped spaces
# and joined fields with underscores; the trailing 3 underscore-separated tokens
# are always (cell, time, dose).
parse_geneset <- function(g) {
  # Strip leading "<digits>." prefix added by the prop generator's cbind step.
  g_clean <- str_remove(g, '^\\d+\\.')
  parts <- strsplit(g_clean, '_', fixed = TRUE)[[1]]
  n <- length(parts)
  if(n < 4) return(c(g_clean, NA_character_, NA_character_, NA_character_))
  cmap_name <- paste(parts[seq_len(n - 3)], collapse = '_')
  cell      <- parts[n - 2]
  ttime     <- parts[n - 1]
  dose      <- parts[n]
  c(cmap_name, cell, ttime, dose)
}
parsed <- do.call(rbind, lapply(res$GeneSet, parse_geneset))
res$cmap_name  <- parsed[, 1]
res$cell_iname <- parsed[, 2]
res$pert_itime <- parsed[, 3]
res$pert_idose <- parsed[, 4]

# Attach MOA from compoundinfo_beta.txt (case-insensitive cmap_name match).
cmpd <- fread(cmap_compoundinfo)
if(!('cmap_name' %in% names(cmpd)) || !('moa' %in% names(cmpd))) {
  stop("cmap_compoundinfo must contain cmap_name and moa columns.")
}
cmpd <- unique(cmpd[, .(cmap_name_lc = tolower(cmap_name), moa)])
res$moa <- cmpd$moa[match(tolower(res$cmap_name), cmpd$cmap_name_lc)]

# ---- Per-signature CSV --------------------------------------------------
drug_out <- res[, .(cmap_name, cell_iname, pert_itime, pert_idose, moa,
                    N_Mem_Avail, Estimate, SE, Z, P, P.FDR, Panel)]
setorder(drug_out, P)
drug_out_path <- paste0(outdir, '/results/', opt$twas, '/twas/cmap/twas_gsea_cmap_', opt$panel, '_drug_res.csv')
fwrite(drug_out, drug_out_path)
message('Wrote per-signature CSV: ', drug_out_path)

# ---- Per-MOA enrichment CSV (per cell line) ------------------------------
res_moa <- res[!is.na(moa) & moa != '', ]

cell_lines <- unique(res_moa$cell_iname)
moa_enrich_all <- list()
idx <- 0L

for(cl in cell_lines) {
  res_cl   <- res_moa[cell_iname == cl, ]
  ranked_T <- rank(res_cl$T)
  moa_terms <- unique(res_cl$moa)

  for(m in moa_terms) {
    in_idx  <- which(res_cl$moa == m)
    if(length(unique(res_cl$cmap_name[in_idx])) < 3) next
    out_idx <- setdiff(seq_len(nrow(res_cl)), in_idx)
    wil <- tryCatch(
      wilcox.test(ranked_T[in_idx], ranked_T[out_idx], conf.int = TRUE, alternative = 'two.sided'),
      error = function(e) NULL
    )
    if(is.null(wil)) next
    idx <- idx + 1L
    moa_enrich_all[[idx]] <- data.frame(
      MOA              = m,
      Cell_Line        = cl,
      N                = length(unique(res_cl$cmap_name[in_idx])),
      Estimate         = as.numeric(wil$estimate),
      Class_Median     = median(res_cl$Estimate[in_idx]),
      Non_Class_Median = median(res_cl$Estimate[out_idx]),
      P                = wil$p.value,
      Panel            = opt$panel,
      stringsAsFactors = FALSE
    )
  }
}

moa_enrich <- rbindlist(moa_enrich_all[seq_len(idx)])
if(nrow(moa_enrich) > 0) {
  moa_enrich[, P.FDR := p.adjust(P, method = 'fdr')]
  setorder(moa_enrich, P)
}
moa_out_path <- paste0(outdir, '/results/', opt$twas, '/twas/cmap/twas_gsea_cmap_', opt$panel, '_moa_res.csv')
fwrite(moa_enrich, moa_out_path)
message('Wrote per-MOA CSV: ', moa_out_path, ' (', nrow(moa_enrich), ' MOA x cell-line rows)')
