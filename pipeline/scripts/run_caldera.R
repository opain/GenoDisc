#!/usr/bin/Rscript
# CALDERA effector-gene probabilities per SuSiE credible set (Heilbron et al. 2024).
# Uses the published PRETRAINED model shipped in the repo (no retraining). Inputs:
#   - PoPS .preds  (ENSGID, PoPS_Score)                  <- Phase 2
#   - per-variant credible-set table (locus, cs, SNP, CHR, BP, PIP) <- Phase 0.5
# CALDERA's own bundled coding-SNP list drives its non-synonymous-PIP feature
# (matched by rsid), so this is independent of the Phase 1 AlphaMissense output;
# the Phase 1 per-cs summary is joined on only for display when available.
#
# caldera(pops_file, cs_file, assembly=37, caldera_path) returns one row per gene
# per locus: locus, locus_pos, gene, caldera (prob 0-1), multi, n_genes, dist,
# pops, coding, impute_pops, ensgid. We add cs, rank, and a nearest-gene baseline.

suppressMessages(library("optparse"))
option_list <- list(
  make_option("--gwas", default = NA, type = "character"),
  make_option("--config_file", default = NA, type = "character"),
  make_option("--pipeline_dir", default = NA, type = "character"),
  make_option("--caldera_path", default = NA, type = "character"),
  make_option("--pops_preds", default = NA, type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
options(pipeline_dir = opt$pipeline_dir)
suppressMessages(library(data.table))
source(file.path(opt$pipeline_dir, "scripts", "functions", "utils_functions.R"))

outdir <- read_param(config = opt$config_file, param = "outdir", return_obj = FALSE)
fm_dir <- paste0(outdir, "/results/", opt$gwas, "/finemap/")
cal_dir <- paste0(outdir, "/results/", opt$gwas, "/caldera/")
dir.create(cal_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- paste0(cal_dir, opt$gwas, ".GW.caldera.csv")

out_cols <- c("locus", "cs", "L", "ensgid", "gene", "caldera", "multi", "rank",
              "nearest_gene", "n_genes", "dist", "pops", "coding", "impute_pops")
write_empty <- function() fwrite(setNames(data.table(matrix(character(0), ncol = length(out_cols))), out_cols), out_path)

cv <- tryCatch(fread(paste0(fm_dir, opt$gwas, ".GW.credset_variants.csv")), error = function(e) NULL)
if (is.null(cv) || nrow(cv) == 0) { write_empty(); quit(save = "no") }

# Run CALDERA SEPARATELY for the L=1 (single-signal) and L=10 (multi-signal)
# SuSiE credible sets, tagging each result with its L so the viewer can switch
# between them. L=1 covers every fine-mapped locus; L=10 covers only loci that
# produced a multi-signal decomposition. Source the shipped model once (chdir so
# its relative source()s + data/ + trained_models/ paths resolve).
source(file.path(opt$caldera_path, "z_caldera.R"), chdir = TRUE)

run_one <- function(cv_L, Lval) {
  if (is.null(cv_L) || nrow(cv_L) == 0) return(NULL)
  cv_L <- data.table::copy(cv_L)
  cv_L[, locus_id := paste0(locus, "_cs", cs)]
  cs_dt <- cv_L[, .(locus = locus_id, chr = CHR, bp = BP, pip = PIP, rsid = SNP)]
  cs_file <- tempfile(fileext = ".tsv")
  fwrite(cs_dt, cs_file, sep = "\t")
  res <- tryCatch(caldera(pops_file = opt$pops_preds, cs_file = cs_file,
                          assembly = 37, caldera_path = opt$caldera_path),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res <- as.data.table(res)
  res[, cs := as.integer(sub(".*_cs", "", locus))]
  res[, locus := sub("_cs[0-9]+$", "", locus)]
  res[, L := Lval]
  res
}

res <- data.table::rbindlist(
  list(run_one(cv[L == 1], 1L), run_one(cv[L == 10], 10L)),
  use.names = TRUE, fill = TRUE)
if (is.null(res) || nrow(res) == 0) { write_empty(); quit(save = "no") }

# Rank by CALDERA probability within each credible set (per L); nearest-gene baseline.
setorder(res, L, locus, cs, -caldera)
res[, rank := seq_len(.N), by = .(L, locus, cs)]
res[, nearest_gene := dist == min(dist, na.rm = TRUE), by = .(L, locus, cs)]

for (cc in setdiff(out_cols, names(res))) res[, (cc) := NA]
res <- res[, ..out_cols]

# Join the Phase 1 coding-annotation per-cs summary when it exists (display only).
coding_summ <- paste0(fm_dir, opt$gwas, ".GW.credset_coding_summary.csv")
if (file.exists(coding_summ)) {
  cs_sum <- tryCatch(fread(coding_summ), error = function(e) NULL)
  if (!is.null(cs_sum) && nrow(cs_sum) > 0) {
    keep <- intersect(c("locus", "cs", "L", "pip_protein_altering", "frac_protein_altering",
                        "top_coding_gene", "top_coding_consequence", "top_coding_am_class",
                        "flag_protein_altering"), names(cs_sum))
    # Key on L too when the coding summary carries it (credible sets differ by L),
    # else rows cross-join across the two L values.
    merge_key <- intersect(c("locus", "cs", "L"), keep)
    res <- merge(res, cs_sum[, ..keep], by = merge_key, all.x = TRUE, sort = FALSE)
  }
}

setorder(res, L, locus, cs, rank)
fwrite(res, out_path)
