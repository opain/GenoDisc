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

out_cols <- c("locus", "cs", "ensgid", "gene", "caldera", "multi", "rank",
              "nearest_gene", "n_genes", "dist", "pops", "coding", "impute_pops")
write_empty <- function() fwrite(setNames(data.table(matrix(character(0), ncol = length(out_cols))), out_cols), out_path)

cv <- tryCatch(fread(paste0(fm_dir, opt$gwas, ".GW.credset_variants.csv")), error = function(e) NULL)
if (is.null(cv) || nrow(cv) == 0) { write_empty(); quit(save = "no") }

# CALDERA takes one credible set per "locus". Use the L=10 (multi-signal) sets
# where a locus has them, but fall back to its L=1 (single-signal) set otherwise,
# so EVERY fine-mapped locus is covered — not just the few that produced an L=10
# decomposition (previously `cv[L == 10]` dropped all L=1-only loci). Make a
# distinct locus id per (lead, cs).
cv[, .has10 := any(L == 10), by = locus]
cv <- cv[(.has10 & L == 10) | (!.has10 & L == 1)]
cv[, .has10 := NULL]
if (nrow(cv) == 0) { write_empty(); quit(save = "no") }
cv[, locus_id := paste0(locus, "_cs", cs)]
cs_dt <- cv[, .(locus = locus_id, chr = CHR, bp = BP, pip = PIP, rsid = SNP)]
cs_file <- tempfile(fileext = ".tsv")
fwrite(cs_dt, cs_file, sep = "\t")

# Run CALDERA with the shipped pretrained model (chdir so its relative source()s
# and data/ + trained_models/ paths resolve).
source(file.path(opt$caldera_path, "z_caldera.R"), chdir = TRUE)
res <- caldera(pops_file = opt$pops_preds, cs_file = cs_file, assembly = 37,
               caldera_path = opt$caldera_path)
res <- as.data.table(res)
if (nrow(res) == 0) { write_empty(); quit(save = "no") }

# Recover lead + cs from the composite locus id.
res[, cs := as.integer(sub(".*_cs", "", locus))]
res[, locus := sub("_cs[0-9]+$", "", locus)]

# Rank by CALDERA probability within each credible set; nearest-gene baseline.
setorder(res, locus, cs, -caldera)
res[, rank := seq_len(.N), by = .(locus, cs)]
res[, nearest_gene := dist == min(dist, na.rm = TRUE), by = .(locus, cs)]

for (cc in setdiff(out_cols, names(res))) res[, (cc) := NA]
res <- res[, ..out_cols]

# Join the Phase 1 coding-annotation per-cs summary when it exists (display only).
coding_summ <- paste0(fm_dir, opt$gwas, ".GW.credset_coding_summary.csv")
if (file.exists(coding_summ)) {
  cs_sum <- tryCatch(fread(coding_summ), error = function(e) NULL)
  if (!is.null(cs_sum) && nrow(cs_sum) > 0) {
    keep <- intersect(c("locus", "cs", "pip_protein_altering", "frac_protein_altering",
                        "top_coding_gene", "top_coding_consequence", "top_coding_am_class",
                        "flag_protein_altering"), names(cs_sum))
    res <- merge(res, cs_sum[, ..keep], by = c("locus", "cs"), all.x = TRUE, sort = FALSE)
  }
}

setorder(res, locus, cs, rank)
fwrite(res, out_path)
