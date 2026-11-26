#!/usr/bin/Rscript
# Summarise genome-wide significant variants removed during GWAS sumstat QC ("lost signals").
#
# Reads the two per-variant sidecars written by the GenoUtils cleaner, both in the original
# user/input build: the removed-variant file (<gwas>.cleaned.removed.gz: SNP, CHR, BP, P, reason) and
# the retained-variant file (<gwas>.cleaned.retained.gz: SNP, CHR, BP, P). It produces
# <gwas>.lost_signals.rds containing:
#   - counts of removed GW-significant (P < 5e-8) variants by reason,
#   - approximate loci (greedy 1 Mb windows around each lead/lowest-P removed variant),
#   - for each locus, whether any variant was RETAINED within the window (and whether a retained
#     variant there is itself GW-significant), flagging loci whose entire GW-significant signal was
#     removed.
# Positions are compared on a single build: the cleaned sumstats (<gwas>.cleaned.gz) are remapped to
# GRCh37 by the cleaner, so we use the retained sidecar (user build) to match the removed sidecar's
# build -- otherwise a non-GRCh37 input mislabels nearby retained signals as lost. For bundles made
# before the retained sidecar existed we fall back to cleaned.gz (GRCh37) and flag build_reconciled
# = FALSE. No LD reference is used; loci are purely positional. This script is additive and reads only.

suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

suppressMessages(library(data.table))
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'lost_signals_functions.R'))

outdir <- read_param(config = opt$config_file, param = 'outdir', return_obj = F)
removed_p_thresh <- suppressWarnings(as.numeric(read_param(config = opt$config_file, param = 'removed_p_thresh', return_obj = F)))
if (is.na(removed_p_thresh)) removed_p_thresh <- 1e-5

GW_THRESH <- 5e-8      # genome-wide significance
WINDOW_BP <- 1e6       # locus window width (lead +/- WINDOW_BP/2)
HALF_WIN  <- WINDOW_BP / 2

ss_dir   <- file.path(outdir, 'results', opt$gwas, 'gwas_sumstat')
rem_path <- file.path(ss_dir, paste0(opt$gwas, '.cleaned.removed.gz'))
ret_path <- file.path(ss_dir, paste0(opt$gwas, '.cleaned.retained.gz'))
cln_path <- file.path(ss_dir, paste0(opt$gwas, '.cleaned.gz'))
out_path <- file.path(ss_dir, paste0(opt$gwas, '.lost_signals.rds'))

removed <- fread(rem_path)
has_coords <- all(c('CHR','BP') %in% names(removed)) && 'P' %in% names(removed)

# GW-significant removed variants (removed.gz is already pre-filtered to P < removed_p_thresh,
# which contains all P < 5e-8 as long as removed_p_thresh >= 5e-8).
removed_sig <- if ('P' %in% names(removed)) removed[!is.na(P) & P < GW_THRESH] else removed[0]

# Counts by reason among removed GW-significant variants.
counts_by_reason <- if (nrow(removed_sig) > 0) {
  tmp <- removed_sig[, .N, by = reason][order(-N)]
  setnames(tmp, c('reason', 'n'))
  tmp
} else {
  data.table(reason = character(0), n = integer(0))
}

# Group removed GW-significant variants into approximate loci and compare to retained variants.
# build_reconciled: TRUE when removed and retained were compared on a common (user) build via the
# retained sidecar; FALSE when we fell back to cleaned.gz (GRCh37) for a pre-sidecar bundle; NA when
# the check did not run (no coordinates / no removed GW-sig variants).
loci <- compute_lost_signal_loci(NULL, NULL)   # empty default
build_reconciled <- NA

if (has_coords && nrow(removed_sig) > 0) {
  # Retained variants in the SAME build as the removed file. cleaned.gz has BP remapped to GRCh37, so
  # prefer the user-build retained sidecar; fall back to cleaned.gz for bundles made before it existed.
  if (file.exists(ret_path)) {
    ret_cols <- names(fread(ret_path, nrows = 0))
    retained <- fread(ret_path, select = intersect(c('CHR','BP','SNP','P'), ret_cols))
    build_reconciled <- TRUE
  } else {
    cln_cols <- names(fread(cln_path, nrows = 0))
    retained <- fread(cln_path, select = intersect(c('CHR','BP','SNP','P'), cln_cols))
    build_reconciled <- FALSE
  }
  loci <- compute_lost_signal_loci(removed_sig, retained, GW_THRESH, HALF_WIN)
}

lost_signals <- list(
  gwas                 = opt$gwas,
  gw_threshold         = GW_THRESH,
  removed_p_thresh     = removed_p_thresh,
  window_bp            = WINDOW_BP,
  has_coords           = has_coords,
  build_reconciled     = build_reconciled,
  n_removed_sig        = nrow(removed_sig),
  counts_by_reason     = counts_by_reason,
  loci                 = loci,
  n_lost_signal_loci   = if (nrow(loci) > 0) sum(loci$lost_signal) else 0L
)

saveRDS(lost_signals, out_path)

cat(paste0('lost_signals: ', nrow(removed_sig), ' removed variants with P < ', GW_THRESH,
           ' across ', nrow(loci), ' loci (', lost_signals$n_lost_signal_loci,
           ' with no retained GW-significant variant). Written to ', out_path, '\n'))
