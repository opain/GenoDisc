#!/usr/bin/Rscript
# Summarise genome-wide significant variants removed during GWAS sumstat QC ("lost signals").
#
# Reads the per-variant removed-variant file written by the GenoUtils cleaner
# (<gwas>.cleaned.removed.gz: original user-build SNP, CHR, BP, P, reason) and the cleaned
# sumstats (<gwas>.cleaned.gz), and produces <gwas>.lost_signals.rds containing:
#   - counts of removed GW-significant (P < 5e-8) variants by reason,
#   - approximate loci (greedy 1 Mb windows around each lead/lowest-P removed variant),
#   - for each locus, whether any variant was RETAINED in the cleaned data within the window
#     (and whether a retained variant there is itself GW-significant), flagging loci whose
#     entire GW-significant signal was removed.
# No LD reference is used; loci are purely positional. This script is additive and reads only.

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

outdir <- read_param(config = opt$config_file, param = 'outdir', return_obj = F)
removed_p_thresh <- suppressWarnings(as.numeric(read_param(config = opt$config_file, param = 'removed_p_thresh', return_obj = F)))
if (is.na(removed_p_thresh)) removed_p_thresh <- 1e-5

GW_THRESH <- 5e-8      # genome-wide significance
WINDOW_BP <- 1e6       # locus window width (lead +/- WINDOW_BP/2)
HALF_WIN  <- WINDOW_BP / 2

ss_dir   <- file.path(outdir, 'results', opt$gwas, 'gwas_sumstat')
rem_path <- file.path(ss_dir, paste0(opt$gwas, '.cleaned.removed.gz'))
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
loci <- data.table(CHR = character(0), lead_SNP = character(0), lead_BP = numeric(0),
                   lead_P = numeric(0), reason = character(0), n_removed = integer(0),
                   n_retained_in_window = integer(0), n_retained_sig_in_window = integer(0),
                   lost_signal = logical(0))

if (has_coords && nrow(removed_sig) > 0) {
  # Retained (cleaned) variant positions, for the "is anything still here?" check.
  cln_cols <- names(fread(cln_path, nrows = 0))
  cleaned <- fread(cln_path, select = intersect(c('CHR','BP','SNP','P'), cln_cols))
  cleaned[, CHR := as.character(CHR)]

  rs <- copy(removed_sig)
  rs[, CHR := as.character(CHR)]
  rs[, BP := as.numeric(BP)]
  rs <- rs[order(P)]            # lowest-P first => lead of its locus
  rs[, assigned := FALSE]

  loci_list <- list()
  for (i in seq_len(nrow(rs))) {
    if (rs$assigned[i]) next
    lead <- rs[i]
    in_win <- rs$CHR == lead$CHR & abs(rs$BP - lead$BP) <= HALF_WIN & !rs$assigned
    rs[in_win, assigned := TRUE]

    cl_win <- cleaned[CHR == lead$CHR & abs(BP - lead$BP) <= HALF_WIN]
    n_ret      <- nrow(cl_win)
    n_ret_sig  <- if ('P' %in% names(cl_win)) sum(!is.na(cl_win$P) & cl_win$P < GW_THRESH) else 0L

    loci_list[[length(loci_list) + 1]] <- data.table(
      CHR = lead$CHR,
      lead_SNP = if ('SNP' %in% names(lead)) as.character(lead$SNP) else NA_character_,
      lead_BP = lead$BP,
      lead_P = lead$P,
      reason = lead$reason,
      n_removed = sum(in_win),
      n_retained_in_window = n_ret,
      n_retained_sig_in_window = n_ret_sig,
      # Lost signal = no retained variant in the window reaches GW significance, i.e. the
      # GW-significant signal at this locus is not represented in the cleaned data.
      lost_signal = (n_ret_sig == 0))
  }
  loci <- rbindlist(loci_list)
  loci <- loci[order(lead_P)]
}

lost_signals <- list(
  gwas                 = opt$gwas,
  gw_threshold         = GW_THRESH,
  removed_p_thresh     = removed_p_thresh,
  window_bp            = WINDOW_BP,
  has_coords           = has_coords,
  n_removed_sig        = nrow(removed_sig),
  counts_by_reason     = counts_by_reason,
  loci                 = loci,
  n_lost_signal_loci   = if (nrow(loci) > 0) sum(loci$lost_signal) else 0L
)

saveRDS(lost_signals, out_path)

cat(paste0('lost_signals: ', nrow(removed_sig), ' removed variants with P < ', GW_THRESH,
           ' across ', nrow(loci), ' loci (', lost_signals$n_lost_signal_loci,
           ' with no retained GW-significant variant). Written to ', out_path, '\n'))
