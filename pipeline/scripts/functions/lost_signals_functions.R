# Core logic for the "lost signal" QC summary (see scripts/summarise_lost_signals.R).
# Kept as a pure, I/O-free function so it can be unit-tested (tests/testthat/test-lost-signals.R).

suppressMessages(library(data.table))

# Group removed genome-wide-significant variants into approximate loci (greedy ±half_win windows
# around each lead/lowest-P removed variant) and, for each locus, count retained variants nearby and
# flag loci whose entire GW-significant signal was removed.
#
# IMPORTANT: removed_sig and retained MUST carry positions in the SAME genome build. The cleaned
# sumstats are remapped to GRCh37 by the GenoUtils cleaner, whereas the removed/retained sidecars are
# in the original user/input build; this function must be fed the two user-build sidecars so the
# positional window comparison is valid for non-GRCh37 inputs (otherwise nearby retained signals are
# missed and loci are falsely flagged lost).
#
# removed_sig: data.table with CHR, BP, P, and optionally SNP, reason (removed GW-sig variants).
# retained:    data.table with CHR, BP, and optionally P (retained variants, user build).
# Returns a data.table of loci (one row per lead removed variant), ordered by lead_P.
compute_lost_signal_loci <- function(removed_sig, retained, gw_thresh = 5e-8, half_win = 5e5) {
  empty <- data.table(CHR = character(0), lead_SNP = character(0), lead_BP = numeric(0),
                      lead_P = numeric(0), reason = character(0), n_removed = integer(0),
                      n_retained_in_window = integer(0), n_retained_sig_in_window = integer(0),
                      lost_signal = logical(0))
  if (is.null(removed_sig) || nrow(removed_sig) == 0) return(empty)

  rs <- copy(as.data.table(removed_sig))
  rs[, CHR := as.character(CHR)]
  rs[, BP := as.numeric(BP)]
  rs <- rs[order(P)]            # lowest-P first => lead of its locus
  rs[, assigned := FALSE]

  cleaned <- copy(as.data.table(retained))
  if (nrow(cleaned) > 0) cleaned[, CHR := as.character(CHR)]

  loci_list <- list()
  for (i in seq_len(nrow(rs))) {
    if (rs$assigned[i]) next
    lead <- rs[i]
    in_win <- rs$CHR == lead$CHR & abs(rs$BP - lead$BP) <= half_win & !rs$assigned
    rs[in_win, assigned := TRUE]

    cl_win <- if (nrow(cleaned) > 0) cleaned[CHR == lead$CHR & abs(BP - lead$BP) <= half_win] else cleaned
    n_ret      <- nrow(cl_win)
    n_ret_sig  <- if ('P' %in% names(cl_win)) sum(!is.na(cl_win$P) & cl_win$P < gw_thresh) else 0L

    loci_list[[length(loci_list) + 1]] <- data.table(
      CHR = lead$CHR,
      lead_SNP = if ('SNP' %in% names(lead)) as.character(lead$SNP) else NA_character_,
      lead_BP = lead$BP,
      lead_P = lead$P,
      reason = if ('reason' %in% names(lead)) lead$reason else NA_character_,
      n_removed = sum(in_win),
      n_retained_in_window = n_ret,
      n_retained_sig_in_window = n_ret_sig,
      # Lost signal = no retained variant in the window reaches GW significance, i.e. the
      # GW-significant signal at this locus is not represented in the cleaned data.
      lost_signal = (n_ret_sig == 0))
  }
  loci <- rbindlist(loci_list)
  loci[order(lead_P)]
}
