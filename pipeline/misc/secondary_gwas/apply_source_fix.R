#!/usr/bin/env Rscript
# Apply a known-source pre-fix to a raw sumstats file, writing the fixed
# copy to a new path. Called from munge_one.sh only for codes listed in
# fixes.csv.
#
# Usage:
#   Rscript apply_source_fix.R <src.gz> <out.gz> <fix>
#
# Supported fixes:
#   drop_freq  Strip the FREQ/FRQ/EAF/MAF column so sumstat_cleaner
#              skips its allele-freq-vs-reference discord filter. Used
#              when the source's FREQ column has been silently swapped
#              (freq of A2, or freq of ref-A1) and can't be trusted.
#              Safe for LDSC purposes: LDSC's --rg only reads
#              SNP/A1/A2/Z/N and ignores FREQ.
#   split_snp  Split a "CHR:POS" or "CHR:POS:INDEL" SNP column into
#              separate CHR + BP integer columns and drop indels; used
#              when the source has neither CHR/BP columns nor rs IDs
#              (ref_harmonise then has nothing to merge on). Drops the
#              chr:pos SNP so ref_harmonise regenerates canonical rs
#              IDs from the reference during CHR/BP merge.
#
# Adding a new fix: add a case below, extend fixes.csv, document it.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
  stop("usage: apply_source_fix.R <src.gz> <out.gz> <fix>")
src <- args[1]; out <- args[2]; fix <- args[3]

suppressMessages(library(data.table))
d <- fread(src)
cat("input rows: ", nrow(d), "  cols: ", paste(names(d), collapse = ","),
    "\n", sep = "")

if (fix == "drop_freq") {
  freq_cols <- grep("^(FREQ|FRQ|EAF|MAF|A1_FREQ|A1_FRQ)$",
                    names(d), value = TRUE)
  if (length(freq_cols) == 0L)
    stop("drop_freq: no FREQ column found; nothing to strip.")
  d[, (freq_cols) := NULL]
  cat("dropped: ", paste(freq_cols, collapse = ","), "\n", sep = "")

} else if (fix == "split_snp") {
  if (!"SNP" %in% names(d))
    stop("split_snp: no SNP column present.")
  parts <- tstrsplit(d$SNP, ":", fixed = TRUE)
  d[, CHR := as.integer(parts[[1]])]
  d[, BP  := as.integer(parts[[2]])]
  # Drop indels (rows with a 3rd :INDEL segment).
  tag <- if (length(parts) >= 3L) parts[[3]] else rep(NA_character_, nrow(d))
  before <- nrow(d)
  d <- d[(is.na(tag) | tag == "") & !is.na(CHR) & !is.na(BP) &
         CHR >= 1L & CHR <= 22L]
  cat("kept ", nrow(d), " / ", before,
      " after indel + autosome filter\n", sep = "")
  # Drop the CHR:POS SNP; ref_harmonise will attach canonical rs IDs
  # from the reference during the CHR/BP merge.
  d[, SNP := NULL]

} else {
  stop("unknown fix: ", fix)
}

cat("output rows: ", nrow(d), "  cols: ", paste(names(d), collapse = ","),
    "\n", sep = "")
fwrite(d, out, sep = "\t", na = "NA", quote = FALSE)
