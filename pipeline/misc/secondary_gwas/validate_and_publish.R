#!/usr/bin/env Rscript
# Validate staged .sumstats.gz files, publish passing ones to the live
# secondary_gwas dir, and rebuild sumstat_repo.csv with a slim schema.
#
# Inputs (in this dir):
#   targets.csv           — driver table (from build_targets.R)
#   targets_excluded.csv  — codes excluded pre-munge (non-EUR / no source)
#   staging/<code>.sumstats.gz — outputs from munge_one.sh
#   logs/munge-<code>.out — per-code worker log
#
# Outputs:
#   rebuild_log.csv                                  — status per input code
#   <SGWAS>/<code>.sumstats.gz                       — published files
#   <SGWAS>/sumstat_repo.csv                         — rebuilt (slim) catalog
#   <SGWAS>/sumstat_repo.csv.bak-preRebuild-<date>   — backup of prior CSV
#
# Validation: header must be CHR SNP BP A1 A2 Z N (same as primary GWAS's
# .cleaned.munged.mergealleles.sumstats.gz), rows >= MIN_SNPS, all SNPs
# subset of w_hm3.snplist.
suppressMessages(library(data.table))

MISC   <- "/scratch/prj/neurohackpain/GenoDisc/repo/current/pipeline/misc/secondary_gwas"
SGWAS  <- "/scratch/prj/neurohackpain/GenoDisc/pipeline_resources/data/secondary_gwas"
WHM3   <- "/scratch/prj/neurohackpain/GenoDisc/pipeline_resources/data/ldsc/w_hm3.snplist"
# Staging lives under SGWAS (not MISC) to keep the ~10GB intermediate
# working directory out of the repo tree.
STAGING <- file.path(SGWAS, "staging")
MIN_SNPS <- 500000L
EXPECTED_HDR <- c("CHR","SNP","BP","A1","A2","Z","N")

targets  <- fread(file.path(MISC, "targets.csv"))
excluded <- fread(file.path(MISC, "targets_excluded.csv"))
whm3     <- fread(WHM3, select = "SNP")
whm3_set <- whm3$SNP

# --- Validation pass ---------------------------------------------------------

validate_one <- function(code) {
  staged <- file.path(STAGING, paste0(code, ".sumstats.gz"))
  if (!file.exists(staged) || file.size(staged) == 0) {
    return(list(status = "munge_failed", n_snps_out = 0L,
                note = "staging file missing / empty"))
  }
  d <- tryCatch(fread(staged), error = function(e)
    list(err = conditionMessage(e)))
  if (is.list(d) && !is.data.table(d)) {
    return(list(status = "munge_failed", n_snps_out = 0L,
                note = paste0("fread: ", d$err)))
  }
  if (!identical(names(d), EXPECTED_HDR)) {
    return(list(status = "munge_failed", n_snps_out = nrow(d),
                note = paste0("bad header: ", paste(names(d), collapse = ","))))
  }
  n <- nrow(d)
  if (n < MIN_SNPS) {
    return(list(status = "low_snp_count", n_snps_out = n,
                note = paste0("only ", n, " variants (< ", MIN_SNPS, ")")))
  }
  n_off <- sum(!(d$SNP %in% whm3_set))
  if (n_off > 0L) {
    return(list(status = "not_hm3", n_snps_out = n,
                note = paste0(n_off, " SNPs outside w_hm3")))
  }
  list(status = "shipped", n_snps_out = n, note = "")
}

check <- lapply(targets$code, validate_one)
targets[, status     := vapply(check, `[[`, "", "status")]
targets[, n_snps_out := vapply(check, `[[`, 0L, "n_snps_out")]
targets[, note       := vapply(check, `[[`, "", "note")]

# --- Build rebuild_log.csv ---------------------------------------------------

log_kept <- targets[, .(code, source_path, status, n_snps_out, note)]
log_excl <- if (nrow(excluded)) excluded[, .(
  code, source_path, status,
  n_snps_out = 0L,
  note = fifelse(!is.na(ancestry) & ancestry != "",
                 paste0("ancestry='", ancestry, "'"), "")
)] else data.table(code=character(), source_path=character(),
                   status=character(), n_snps_out=integer(), note=character())
rebuild_log <- rbindlist(list(log_kept, log_excl), fill = TRUE)
setorder(rebuild_log, status, code)
fwrite(rebuild_log, file.path(MISC, "rebuild_log.csv"))

# --- Publish -----------------------------------------------------------------

shipped <- targets[status == "shipped", code]
cat(sprintf("validation: %d shipped, %d failed (of %d munged) + %d excluded pre-munge\n",
            length(shipped), nrow(targets) - length(shipped),
            nrow(targets), nrow(excluded)))
cat(sprintf("status breakdown: %s\n",
            paste(paste0(names(table(rebuild_log$status)), "=",
                         as.integer(table(rebuild_log$status))),
                  collapse = ", ")))

if (length(shipped) == 0L) {
  stop("No codes passed validation — refusing to touch the live panel.")
}

# --- Backup old CSV, then replace files --------------------------------------

old_csv <- file.path(SGWAS, "sumstat_repo.csv")
bak_csv <- file.path(SGWAS, paste0("sumstat_repo.csv.bak-preRebuild-",
                                   format(Sys.Date(), "%Y%m%d")))
if (file.exists(old_csv) && !file.exists(bak_csv)) {
  file.copy(old_csv, bak_csv)
  cat("backed up old CSV -> ", basename(bak_csv), "\n", sep="")
}

# Copy passing staging files. Failing codes that used to ship get their
# stale .sumstats.gz REMOVED so downstream jobs can't pick up broken data.
existing <- sub("\\.sumstats\\.gz$", "",
                list.files(SGWAS, pattern = "\\.sumstats\\.gz$"))
to_remove <- setdiff(existing, shipped)
for (c in to_remove) {
  fn <- file.path(SGWAS, paste0(c, ".sumstats.gz"))
  if (file.exists(fn)) {
    file.remove(fn)
    cat("removed stale: ", basename(fn), "\n", sep = "")
  }
}
for (c in shipped) {
  src <- file.path(STAGING, paste0(c, ".sumstats.gz"))
  dst <- file.path(SGWAS,   paste0(c, ".sumstats.gz"))
  file.copy(src, dst, overwrite = TRUE)
}
cat("published ", length(shipped), " files to ", SGWAS, "\n", sep = "")

# --- Rebuild sumstat_repo.csv (slim schema) ----------------------------------

auto_short <- function(label) {
  # Strip trailing "(...)", trailing ": <model>" clauses, and squash spaces.
  # Return input unchanged if nothing to strip.
  s <- label
  s <- gsub("\\s*\\([^)]*\\)\\s*", " ", s)
  s <- sub("\\s*:\\s*[^:]+$", "", s)
  s <- gsub("\\s+", " ", s)
  trimws(s)
}

titlecase <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  vapply(x, function(v) {
    if (is.na(v) || v == "") return(v)
    words <- strsplit(v, "\\s+")[[1]]
    paste(paste0(toupper(substr(words, 1, 1)),
                 tolower(substr(words, 2, nchar(words)))), collapse = " ")
  }, character(1))
}

ship_meta <- targets[status == "shipped"]
ship_meta[, trait_category := titlecase(trimws(as.character(category)))]

# short_label: keep curated (non-empty, <=40 chars); else auto-shorten trait_label
sl_raw <- ship_meta$short_label
tl     <- ship_meta$trait_label
short_final <- ifelse(!is.na(sl_raw) & nzchar(trimws(sl_raw)) & nchar(trimws(sl_raw)) <= 40,
                      trimws(sl_raw),
                      auto_short(tl))
ship_meta[, short_label := short_final]

slim <- ship_meta[, .(
  code,
  trait_category,
  trait_label,
  short_label,
  ancestry,
  phenotype_type,
  n_cases,
  n_controls,
  sample_size_discovery,
  pmid
)]
setorder(slim, trait_category, short_label, code)
fwrite(slim, old_csv)
cat("rebuilt ", basename(old_csv), " with ", nrow(slim), " rows / ",
    ncol(slim), " cols\n", sep = "")
