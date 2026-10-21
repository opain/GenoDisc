#!/usr/bin/Rscript
# Run LDAK SumHer (BLD-LDAK model) for a single GWAS and extract SNP-heritability.
# This is a THIRD, independent SNP-h2 estimator alongside LDSC (scripts) and the
# two SBayesRC routes (scripts/sbayesrc.R, scripts/sbayesrc_r.R). It drives the
# LDAK binary's --sum-hers over a pre-computed BLD-LDAK HapMap3 tagging file.
#
# Observed-scale only (no --prevalence / --ascertainment: GenoDisc N may be total
# N or Neff depending on the upload, so ascertainment is ambiguous). No --intercept
# / --genomic-control (per LDAK's SumHer recommendation). Large-effect predictors
# are down-weighted via --cutoff. SumHer assumes a different heritability model
# from LDSC/SBayesRC, so its estimate is not expected to match theirs.

suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]"),
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]"),
  make_option("--ldak", action="store", default=NA, type='character',
              help="Path to the LDAK executable [required]"),
  make_option("--tagfile", action="store", default=NA, type='character',
              help="Path to the BLD-LDAK tagging file for this population [required]"),
  make_option("--cutoff", action="store", default=0.01, type='double',
              help="LDAK --cutoff (exclude large-effect predictors) [default: 0.01]"),
  make_option("--max_missing", action="store", default=0.2, type='double',
              help="Warn if the fraction of tagging predictors missing exceeds this [default: 0.2]"),
  make_option("--threads", action="store", default=1, type='integer',
              help="LDAK --max-threads [default: 1]")
)

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

# Read config / identify outdir
config <- readLines(opt$config_file)
outdir <- gsub('outdir: ', '', config[grepl('outdir: ', config)])

gwas       <- opt$gwas
sumher_dir <- paste0(outdir, '/results/', gwas, '/sumher')
tmp_dir    <- paste0(sumher_dir, '/tmp')
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
# Start from a clean tmp dir so a re-run can't reuse a previous run's LDAK outputs.
unlink(list.files(tmp_dir, full.names = TRUE), recursive = TRUE, force = TRUE)

# --- Build the LDAK summary file: Predictor A1 A2 n Z A1Freq -----------------
# Predictor = rsID (the pre-computed HapMap3 taggings are rsID-keyed), n = per-SNP
# sample size, Z = BETA/SE, A1Freq = effect-allele frequency. LDAK matches by the
# Predictor name and reads columns by header, so this exact header is what it wants.
ss <- fread(paste0(outdir, '/results/', gwas, '/gwas_sumstat/', gwas, '.cleaned.gz'))
if (!any(names(ss) == 'BETA')) ss$BETA <- log(ss$OR)

n_raw <- nrow(ss)
# Defensive QC: single-character alleles, finite stats, positive SE/N, and no
# duplicate predictors (LDAK errors on duplicates).
ss <- ss[nchar(A1) == 1 & nchar(A2) == 1 &
         is.finite(BETA) & is.finite(SE) & SE > 0 &
         is.finite(N) & N > 0 & is.finite(FREQ)]
ss <- ss[!duplicated(SNP)]
n_kept <- nrow(ss)
cat('Summary file: kept', n_kept, 'of', n_raw, 'variants after QC\n')
if (n_kept == 0) stop('No variants left after QC for ', gwas)

summ <- ss[, .(Predictor = SNP, A1 = A1, A2 = A2, n = N, Z = BETA / SE, A1Freq = FREQ)]
summ_predictors <- summ$Predictor
summ_file <- paste0(tmp_dir, '/', gwas, '.sumher.summary')
fwrite(summ, summ_file, sep = ' ', na = 'NA', quote = FALSE)
rm(ss, summ); gc()

# --- Record the LDAK version line for provenance -----------------------------
# LDAK's no-arg banner prints the version on its own line ("Version 6 - Help
# pages at www.dougspeed.com"); grab that line specifically.
ldak_banner <- tryCatch({
  bl <- suppressWarnings(system2(opt$ldak, args = character(0), stdout = TRUE, stderr = TRUE))
  vl <- bl[grepl('Version', bl)]
  if (length(vl)) trimws(vl[1]) else NA_character_
}, error = function(e) NA_character_)
cat('LDAK version:', ldak_banner, '\n')

# LDAK writes -9999 as an "undefined" sentinel (e.g. an SE when the fit did not
# converge); treat it as NA.
na_if_sentinel <- function(x) ifelse(is.finite(x) & x == -9999, NA_real_, x)

# --- Run SumHer --------------------------------------------------------------
out_prefix <- paste0(tmp_dir, '/sumher')
ldak_args <- c('--sum-hers', out_prefix,
               '--summary', summ_file,
               '--tagfile', opt$tagfile,
               '--check-sums', 'NO',
               '--cutoff', format(opt$cutoff, scientific = FALSE),
               '--max-threads', opt$threads)
cat('Running:', opt$ldak, paste(ldak_args, collapse = ' '), '\n')
ldak_out <- suppressWarnings(system2(opt$ldak, args = ldak_args, stdout = TRUE, stderr = TRUE))
cat(ldak_out, sep = '\n'); cat('\n')
ldak_status <- attr(ldak_out, 'status')
if (!is.null(ldak_status) && ldak_status != 0)
  stop('LDAK --sum-hers failed (exit ', ldak_status, ') for ', gwas)

hers_file  <- paste0(out_prefix, '.hers')
extra_file <- paste0(out_prefix, '.extra')
enrich_file<- paste0(out_prefix, '.enrich')
if (!file.exists(hers_file)) stop('LDAK did not write ', hers_file)

# --- Parse SNP-h2 + SD from the final (all-SNPs) row of .hers ----------------
# .hers has one row per annotation category then a genome-wide total row
# ("Her_All" / "...All..."); the total row's heritability + SD is the overall
# SNP-h2. Parse robustly by header name (confirmed against real output at test).
hers <- fread(hers_file)
comp_col <- names(hers)[1]
total_row <- hers[grepl('All', get(comp_col), ignore.case = TRUE)]
if (nrow(total_row) == 0) total_row <- hers[nrow(hers)]
h2_col <- names(hers)[grepl('^herit', names(hers), ignore.case = TRUE)]
sd_col <- names(hers)[grepl('^sd$|^std|^se$', names(hers), ignore.case = TRUE)]
if (length(h2_col) == 0) h2_col <- names(hers)[2]
if (length(sd_col) == 0) sd_col <- names(hers)[min(3, ncol(hers))]
h2    <- na_if_sentinel(as.numeric(total_row[[h2_col[1]]][1]))
h2_sd <- na_if_sentinel(as.numeric(total_row[[sd_col[1]]][1]))
if (length(total_row) == 0 || is.null(total_row[[h2_col[1]]])) stop('Could not parse SNP-h2 from ', hers_file)

# --- .extra: null/alt log-likelihoods + convergence --------------------------
# .extra is "Key Value" per line, e.g. "Null_logl -1303.28", "Alt_logl -1283.35",
# "Converged NO". Parse by key.
extra_val <- function(lines, key) {
  hit <- lines[grepl(paste0('^', key, '\\b'), lines)]
  if (length(hit) == 0) return(NA_character_)
  trimws(sub(paste0('^', key), '', hit[1]))
}
extra_lines <- if (file.exists(extra_file)) readLines(extra_file) else character(0)
logl_alt  <- suppressWarnings(as.numeric(extra_val(extra_lines, 'Alt_logl')))
logl_null <- suppressWarnings(as.numeric(extra_val(extra_lines, 'Null_logl')))
converged <- toupper(extra_val(extra_lines, 'Converged')) == 'YES'
num_datapoints <- suppressWarnings(as.integer(extra_val(extra_lines, 'Num_Datapoints')))

# Predictors excluded by --cutoff, from LDAK's screen output. LDAK prints either
# "There are no predictors that explain at least X ..." (0 excluded) or
# "Excluding N predictors that explain ...".
if (any(grepl('no predictors that explain at least', ldak_out, ignore.case = TRUE))) {
  n_excluded_cutoff <- 0L
} else {
  excl <- ldak_out[grepl('[Ee]xcluding [0-9]+ predictors', ldak_out)]
  n_excluded_cutoff <- if (length(excl)) as.integer(sub('.*[Ee]xcluding ([0-9]+) predictors.*', '\\1', excl[1])) else NA_integer_
}

# --- Missingness: tagging predictors absent from the summary file ------------
# We compute n_matched ourselves (summary predictors present in the tagging set)
# and take n_tagging from LDAK's authoritative "Reading details for N predictors"
# line (the raw tagging-file row count includes ~66 trailing annotation-category
# rows). prop_missing = 1 - n_matched/n_tagging. Both are cross-checked against
# LDAK's own counts (Num_Datapoints in .extra; the "only N of the M predictors"
# warning) and any mismatch is logged.
tag_pred <- tryCatch(fread(opt$tagfile, select = 1L)[[1]], error = function(e) character(0))
n_matched <- sum(summ_predictors %in% tag_pred)
read_line <- ldak_out[grepl('Reading details for [0-9]+ predictors', ldak_out)]
n_tagging <- if (length(read_line)) as.integer(sub('.*for ([0-9]+) predictors.*', '\\1', read_line[1])) else length(tag_pred)
prop_missing <- if (!is.na(n_tagging) && n_tagging > 0) (n_tagging - n_matched) / n_tagging else NA_real_
missing_warning <- !is.na(prop_missing) && prop_missing > opt$max_missing
cat(sprintf('Missingness: n_tagging=%s n_matched=%s prop_missing=%s warning=%s converged=%s\n',
            n_tagging, n_matched, signif(prop_missing, 4), missing_warning, converged))
# Cross-check against LDAK's own counts.
if (!is.na(num_datapoints) && num_datapoints != n_matched)
  cat(sprintf('NOTE: LDAK Num_Datapoints=%s differs from our n_matched=%s\n', num_datapoints, n_matched))
only_line <- ldak_out[grepl('only [0-9]+ of the [0-9]+ predictors', ldak_out)]
if (length(only_line)) cat('LDAK coverage: ', trimws(only_line[1]), '\n')

# --- Write the clean.csv (flat val + provenance) -----------------------------
out <- data.frame(
  gwas              = gwas,
  h2                = h2,
  h2_sd             = h2_sd,
  model             = 'BLD-LDAK',
  n_tagging         = n_tagging,
  n_matched         = n_matched,
  prop_missing      = prop_missing,
  missing_warning   = missing_warning,
  converged         = converged,
  n_excluded_cutoff = n_excluded_cutoff,
  logl_null         = logl_null,
  logl_alt          = logl_alt,
  cutoff            = opt$cutoff,
  ldak_banner       = ldak_banner,
  stringsAsFactors  = FALSE
)
cat('Parsed SumHer result:\n'); print(out)

clean_csv <- paste0(sumher_dir, '/', gwas, '.sumher.clean.csv')
fwrite(out, clean_csv)

# --- Keep .hers / .extra / .enrich + a provenance log; drop the summary file ---
for (f in c(hers_file, extra_file, enrich_file,
            paste0(out_prefix, '.cats'), paste0(out_prefix, '.share'))) {
  if (file.exists(f)) file.copy(f, paste0(sumher_dir, '/', gwas, '.', sub('^sumher\\.', '', basename(f))), overwrite = TRUE)
}
writeLines(c(paste('LDAK:', opt$ldak), paste('banner:', ldak_banner),
             paste('command:', opt$ldak, paste(ldak_args, collapse = ' ')),
             '', ldak_out),
           paste0(sumher_dir, '/', gwas, '.ldak.log'))
if (file.exists(summ_file)) file.remove(summ_file)
cat('Wrote', clean_csv, '\n')
