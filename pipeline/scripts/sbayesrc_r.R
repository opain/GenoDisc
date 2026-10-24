#!/usr/bin/Rscript
# Run SBayesRC via the SBayesRC R package (the GenoPred route) for a single GWAS
# and extract SNP-heritability + polygenicity. This is the ALTERNATIVE SBayesRC
# implementation: it uses the SBayesRC R package (tidy -> impute -> sbayesrc),
# which is multi-threaded (OpenMP) over a precomputed eigen-LD panel and is
# typically faster than driving the gctb binary from scratch (scripts/sbayesrc.R).
#
# Unlike GenoPred (which keeps only the per-SNP PGS weights and discards the
# parameter results), we RETAIN the <prefix>.parRes written by SBayesRC::sbayesrc
# and report the model-level SNP-h2 (posterior mean + SD) and polygenicity
# (non-zero-effect variant count + pi). No PGS weights are kept.

suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]"),
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]"),
  make_option("--ldref", action="store", default=NA, type='character',
              help="Per-population SBayesRC eigen-LD reference directory [required]"),
  make_option("--annot", action="store", default=NA, type='character',
              help="SBayesRC functional annotation file [required]"),
  make_option("--threads", action="store", default=1, type='integer',
              help="Number of threads (OpenMP) [default: 1]"),
  make_option("--seed", action="store", default=1, type='integer',
              help="Random seed [default: 1]")
)

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

# Read config / identify outdir
config <- readLines(opt$config_file)
outdir <- gsub('outdir: ', '', config[grepl('outdir: ', config)])

gwas     <- opt$gwas
sbrc_dir <- paste0(outdir, '/results/', gwas, '/sbayesrc_r')
tmp_dir  <- paste0(sbrc_dir, '/tmp')
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
# Start from a clean tmp dir: the SBayesRC package reuses existing intermediates
# (tidy/impute/sbrc.*) if present, so a re-run would otherwise silently reuse a
# previous run's results instead of recomputing.
unlink(list.files(tmp_dir, full.names = TRUE), recursive = TRUE, force = TRUE)

# The SBayesRC package parallelises its (RcppEigen) linear algebra via OpenMP;
# pin the thread count so it uses `threads` cores.
threads <- opt$threads
Sys.setenv(OMP_NUM_THREADS = threads,
           OPENBLAS_NUM_THREADS = threads,
           MKL_NUM_THREADS = threads)

set.seed(opt$seed)

# --- Build COJO/.ma input (SNP A1 A2 freq b se p N) from the cleaned sumstats ---
# (Same input the gctb-binary method builds; the R package's tidy() step then QCs
#  it against the LD reference.)
ss <- fread(paste0(outdir, '/results/', gwas, '/gwas_sumstat/', gwas, '.cleaned.gz'))
if (!any(names(ss) == 'BETA')) ss$BETA <- log(ss$OR)
# Allele frequency is 'FREQ' when the GWAS supplied its own, else 'REF.FREQ'
# (reference-panel frequency the cleaner backfills). Accept either, mirroring
# process_cojo.R, so a GWAS without a frequency column still works.
freq_col <- if ('FREQ' %in% names(ss)) 'FREQ' else 'REF.FREQ'
if (!freq_col %in% names(ss)) stop("No FREQ or REF.FREQ column in cleaned sumstats for ", gwas)
ma <- ss[, .(SNP, A1, A2, freq = get(freq_col), b = BETA, se = SE, p = P, N)]
ma_file <- paste0(tmp_dir, '/', gwas, '.ma')
fwrite(ma, ma_file, sep = ' ', na = 'NA', quote = FALSE)
cat('Wrote', nrow(ma), 'variants to', ma_file, '\n')
rm(ss, ma); gc()

suppressMessages(library(SBayesRC))

tidy_ma <- paste0(tmp_dir, '/tidy.ma')
imp_ma  <- paste0(tmp_dir, '/tidy.imp.ma')
prefix  <- paste0(tmp_dir, '/sbrc')

# --- Step 1: tidy (QC the sumstats to the LD reference SNP set) ---
cat('SBayesRC::tidy ...\n')
SBayesRC::tidy(mafile = ma_file, LDdir = opt$ldref, output = tidy_ma, log2file = TRUE)

# --- Step 2: impute summary stats to the LD reference ---
cat('SBayesRC::impute ...\n')
SBayesRC::impute(mafile = tidy_ma, LDdir = opt$ldref, output = imp_ma, log2file = TRUE)

# --- Step 3: SBayesRC model (with GenoPred's tuning fallbacks) ---
cat('SBayesRC::sbayesrc ...\n')
run_sbrc <- function(...) SBayesRC::sbayesrc(mafile = imp_ma, LDdir = opt$ldref,
                                             outPrefix = prefix, annot = opt$annot,
                                             log2file = FALSE, ...)
tryCatch(
  run_sbrc(),
  error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("All correlations are negative", msg)) {
      cat('Retry: all correlations negative -> bTune = FALSE\n'); run_sbrc(bTune = FALSE)
    } else if (grepl("the best parameter is the minimumn threshold", msg)) {
      cat('Retry: best param at min threshold -> tuneStep set\n')
      run_sbrc(tuneStep = c(0.995, 0.9, 0.8, 0.7, 0.6))
    } else if (grepl("Invalid tune outputs, all correlations are invalid number", msg)) {
      cat('Retry: invalid tune outputs -> bTune = FALSE\n'); run_sbrc(bTune = FALSE)
    } else {
      stop("An unexpected error occurred in SBayesRC::sbayesrc: ", msg)
    }
  }
)

# --- SNP-h2 + polygenicity from the SBayesRC parameter file (.par) ---
# The SBayesRC R package writes `<prefix>.par` (tab-separated: Item / Mean / SD),
# with rows hsq, nnz (non-zero-effect variant count), SigmaSq, NumSnp2..5, etc.
# (Different filename/layout from the gctb-binary method's .parRes.)
par_file <- paste0(prefix, '.par')
if (!file.exists(par_file)) stop('SBayesRC .par not found (did not complete): ', par_file)
par_lines <- readLines(par_file)
get_par <- function(name) {
  hit <- par_lines[grepl(paste0('^\\s*', name, '\\b'), par_lines)]
  if (length(hit) == 0) return(c(est = NA_real_, sd = NA_real_))
  nums <- suppressWarnings(as.numeric(strsplit(trimws(hit[1]), '\\s+')[[1]][-1]))
  nums <- nums[!is.na(nums)]
  c(est = if (length(nums) >= 1) nums[1] else NA_real_,
    sd  = if (length(nums) >= 2) nums[2] else NA_real_)
}
hsq       <- get_par('hsq')
n_nonzero <- unname(get_par('nnz')['est'])
# fall back to summing the non-null mixture components if `nnz` is absent
if (is.na(n_nonzero)) {
  n_nonzero <- sum(vapply(2:5, function(i) unname(get_par(paste0('NumSnp', i))['est']), numeric(1)), na.rm = TRUE)
}
# The .par has no total-SNP count; use the model SNP set (rows of the imputed .ma).
n_snp_total <- tryCatch(nrow(data.table::fread(imp_ma, select = 1L)), error = function(e) NA_real_)
pi_est      <- if (!is.na(n_snp_total) && n_snp_total > 0) n_nonzero / n_snp_total else NA_real_

if (is.na(hsq['est'])) stop('Could not parse hsq from ', par_file)

out <- data.frame(
  gwas        = gwas,
  h2_mean     = unname(hsq['est']),
  h2_sd       = unname(hsq['sd']),
  n_nonzero   = n_nonzero,
  n_snp_total = n_snp_total,
  pi_est      = pi_est,
  stringsAsFactors = FALSE
)
cat('Parsed SBayesRC (R package) parameters:\n'); print(out)

clean_csv <- paste0(sbrc_dir, '/', gwas, '.GW.sbayesrc_r.clean.csv')
fwrite(out, clean_csv)

# --- Keep the parameter file + annotation enrichment + logs for provenance;
#     drop the large intermediates (imputed .ma, per-SNP effects, rds). ---
if (file.exists(par_file)) file.copy(par_file, paste0(sbrc_dir, '/', gwas, '.par'), overwrite = TRUE)
enrich <- paste0(prefix, '.AnnoPerSnpHsqEnrichment')
if (file.exists(enrich)) file.copy(enrich, paste0(sbrc_dir, '/', gwas, '.AnnoPerSnpHsqEnrichment'), overwrite = TRUE)
for (lg in c(paste0(tidy_ma, '.log'), paste0(imp_ma, '.log'))) {
  if (file.exists(lg)) file.copy(lg, paste0(sbrc_dir, '/', gwas, '.', basename(lg)), overwrite = TRUE)
}
for (f in c(ma_file, tidy_ma, paste0(tidy_ma, '.full'), imp_ma,
            paste0(prefix, '.txt'), paste0(prefix, '.rds'), paste0(prefix, '.snpRes'))) {
  if (file.exists(f)) file.remove(f)
}
cat('Wrote', clean_csv, '\n')
