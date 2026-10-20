#!/usr/bin/Rscript
# Run GCTB SBayesRC for a single GWAS and extract SNP-heritability + polygenicity.
#
# Two-step official gctb recipe (gctb 2.5.5 binary):
#   1. gctb --ldm-eigen <ldref> --gwas-summary <gwas>.ma --impute-summary --out <prefix>
#   2. gctb --ldm-eigen <ldref> --gwas-summary <prefix>.imputed.ma --sbayes RC --annot <annot>
#        --out <prefix> --seed S --chain-length C --burn-in B --out-freq 1
# SNP-h2 is summarised from the post-burn-in MCMC trace of hsq (posterior mode +
# 95% HPD + SD) parsed from gctb's stdout; polygenicity (non-zero-effect variant
# count) comes from .parRes. Only model-level parameters are kept; the per-SNP
# .snpRes is discarded.

suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]"),
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]"),
  make_option("--gctb", action="store", default=NA, type='character',
              help="Path to the gctb binary [required]"),
  make_option("--ldref", action="store", default=NA, type='character',
              help="Per-population eigen-decomposed LD reference directory [required]"),
  make_option("--annot", action="store", default=NA, type='character',
              help="SBayesRC functional annotation file [required]"),
  make_option("--threads", action="store", default=1, type='integer',
              help="Number of threads [default: 1]"),
  make_option("--seed", action="store", default=1, type='integer',
              help="MCMC random seed [default: 1]"),
  make_option("--chain_length", action="store", default=3000, type='integer',
              help="MCMC chain length [default: 3000]"),
  make_option("--burn_in", action="store", default=1000, type='integer',
              help="MCMC burn-in iterations to discard [default: 1000]")
)

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

# Read config / identify outdir
config <- readLines(opt$config_file)
outdir <- gsub('outdir: ', '', config[grepl('outdir: ', config)])

gwas     <- opt$gwas
sbrc_dir <- paste0(outdir, '/results/', gwas, '/sbayesrc')
tmp_dir  <- paste0(sbrc_dir, '/tmp')
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# --- Build COJO/.ma input (SNP A1 A2 freq b se p N) from the cleaned sumstats ---
ss <- fread(paste0(outdir, '/results/', gwas, '/gwas_sumstat/', gwas, '.cleaned.gz'))
if (!any(names(ss) == 'BETA')) ss$BETA <- log(ss$OR)
ma <- ss[, .(SNP, A1, A2, freq = FREQ, b = BETA, se = SE, p = P, N)]
ma_file <- paste0(tmp_dir, '/', gwas, '.ma')
fwrite(ma, ma_file, sep = ' ', na = 'NA', quote = FALSE)
cat('Wrote', nrow(ma), 'variants to', ma_file, '\n')

prefix  <- paste0(tmp_dir, '/', gwas)
threads <- opt$threads

# gctb parallelises the heavy eigen (LD) operations only when given
# --multi-thread-eigen (added to both gctb calls below); without it the threads
# are largely idle. Also pin the OpenMP / BLAS thread counts so every parallel
# region uses `threads` (gctb is statically linked; OMP_NUM_THREADS still governs it).
Sys.setenv(OMP_NUM_THREADS = threads,
           OPENBLAS_NUM_THREADS = threads,
           MKL_NUM_THREADS = threads)

run <- function(cmd) {
  cat('RUN:', cmd, '\n')
  st <- system(cmd)
  if (st != 0) stop('Command failed (exit ', st, '): ', cmd)
}

# --- Step 1: QC + impute summary stats to the LD reference SNP set ---
run(paste(shQuote(opt$gctb),
          '--ldm-eigen', shQuote(opt$ldref),
          '--gwas-summary', shQuote(ma_file),
          '--impute-summary',
          '--out', shQuote(prefix),
          '--multi-thread-eigen',
          '--thread', threads))

imputed_ma <- paste0(prefix, '.imputed.ma')
if (!file.exists(imputed_ma)) stop('Expected imputed summary file not found: ', imputed_ma)

# --- Step 2: SBayesRC model ---
# --out-freq 1 prints the per-iteration MCMC trace (incl. hsq) to stdout for
# every iteration, which we capture to parse the SNP-h2 posterior. (This avoids
# gctb's --write-mcmc-txt, which also dumps multi-GB per-SNP sample files.)
sbrc_log <- paste0(prefix, '.sbrc.stdout.log')
run(paste(shQuote(opt$gctb),
          '--ldm-eigen', shQuote(opt$ldref),
          '--gwas-summary', shQuote(imputed_ma),
          '--sbayes RC',
          '--annot', shQuote(opt$annot),
          '--out', shQuote(prefix),
          '--seed', opt$seed,
          '--chain-length', opt$chain_length,
          '--burn-in', opt$burn_in,
          '--out-freq 1',
          '--multi-thread-eigen',
          '--thread', threads,
          '>', shQuote(sbrc_log), '2>&1'))

# --- Polygenicity from the parameter results (.parRes) ---
# Header `Parameter Mean SD`; rows include NumSnp1..5 (variants per mixture
# component; component 1 = zero effect), Vg1..5, SigmaSq, hsq, ResVar.
# polygenicity = non-zero-effect variants = sum(NumSnp2..NumSnp5) (posterior mean);
# pi = that / sum(NumSnp1..5). (No Pi/NnzSnp row exists.)
par_file <- paste0(prefix, '.parRes')
if (!file.exists(par_file)) stop('.parRes not found: ', par_file)
par_lines <- readLines(par_file)
get_par <- function(name) {
  hit <- par_lines[grepl(paste0('^\\s*', name, '\\b'), par_lines)]
  if (length(hit) == 0) return(c(est = NA_real_, sd = NA_real_))
  nums <- suppressWarnings(as.numeric(strsplit(trimws(hit[1]), '\\s+')[[1]][-1]))
  nums <- nums[!is.na(nums)]
  c(est = if (length(nums) >= 1) nums[1] else NA_real_,
    sd  = if (length(nums) >= 2) nums[2] else NA_real_)
}
numsnp    <- vapply(1:5, function(i) unname(get_par(paste0('NumSnp', i))['est']), numeric(1))
total_snp <- sum(numsnp, na.rm = TRUE)
n_nonzero <- sum(numsnp[-1], na.rm = TRUE)
pi_est    <- if (total_snp > 0) n_nonzero / total_snp else NA_real_

# --- SNP-h2 posterior from the MCMC trace (post-burn-in) ---
# The posterior of SNP-h2 is skewed when the true value is near zero, so we
# report the posterior MODE as the point estimate and a 95% HPD credible
# interval (per J. Zeng). The posterior SD (= SD of the MCMC samples; no
# sqrt(N) division) is the SE. Burn-in iterations (the known --burn-in value)
# are discarded.
loglines <- readLines(sbrc_log)
burn_in <- opt$burn_in

# Locate the Iter/.../hsq progress table and pull iter + hsq from each MCMC row
# (16 whitespace fields, integer iter, trailing h:m:s TimeLeft).
hdr_idx <- grep('\\<Iter\\>.*\\<hsq\\>', loglines)
h2_samples <- numeric(0)
if (length(hdr_idx) >= 1) {
  hdr <- strsplit(trimws(loglines[hdr_idx[1]]), '\\s+')[[1]]
  nf <- length(hdr); hsq_col <- which(hdr == 'hsq')
  for (ln in loglines[(hdr_idx[1] + 1):length(loglines)]) {
    f <- strsplit(trimws(ln), '\\s+')[[1]]
    if (length(f) != nf) next
    if (!grepl('^[0-9]+$', f[1])) next
    if (!grepl('^[0-9]+:[0-9]+:[0-9]+$', f[nf])) next
    it <- as.integer(f[1]); hv <- suppressWarnings(as.numeric(f[hsq_col]))
    if (!is.na(hv) && it > burn_in) h2_samples <- c(h2_samples, hv)
  }
}

# Fail hard rather than silently falling back — a missing trace means something
# is wrong with the run (flags, gctb version, truncated log) and the h2 summary
# would otherwise be unreliable.
if (length(h2_samples) < 10) {
  stop('Could not parse the post-burn-in hsq MCMC trace from ', sbrc_log,
       ' (found ', length(h2_samples), ' usable iterations; check --out-freq / --burn-in',
       ' and that the Iter...hsq table is present).')
}

hpd_interval <- function(x, prob = 0.95) {
  x <- sort(x); n <- length(x); k <- floor(prob * n)
  if (n < 2 || k < 1) return(c(NA_real_, NA_real_))
  w <- x[(k + 1):n] - x[1:(n - k)]; i <- which.min(w)
  c(x[i], x[i + k])
}

# Effective sample size: n / (1 + 2*sum of the initial positive autocorrelations).
# Computed manually (coda is not in the pipeline env). Approximate but adequate
# for flagging poor mixing.
eff_sample_size <- function(x) {
  n <- length(x); if (n < 10 || var(x) == 0) return(NA_real_)
  ac <- as.numeric(acf(x, lag.max = min(n - 1L, 1000L), plot = FALSE, demean = TRUE)$acf)[-1]
  rho_sum <- 0
  for (r in ac) { if (is.na(r) || r <= 0) break; rho_sum <- rho_sum + r }
  max(1, min(n, n / (1 + 2 * rho_sum)))
}

dens    <- density(h2_samples, from = 0)          # h2 >= 0: respect the boundary
h2_mode <- dens$x[which.max(dens$y)]
h2_mean <- mean(h2_samples); h2_sd <- sd(h2_samples)
hpd     <- hpd_interval(h2_samples, 0.95)
n_ess   <- eff_sample_size(h2_samples)

out <- data.frame(
  gwas        = gwas,
  h2_mode     = h2_mode,
  h2_mean     = h2_mean,
  h2_sd       = h2_sd,
  h2_hpd_low  = hpd[1],
  h2_hpd_high = hpd[2],
  n_post_samples = length(h2_samples),
  n_ess       = n_ess,
  n_nonzero   = n_nonzero,
  n_snp_total = total_snp,
  pi_est      = pi_est,
  stringsAsFactors = FALSE
)
cat('Parsed SBayesRC parameters:\n'); print(out)

clean_csv <- paste0(sbrc_dir, '/', gwas, '.GW.sbayesrc.clean.csv')
fwrite(out, clean_csv)
# Persist the post-burn-in hsq samples for the viewer's posterior-density plot.
samples_file <- paste0(sbrc_dir, '/', gwas, '.GW.sbayesrc.h2_samples.txt')
writeLines(format(h2_samples, scientific = FALSE, trim = TRUE), samples_file)
cat('Wrote', clean_csv, 'and', samples_file, '\n')

# --- Tidy large intermediates; keep model-parameter outputs + stdout log for provenance ---
# .parRes (summary), .enrich (per-annotation h2 enrichment) and .parSetRes (per-annotation x
# per-component parameters) are retained in the results dir for provenance. They are NOT packaged
# into the bundle or shown in the viewer (heritability + polygenicity only).
for (ext in c('.parRes', '.enrich', '.parSetRes')) {
  src <- paste0(prefix, ext)
  if (file.exists(src)) file.copy(src, paste0(sbrc_dir, '/', gwas, ext), overwrite = TRUE)
}
if (file.exists(sbrc_log)) file.copy(sbrc_log, paste0(sbrc_dir, '/', gwas, '.sbrc.stdout.log'), overwrite = TRUE)
for (f in c(ma_file, imputed_ma, paste0(prefix, '.snpRes'))) {
  if (file.exists(f)) file.remove(f)
}
