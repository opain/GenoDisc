#!/usr/bin/Rscript
library("optparse")

option_list <- list(
  make_option("--gwas", action = "store", default = NA, type = "character",
              help = "Primary GWAS name [required]"),
  make_option("--config_file", action = "store", default = NA, type = "character",
              help = "Path to config file [required]"),
  make_option("--pipeline_dir", action = "store", default = NA, type = "character",
              help = "Path to the pipeline directory [required]")
)

opt <- parse_args(OptionParser(option_list = option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

outdir      <- read_param(config = opt$config_file, param = 'outdir',           return_obj = F)
gencor_path <- read_param(config = opt$config_file, param = 'gencor_gwas_list', return_obj = F)

if (is.null(gencor_path) || is.na(gencor_path) || gencor_path == '') {
  stop("process_ldsc_gencor.R invoked but gencor_gwas_list is not set in the config.")
}

secondary <- fread(gencor_path, sep = ' ', header = TRUE)
gencor_dir <- file.path(outdir, 'results', opt$gwas, 'gencor')

# Any column in gencor_gwas_list beyond the known system columns is treated
# as facet-able metadata (e.g. category) and preserved through to the
# results CSV so the Shiny app can offer it as a "Facet by" option.
KNOWN_GENCOR_LIST_COLS <- c("name", "label", "path", "population", "n",
                             "sampling", "prevalence", "mean", "sd")
extra_cols <- setdiff(names(secondary), KNOWN_GENCOR_LIST_COLS)

# Parse rg / SE / p / gcov_int / n_snps from a single LDSC --rg log
parse_pair_log <- function(log_file) {
  out <- list(rg = NA_real_, rg_se = NA_real_, rg_p = NA_real_,
              gcov_int = NA_real_, n_snps = NA_integer_)
  if (!file.exists(log_file)) return(out)
  lines <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) character(0))
  if (length(lines) == 0) return(out)
  # Explicit failure markers: shell-level (||) and LDSC's own per-pair error.
  # LDSC exits 0 even when it cannot read the secondary file — it catches the
  # IOError, prints a NA-filled summary row, and continues — so we have to look
  # for "ERROR computing rg" to distinguish a real result from a swallowed
  # failure. Without this check, the "After merging with regression SNP LD"
  # line earlier in the log (which refers to the primary alone) bleeds into
  # n_snps and gives a misleading row.
  if (any(grepl('^GENCOR_PAIR_FAILED', lines))) return(out)
  if (any(grepl('ERROR computing rg', lines)))  return(out)

  # Summary table emitted by ldsc.py --rg. Header line:
  # "p1                 p2                 rg     se     z     p     ..."
  hdr_idx <- grep('^\\s*p1\\s+p2\\s+rg\\s+se\\s+z\\s+p', lines)
  if (length(hdr_idx) == 1 && hdr_idx < length(lines)) {
    header <- strsplit(trimws(lines[hdr_idx]), '\\s+')[[1]]
    vals_line <- lines[hdr_idx + 1]
    vals <- strsplit(trimws(vals_line), '\\s+')[[1]]
    if (length(vals) >= length(header)) {
      get_num <- function(col) {
        i <- which(header == col)
        if (length(i) != 1) return(NA_real_)
        v <- suppressWarnings(as.numeric(vals[i]))
        if (is.na(v)) NA_real_ else v
      }
      out$rg    <- get_num('rg')
      out$rg_se <- get_num('se')
      out$rg_p  <- get_num('p')
      out$gcov_int <- get_num('gcov_int')
    }
  }

  # SNP count: LDSC logs include lines like
  # "After merging with regression SNP LD, NNN SNPs remain."
  snp_line <- lines[grepl('After merging with regression SNP LD,.*SNPs remain', lines)]
  if (length(snp_line) > 0) {
    n <- suppressWarnings(as.integer(gsub('.*,\\s*([0-9]+)\\s+SNPs.*', '\\1', snp_line[length(snp_line)])))
    if (!is.na(n)) out$n_snps <- n
  }

  out
}

rows <- vector('list', nrow(secondary))
for (i in seq_len(nrow(secondary))) {
  nm <- secondary$name[i]
  log_file <- file.path(gencor_dir, paste0(opt$gwas, '__', nm, '.log'))
  parsed <- tryCatch(parse_pair_log(log_file),
                     error = function(e) list(rg = NA_real_, rg_se = NA_real_, rg_p = NA_real_,
                                              gcov_int = NA_real_, n_snps = NA_integer_))
  rows[[i]] <- data.table(
    name     = nm,
    label    = secondary$label[i],
    rg       = parsed$rg,
    rg_se    = parsed$rg_se,
    rg_p     = parsed$rg_p,
    n_snps   = parsed$n_snps,
    gcov_int = parsed$gcov_int
  )
}
res <- rbindlist(rows)

# BH-FDR across non-NA p-values
res[, rg_p_fdr := NA_real_]
valid <- !is.na(res$rg_p)
if (any(valid)) {
  res[valid, rg_p_fdr := p.adjust(rg_p, method = 'BH')]
}

# Attach preserved extras by matching on name (defensive; res / secondary
# should already be in the same order, but cbind-by-position would silently
# misalign if that ever changed).
if (length(extra_cols) > 0) {
  idx <- match(res$name, secondary$name)
  extras <- secondary[idx, ..extra_cols]
  res <- cbind(res, extras)
}

core_cols <- c('name', 'label', 'rg', 'rg_se', 'rg_p', 'rg_p_fdr', 'n_snps', 'gcov_int')
setcolorder(res, c(core_cols, setdiff(names(res), core_cols)))

dir.create(gencor_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(gencor_dir, paste0(opt$gwas, '_gencor_res.csv'))
fwrite(res, out_csv)
