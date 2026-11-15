#!/usr/bin/Rscript
# Builds the per-locus data the viewer's interactive locus plot needs for every
# clump lead (P < 1e-5) - including suggestive loci that fine-mapping (which
# otherwise supplies this) does not cover:
#   * a full windowed per-SNP association table (every GWAS variant in the clump
#     window, with P and r^2-to-lead) so the association panel is dense and
#     LD-coloured, not restricted to the thinned genome-wide manhattan_data; and
#   * a recombination-rate track for the window.
# The locus plots themselves are drawn in the Shiny app from this packaged data,
# so this step no longer renders static locuszoom PNGs.
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

suppressMessages({
  library(data.table)
})
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

outdir  <- read_param(config = opt$config_file, param = 'outdir',  return_obj = F)
resdir  <- read_param(config = opt$config_file, param = 'resdir',  return_obj = F)
gwas_list <- read_param(config = opt$config_file, param = 'gwas_list', return_obj = T)
population <- gwas_list$population[gwas_list$name == opt$gwas]

# The sentinel keeps its historical name/location so the DAG wiring (report.smk,
# package_results) is unchanged; the directory now only holds the sentinel.
locus_dir <- paste0(outdir, '/results/', opt$gwas, '/locus_plots/')
dir.create(locus_dir, recursive = TRUE, showWarnings = FALSE)
sentinel  <- paste0(locus_dir, opt$gwas, '.locus_plots.done')
clump_dir <- paste0(outdir, '/results/', opt$gwas, '/clump/')

clump_csv <- paste0(clump_dir, opt$gwas, '.GW.clump.clean.csv')
idx <- fread(clump_csv)
idx <- idx[!is.na(P) & P > 0 & P < 1e-5]
setorder(idx, P)
# No cap: the viewer's locus plot can be opened for any clump lead, so every
# suggestive locus needs its own data (previously only the top 50 did, a limit
# that existed to bound the number of static PNGs this step rendered).

if (nrow(idx) == 0) {
  message("No index variants with P < 1e-5; writing sentinel and exiting.")
  file.create(sentinel)
  quit(save = "no", status = 0)
}

# Parse .clumped per-chromosome files to build index SNP -> member SNP map
clumped_files <- Sys.glob(paste0(clump_dir, opt$gwas, '_chr*.clumped'))
clump_members_by_index <- list()
for (f in clumped_files) {
  tmp <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(tmp) || !all(c("SNP", "SP2") %in% names(tmp))) next
  for (i in seq_len(nrow(tmp))) {
    sp2 <- tmp$SP2[i]
    if (is.na(sp2) || sp2 == "NONE") {
      clump_members_by_index[[tmp$SNP[i]]] <- character(0)
    } else {
      members <- unlist(strsplit(sp2, ","))
      members <- gsub("\\(.*$", "", members)
      clump_members_by_index[[tmp$SNP[i]]] <- members
    }
  }
}

ss <- fread(paste0(outdir, '/results/', opt$gwas, '/gwas_sumstat/', opt$gwas, '.cleaned.gz'),
            select = c("CHR", "BP", "SNP", "P"))
ss <- ss[!is.na(P) & P > 0 & P <= 1 & !is.na(BP) & !is.na(CHR)]
ss[, CHR := as.integer(CHR)]

sanitise_id <- function(x) gsub("[/:]", "_", x)

# Per-chr interpolated genetic map -> cM/Mb (mirrors process_finemap.R get_recomb).
recomb_cache <- new.env(parent = emptyenv())
get_recomb <- function(chr_i) {
  key <- as.character(chr_i); hit <- recomb_cache[[key]]
  if (!is.null(hit)) return(if (is.data.table(hit)) hit else NULL)
  path <- sprintf("%s/data/recomb_maps/chr%s.interpolated_genetic_map.gz", resdir, chr_i)
  m <- tryCatch(fread(path, header = FALSE, col.names = c("SNP", "BP", "cM")),
                error = function(e) NULL)
  if (is.null(m) || nrow(m) < 2) { assign(key, list(), envir = recomb_cache); return(NULL) }
  m <- m[order(BP)]; m[, rate := c(diff(cM) / diff(BP) * 1e6, NA_real_)]; m <- m[is.finite(rate)]
  assign(key, m, envir = recomb_cache); m
}

# For one clump lead: the full windowed per-SNP association (every GWAS variant
# in the clump window +-500kb, with P and r^2-to-lead) and the recombination
# track. Accumulates into locus_acc / recomb_acc.
compute_locus <- function(idx_row) {
  snp_id <- idx_row$SNP
  chr_i  <- as.integer(idx_row$CHR)
  bp_i   <- as.integer(idx_row$BP)

  members <- clump_members_by_index[[snp_id]]
  if (is.null(members)) members <- character(0)
  member_bps <- ss[CHR == chr_i & SNP %in% members, BP]
  all_bps <- c(bp_i, member_bps)
  win_lo <- min(all_bps) - 500000L
  win_hi <- max(all_bps) + 500000L

  win <- ss[CHR == chr_i & BP >= win_lo & BP <= win_hi]
  if (nrow(win) == 0) {
    warning(sprintf("Locus %s: no variants in window; skipping.", snp_id))
    return(invisible(NULL))
  }

  # r^2 of every window variant to the lead (PLINK). On failure the association
  # table is still emitted (r^2 = NA -> grey points), so the scatter is complete.
  bfile <- sprintf("%s/data/1kg_hgdp/1KG.Phase3.%s.MAF_001.chr%d", resdir, population, chr_i)
  tmp_prefix <- tempfile(pattern = paste0("ld_", sanitise_id(snp_id), "_"))
  window_kb_half <- ceiling(max(bp_i - win_lo, win_hi - bp_i) / 1000)
  cmd <- sprintf(
    "plink --bfile %s --r2 --ld-snp %s --ld-window-kb %d --ld-window 99999 --ld-window-r2 0 --out %s",
    bfile, shQuote(snp_id), window_kb_half, tmp_prefix
  )
  rc <- system(paste(cmd, "> /dev/null 2>&1"))
  ld_path <- paste0(tmp_prefix, ".ld")
  ld <- if (rc == 0 && file.exists(ld_path))
    tryCatch(fread(ld_path, select = c("SNP_B", "R2")), error = function(e) NULL) else NULL
  unlink(paste0(tmp_prefix, "*"))
  if (!is.null(ld) && nrow(ld) > 0) {
    setnames(ld, c("SNP", "R2"))
    win <- merge(win, ld, by = "SNP", all.x = TRUE, sort = FALSE)
  } else {
    warning(sprintf("Locus %s: no LD output; association kept without r^2.", snp_id))
    win[, R2 := NA_real_]
  }
  locus_acc[[snp_id]] <<- win[, .(locus = snp_id, SNP = SNP, CHR = CHR, BP = BP, P = P, r2 = R2)]

  # Recombination-rate track for the window, thinned to ~500 points.
  rm_map <- get_recomb(chr_i)
  if (!is.null(rm_map)) {
    rw <- rm_map[BP >= win_lo & BP <= win_hi, .(BP, rate)]
    if (nrow(rw) >= 2) {
      if (nrow(rw) > 500) rw <- rw[seq(1, nrow(rw), length.out = 500)]
      rw[, locus := snp_id]
      recomb_acc[[snp_id]] <<- rw[, .(locus, BP, rate)]
    }
  }
  invisible(NULL)
}

locus_acc <- list()
recomb_acc <- list()
for (i in seq_len(nrow(idx))) {
  tryCatch(compute_locus(idx[i]),
           error = function(e) warning(sprintf("Locus %s: %s", idx$SNP[i], conditionMessage(e))))
}

# Full windowed per-SNP association + r^2-to-lead for every clump lead (the file
# keeps its historical ld_to_lead name; it now also carries CHR/BP/P), plus the
# per-locus recombination track - both for the viewer's locus dashboard.
if (length(locus_acc) > 0) {
  fwrite(data.table::rbindlist(locus_acc),
         paste0(clump_dir, opt$gwas, '.GW.clump.ld_to_lead.csv'))
}
if (length(recomb_acc) > 0) {
  fwrite(data.table::rbindlist(recomb_acc),
         paste0(clump_dir, opt$gwas, '.GW.clump.recomb.csv'))
}

file.create(sentinel)
