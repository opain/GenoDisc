#!/usr/bin/Rscript
# MAGMA gene-set x GTEx-tissue interaction analysis.
# Post-hoc replication of de Leeuw et al. 2018 (Nat Commun 9:3768), with
# the deliberate deviation that all tissues in GTEx_v8_tissue.tsv are
# tested (not restricted to marginally-significant ones).
#
# Per GWAS, per tissue t, per gene set S with >= min_set_size genes, fit:
#   Z ~ [MAGMA internal covariates] + Average + Average x S + t + S + t x S
# and test t x S one-sided positive. Average is the row-mean expression
# across tissues in GTEx_v8_tissue.tsv, i.e. the overall expression.
#
# Reuses magma_gene_level.genes.raw and GTEx_v8_tissue.tsv - the MAGMA
# gene analysis is never re-run here.

start.time <- Sys.time()

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(yaml))
suppressMessages(library(jsonlite))
suppressMessages(library(parallel))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--pipeline_dir", type = "character"),
  make_option("--gwas",         type = "character"),
  make_option("--config_file",  type = "character"),
  make_option("--resdir",       type = "character"),
  make_option("--outdir",       type = "character"),
  make_option("--n_cores",      type = "integer", default = 1L,
              help = "Parallel workers for the per-set MAGMA loop [default: 1]"),
  make_option("--tmpdir",       type = "character", default = NA_character_,
              help = "Root for the per-set MAGMA temp workspace. 'auto' or unset -> prefer /dev/shm if writable, else R's tempdir(). On shared-filesystem clusters (CephFS, Lustre etc.) 8+ concurrent MAGMA workers block on I/O and slow ~9x; pointing at a node-local RAM disk (/dev/shm) or SSD restores expected speed.")
)))

# The pipeline sets resdir: NA to mean "resources"; mirror the
# dependencies.smk fallback so this script works standalone too.
if (is.null(opt$resdir) || identical(opt$resdir, "NA")) opt$resdir <- "resources"

cfg <- yaml::read_yaml(opt$config_file)
mi  <- cfg$magma_interaction
if (is.null(mi) || !isTRUE(mi$enabled)) {
  stop("magma_interaction.enabled is not true in ", opt$config_file)
}
min_set_size    <- as.integer(if (is.null(mi$min_set_size))          100L  else mi$min_set_size)
fdr_thr         <- as.numeric(if (is.null(mi$fdr_threshold))         0.05  else mi$fdr_threshold)
followup_p_thr  <- as.numeric(if (is.null(mi$followup_p_threshold))  1e-4  else mi$followup_p_threshold)
bias_p_thr      <- as.numeric(if (is.null(mi$bias_flag_p_threshold)) 1e-3  else mi$bias_flag_p_threshold)
bias_min_hits   <- as.integer(if (is.null(mi$bias_flag_min_tissues)) 5L    else mi$bias_flag_min_tissues)
gmt_dir      <- cfg$pathway_gmt_dir
if (is.null(gmt_dir) || identical(gmt_dir, "NA")) {
  stop("pathway_gmt_dir must be set for magma_interaction")
}

magma_bin  <- file.path(opt$resdir, "software", "magma", "magma")
tissue_tsv <- file.path(opt$resdir, "data", "gtex", "GTEx_v8_tissue.tsv")
gwas_dir   <- file.path(opt$outdir, "results", opt$gwas, "magma")
genes_raw  <- file.path(gwas_dir, "magma_gene_level.genes.raw")
genes_out  <- file.path(gwas_dir, "magma_gene_level.genes.out")
out_dir    <- file.path(gwas_dir, "interaction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_results <- file.path(out_dir, "interaction_results.tsv")
out_calib   <- file.path(out_dir, "interaction_calibration.tsv")
out_cmds    <- file.path(out_dir, "interaction_commands.json")

for (p in c(magma_bin, tissue_tsv, genes_raw, genes_out)) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}

cat(sprintf("magma_interaction: gwas=%s  min_set_size=%d  fdr_threshold=%g\n",
            opt$gwas, min_set_size, fdr_thr))

# ---------------------------------------------------------------------
# Load .genes.out for Z-scores (used by the outlier follow-up test) and
# the intersection universe. .genes.raw and .genes.out are produced by
# the same MAGMA gene-analysis run - .genes.out has the ZSTAT column.
# ---------------------------------------------------------------------
gres <- fread(genes_out)
gene_z <- setNames(as.numeric(gres$ZSTAT), as.character(gres$GENE))

# Tissue covariate. Columns: entrez, 55 tissues, Average.
tv <- fread(tissue_tsv)
setnames(tv, 1, "entrez")
tv[, entrez := as.character(entrez)]
if (!"Average" %in% names(tv)) stop("GTEx tissue file has no Average column")
tissue_cols <- setdiff(names(tv), c("entrez", "Average"))
if (length(tissue_cols) == 0L) stop("No tissue columns found in ", tissue_tsv)

# Interaction analysis drops genes with any missing expression value
# (spec: no median imputation). Restrict to complete rows now.
tv_complete <- tv[complete.cases(tv)]
gene_universe <- intersect(names(gene_z), tv_complete$entrez)
cat(sprintf("magma_interaction: gene universe size = %d (genes in .genes.out AND non-missing across GTEx)\n",
            length(gene_universe)))

# MAGMA v1.10 rejects the `condition-interaction=Average,S` + `interaction-each=S`
# combination with "conditioned-on variables cannot be used for interaction analysis"
# (S can't be both conditioned-on and the interaction target). To include the
# spec's `overall_expression x S` covariate we must therefore build it ourselves:
# an extra gene-covar column Avg_x_S = Average * 1{g in S}, added as
# `condition-hide=Average,Avg_x_S` to every primary MAGMA call. See plan file.
AVG_INT_COL <- "Avg_x_S"

# ---------------------------------------------------------------------
# Parse every GMT under pathway_gmt_dir, filter each set to the gene
# universe, keep sets with >= min_set_size retained genes.
# ---------------------------------------------------------------------
gmt_files <- sort(list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE))
if (length(gmt_files) == 0L) stop("No .gmt files under ", gmt_dir)

parse_gmt <- function(gf) {
  lines <- readLines(gf)
  rows <- lapply(lines, function(ln) {
    fs <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(fs) < 3L) return(NULL)
    list(name = fs[1], genes = fs[-c(1, 2)])
  })
  rows[!vapply(rows, is.null, logical(1))]
}

sets <- list()
for (gf in gmt_files) {
  gm <- tools::file_path_sans_ext(basename(gf))
  parsed <- parse_gmt(gf)
  for (r in parsed) {
    kept <- intersect(r$genes, gene_universe)
    if (length(kept) < min_set_size) next
    sets[[length(sets) + 1L]] <- list(
      gmt      = gm,
      name     = r$name,
      genes    = kept,
      n_in_set = length(kept),
      n_raw    = length(r$genes)
    )
  }
}
cat(sprintf("magma_interaction: %d sets retained across %d gmt file(s) at min_set_size=%d\n",
            length(sets), length(gmt_files), min_set_size))

if (length(sets) == 0L) {
  # Emit empty schema-compliant outputs so the rule succeeds and
  # downstream aggregators do not crash. This is a valid outcome when
  # every set drops below min_set_size (e.g. chr-only test GWAS with
  # aggressive filters).
  empty_results <- data.table(
    gwas = character(), tissue = character(), gene_set = character(), gmt = character(),
    n_genes_in_set = integer(), n_genes_tested = integer(),
    beta = numeric(), se = numeric(),
    p_interaction = numeric(), p_interaction_outliers_removed = numeric(),
    n_outliers = integer(), p_top25 = numeric(),
    fdr = numeric(), p_bonferroni_tissue = numeric(), retained = logical()
  )
  empty_calib <- data.table(
    tissue = character(), n_tests = integer(), lambda_gc = numeric(),
    median_p = numeric(), ks_p_uniform = numeric()
  )
  fwrite(empty_results, out_results, sep = "\t")
  fwrite(empty_calib,   out_calib,   sep = "\t")
  writeLines(jsonlite::toJSON(list(), auto_unbox = TRUE), out_cmds)
  cat("magma_interaction: 0 sets to test - emitted empty outputs.\n")
  quit(status = 0)
}

# ---------------------------------------------------------------------
# Helpers: temp workspace, single-set GMT writer, MAGMA runner.
# ---------------------------------------------------------------------
# Auto-select a fast temp root: node-local /dev/shm (RAM-backed) if
# writable, else R's tempdir(). Users can override via --tmpdir (or the
# magma_interaction.tmpdir config key) - critical on shared-FS clusters
# where CephFS/Lustre I/O contention balloons per-set MAGMA wall time
# from ~8s (local) to ~70s (shared) under 8-way concurrency.
resolve_tmpdir <- function(user_tmpdir) {
  if (!is.null(user_tmpdir) && !is.na(user_tmpdir) && nzchar(user_tmpdir) &&
      !identical(tolower(user_tmpdir), "auto")) {
    dir.create(user_tmpdir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(user_tmpdir) || file.access(user_tmpdir, mode = 2) != 0L) {
      stop("--tmpdir=", user_tmpdir, " is not a writable directory")
    }
    return(user_tmpdir)
  }
  # Auto: prefer /dev/shm; else fall back to R's tempdir().
  for (cand in c("/dev/shm", tempdir())) {
    if (dir.exists(cand) && file.access(cand, mode = 2) == 0L) return(cand)
  }
  tempdir()
}
# Config file may also set the tmpdir; CLI takes precedence.
cfg_tmpdir <- if (!is.null(mi$tmpdir) && !identical(mi$tmpdir, "NA")) mi$tmpdir else NA_character_
tmpdir_root <- resolve_tmpdir(if (!is.null(opt$tmpdir) && !is.na(opt$tmpdir))
                                opt$tmpdir else cfg_tmpdir)
tmp_root <- tempfile("magma_interaction_", tmpdir = tmpdir_root)
dir.create(tmp_root)
on.exit(unlink(tmp_root, recursive = TRUE), add = TRUE)
cat(sprintf("magma_interaction: per-set MAGMA temp workspace = %s\n", tmp_root))

write_gmt <- function(path, sets_list) {
  # sets_list is a list of list(name=, genes=)
  lines <- vapply(sets_list, function(s) paste(c(s$name, "-", s$genes), collapse = "\t"),
                  FUN.VALUE = character(1))
  writeLines(lines, path)
}

# Build a per-set gene-covar TSV with the tissue columns unchanged plus one
# extra column Avg_x_S = Average * 1{g in S}. Uncentered - the interaction
# test is invariant to the S-indicator centering of the interaction column.
write_perset_covar <- function(path, set_genes) {
  is_in <- tv_complete$entrez %in% set_genes
  perset <- copy(tv_complete)
  perset[, (AVG_INT_COL) := ifelse(is_in, get("Average"), 0)]
  fwrite(perset, path, sep = "\t")
}

run_magma <- function(args_vec, prefix) {
  # Prepend --gene-results and append --out; return command + status.
  args_full <- c("--gene-results", genes_raw, args_vec, "--out", prefix)
  status <- suppressWarnings(system2(magma_bin, args = args_full,
                                     stdout = paste0(prefix, ".stdout"),
                                     stderr = paste0(prefix, ".stderr")))
  list(cmd = paste(shQuote(c(magma_bin, args_full)), collapse = " "),
       status = status)
}

read_gsa <- function(prefix) {
  f <- paste0(prefix, ".gsa.out")
  if (!file.exists(f)) return(NULL)
  tryCatch(fread(cmd = paste("grep -v '^#'", shQuote(f))),
           error = function(e) NULL)
}

# analyse=list requires a specific comma-joined form: "list,name1,name2,...".
# The first item ("list") is the mode; the rest are the tested variables.
tissue_analyse <- paste(tissue_cols, collapse = ",")

# Parse INTER-SC row -> tissue name. MAGMA emits FULL_NAME as
# "INTERACT::<set>::<tissue>"; strip the prefix to recover the tissue.
parse_inter_tissue <- function(d) {
  # MAGMA emits FULL_NAME as "INTERACT::<set>::<tissue>"; strip the prefix.
  raw <- if ("FULL_NAME" %in% names(d)) as.character(d$FULL_NAME) else as.character(d$VARIABLE)
  sub("^INTERACT::.+?::", "", raw)
}

# ---------------------------------------------------------------------
# Primary loop: one MAGMA invocation per (retained set). Each invocation
# uses interaction-each={S} to produce one t x S test per tissue.
# The Average main effect AND the Average x S interaction MUST be
# included as covariates - see the spec's non-negotiable check.
# Parallelised over sets via parallel::mclapply when --n_cores > 1.
# ---------------------------------------------------------------------
n_cores <- if (is.null(opt$n_cores)) 1L else max(1L, as.integer(opt$n_cores))
cat(sprintf("magma_interaction: using %d parallel worker(s) for the per-set MAGMA loop\n",
            n_cores))

run_primary <- function(i) {
  s <- sets[[i]]
  set_gmt_path <- file.path(tmp_root, sprintf("set_%05d.gmt", i))
  covar_path   <- file.path(tmp_root, sprintf("cov_%05d.tsv", i))
  write_gmt(set_gmt_path, list(list(name = s$name, genes = s$genes)))
  write_perset_covar(covar_path, s$genes)
  prefix <- file.path(tmp_root, sprintf("prim_%05d", i))

  args <- c(
    "--set-annot",   set_gmt_path,
    "--gene-covar",  covar_path, "missing-values=drop",
    "--model",
      sprintf("interaction-each=%s", s$name),
      sprintf("condition-hide=Average,%s", AVG_INT_COL),
      "direction-interaction=greater",
      sprintf("interaction-sc-size=%d", min_set_size),
      sprintf("analyse=list,%s,%s", s$name, tissue_analyse)
  )
  r <- run_magma(args, prefix)
  # Free the covar file immediately; keep .gmt and .gsa.out for now.
  unlink(covar_path)

  cmd_entry <- list(
    stage = "primary", set_index = i, gmt = s$gmt, set = s$name,
    n_in_set = s$n_in_set, cmd = r$cmd, status = r$status
  )

  # Always emit one row per tissue for this set (spec: do not silently
  # drop anything). Rows where MAGMA did not return an INTER-SC entry
  # are filled with NA.
  skeleton <- data.table(
    gwas           = opt$gwas,
    tissue         = tissue_cols,
    gene_set       = s$name,
    gmt            = s$gmt,
    n_genes_in_set = s$n_in_set,
    n_genes_tested = NA_integer_,
    beta           = NA_real_,
    se             = NA_real_,
    p_interaction  = NA_real_
  )

  if (r$status != 0L) {
    return(list(cmd = cmd_entry, row = skeleton,
                warn = sprintf("MAGMA failed for set %s (status %d); see %s.stderr",
                               s$name, r$status, prefix)))
  }

  d <- read_gsa(prefix)
  if (is.null(d) || nrow(d) == 0L) return(list(cmd = cmd_entry, row = skeleton))
  # INTER-SC = "interaction: set by covariate" per MAGMA manual p.24
  d <- d[TYPE == "INTER-SC"]
  if (nrow(d) == 0L) return(list(cmd = cmd_entry, row = skeleton))
  d[, tissue_name := parse_inter_tissue(d)]
  # Left-merge onto the tissue skeleton: any tissue MAGMA skipped stays NA.
  m <- match(skeleton$tissue, d$tissue_name)
  ok <- !is.na(m)
  skeleton$n_genes_tested[ok] <- as.integer(d$NGENES[m[ok]])
  skeleton$beta[ok]           <- as.numeric(d$BETA[m[ok]])
  skeleton$se[ok]             <- as.numeric(d$SE[m[ok]])
  skeleton$p_interaction[ok]  <- as.numeric(d$P[m[ok]])
  list(cmd = cmd_entry, row = skeleton)
}

par_lapply <- function(x, FUN) {
  if (n_cores > 1L) {
    parallel::mclapply(x, FUN, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(x, FUN)
  }
}

primary_results <- par_lapply(seq_along(sets), run_primary)
commands <- lapply(primary_results, `[[`, "cmd")
for (pr in primary_results) if (!is.null(pr$warn)) warning(pr$warn, call. = FALSE)
primary <- lapply(primary_results, `[[`, "row")

results <- rbindlist(Filter(Negate(is.null), primary), use.names = TRUE, fill = TRUE)
if (nrow(results) == 0L) {
  warning("No interaction rows produced by MAGMA - check per-set stderr files under ", tmp_root)
}

# Drop any tissue row that isn't in our tissue_cols list (safety - MAGMA
# should only emit INTER-SC for the tested pool but be defensive).
results <- results[tissue %in% tissue_cols]

# ---------------------------------------------------------------------
# Multiple-testing: BH-FDR across the full (tissue x set) grid within
# this GWAS (primary), and per-tissue Bonferroni for parity with the
# 2018 paper.
# ---------------------------------------------------------------------
if (nrow(results) > 0L) {
  results[, fdr := p.adjust(p_interaction, method = "BH")]
  results[, p_bonferroni_tissue := pmin(1, p_interaction * .N), by = tissue]
} else {
  results[, fdr := numeric()]
  results[, p_bonferroni_tissue := numeric()]
}

# ---------------------------------------------------------------------
# Calibration summary per tissue. Convert one-sided P to a z under H0
# (upper-tail) and compute a GC-lambda style ratio of median(z^2) to
# the 1-df chi-square median (0.4549). Also report the raw p median
# (0.5 under H0) and a KS test vs Uniform(0,1).
# ---------------------------------------------------------------------
chi_med_1df <- qchisq(0.5, df = 1)
calib <- results[, {
  p_ok <- p_interaction[!is.na(p_interaction) & p_interaction > 0 & p_interaction < 1]
  list(
    n_tests      = .N,
    lambda_gc    = if (length(p_ok) > 0)
                     median((qnorm(1 - p_ok))^2) / chi_med_1df
                   else NA_real_,
    median_p     = if (length(p_ok) > 0) median(p_ok) else NA_real_,
    ks_p_uniform = if (length(p_ok) > 1)
                     suppressWarnings(ks.test(p_ok, "punif")$p.value)
                   else NA_real_
  )
}, by = tissue]
setkey(calib, tissue)

# ---------------------------------------------------------------------
# Follow-up tests, only for rows passing primary FDR. Both are stored
# side-by-side; neither is used to drop rows.
# ---------------------------------------------------------------------

# Helper: standardise a numeric vector to unit SD within a set.
std <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# --- Follow-up 1: top-25% test ---------------------------------------
run_top25 <- function(s, tissue_name) {
  set_genes <- s$genes
  sub <- tv_complete[entrez %in% set_genes,
                     .(entrez,
                       tv = get(tissue_name),
                       av = get("Average"))]
  if (nrow(sub) < 4L) return(NA_real_)
  resid_t <- residuals(lm(tv ~ av, data = sub))
  ord <- order(resid_t, decreasing = TRUE)
  n_top <- max(1L, floor(0.25 * nrow(sub)))
  top_genes <- sub$entrez[ord[seq_len(n_top)]]
  top_name  <- paste0(s$name, "_TOP25")

  two_gmt <- file.path(tmp_root, sprintf("top25_%s_%s.gmt",
                                         substr(gsub("[^A-Za-z0-9]", "_", s$name), 1, 40),
                                         substr(gsub("[^A-Za-z0-9]", "_", tissue_name), 1, 30)))
  write_gmt(two_gmt, list(list(name = s$name,   genes = set_genes),
                          list(name = top_name, genes = top_genes)))
  prefix <- sub("\\.gmt$", "", two_gmt)
  # analyse=sets restricts targets to set variables; condition-hide=S,tissue
  # keeps the full set S and the tissue covariate in the model as hidden
  # covariates, leaving S_top25 as the only tested set.
  args <- c(
    "--set-annot",  two_gmt,
    "--gene-covar", tissue_tsv, "missing-values=drop",
    "--model",
      "analyse=sets",
      sprintf("condition-hide=%s,%s", s$name, tissue_name),
      "direction-sets=greater"
  )
  r <- run_magma(args, prefix)
  cmd_entry <- list(
    stage = "top25", set = s$name, tissue = tissue_name,
    cmd = r$cmd, status = r$status
  )
  if (r$status != 0L) return(list(p = NA_real_, cmd = cmd_entry))
  d <- read_gsa(prefix)
  if (is.null(d) || nrow(d) == 0L) return(list(p = NA_real_, cmd = cmd_entry))
  # The row corresponding to top_name; match on FULL_NAME to sidestep 30-char truncation.
  nm <- if ("FULL_NAME" %in% names(d)) d$FULL_NAME else d$VARIABLE
  hit <- which(as.character(nm) == top_name)
  if (length(hit) == 0L) return(list(p = NA_real_, cmd = cmd_entry))
  list(p = as.numeric(d$P[hit[1]]), cmd = cmd_entry)
}

# --- Follow-up 2: neighbourhood-peel outlier removal -----------------
# A gene g in S is an outlier if:
#   (a) its distance from the origin in the (residualised-tissue, gene-Z)
#       plane, with each axis standardised to unit SD within S, is > 2, AND
#   (b) every non-outlier gene within Euclidean distance 2 of g is at
#       distance >= d_g from the origin (i.e. equally or more extreme).
# Iterate outside-in until stable.
mark_outliers <- function(coords) {
  # coords: n x 2 matrix, columns (x = resid_t_std, y = z_std)
  d0 <- sqrt(coords[, 1]^2 + coords[, 2]^2)
  n  <- nrow(coords)
  outlier <- rep(FALSE, n)
  candidates <- which(d0 > 2)
  if (length(candidates) == 0L) return(outlier)
  # Precompute pairwise distances once (n typically <= a few thousand).
  D <- as.matrix(dist(coords))
  repeat {
    changed <- FALSE
    # Peel outside-in: try the farthest-from-origin candidates first.
    ord <- candidates[order(-d0[candidates])]
    for (g in ord) {
      if (outlier[g]) next
      nb <- which(!outlier & D[g, ] <= 2)
      nb <- setdiff(nb, g)
      if (length(nb) == 0L || all(d0[nb] >= d0[g])) {
        outlier[g] <- TRUE
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  outlier
}

run_outlier_check <- function(s, tissue_name) {
  set_genes <- s$genes
  sub <- tv_complete[entrez %in% set_genes,
                     .(entrez,
                       tv = get(tissue_name),
                       av = get("Average"))]
  sub <- sub[entrez %in% names(gene_z)]
  if (nrow(sub) < 10L) return(list(p = NA_real_, n_out = NA_integer_))
  sub[, z := gene_z[entrez]]
  resid_t <- residuals(lm(tv ~ av, data = sub))
  coords <- cbind(std(resid_t), std(sub$z))
  is_out <- mark_outliers(coords)
  n_out  <- sum(is_out)
  kept   <- sub$entrez[!is_out]
  if (length(kept) < min_set_size) return(list(p = NA_real_, n_out = as.integer(n_out)))

  set_gmt_path <- file.path(tmp_root, sprintf("out_%s_%s.gmt",
                                              substr(gsub("[^A-Za-z0-9]", "_", s$name), 1, 40),
                                              substr(gsub("[^A-Za-z0-9]", "_", tissue_name), 1, 30)))
  covar_path <- paste0(sub("\\.gmt$", "", set_gmt_path), ".covar.tsv")
  write_gmt(set_gmt_path, list(list(name = s$name, genes = kept)))
  write_perset_covar(covar_path, kept)
  prefix <- sub("\\.gmt$", "", set_gmt_path)
  args <- c(
    "--set-annot",  set_gmt_path,
    "--gene-covar", covar_path, "missing-values=drop",
    "--model",
      sprintf("interaction-each=%s", s$name),
      sprintf("condition-hide=Average,%s", AVG_INT_COL),
      "direction-interaction=greater",
      sprintf("interaction-sc-size=%d", min_set_size),
      sprintf("analyse=list,%s,%s", s$name, tissue_analyse)
  )
  r <- run_magma(args, prefix)
  cmd_entry <- list(
    stage = "outlier", set = s$name, tissue = tissue_name,
    n_outliers = as.integer(n_out),
    cmd = r$cmd, status = r$status
  )
  na_ret <- function() list(p = NA_real_, n_out = as.integer(n_out), cmd = cmd_entry)
  if (r$status != 0L) return(na_ret())
  d <- read_gsa(prefix)
  if (is.null(d) || nrow(d) == 0L) return(na_ret())
  d <- d[TYPE == "INTER-SC"]
  if (nrow(d) == 0L) return(na_ret())
  d[, .tissue := parse_inter_tissue(d)]
  hit <- which(d$.tissue == tissue_name)
  if (length(hit) == 0L) return(na_ret())
  list(p = as.numeric(d$P[hit[1]]), n_out = as.integer(n_out), cmd = cmd_entry)
}

# Initialise follow-up columns with NA
results[, p_top25 := NA_real_]
results[, p_interaction_outliers_removed := NA_real_]
results[, n_outliers := NA_integer_]

if (nrow(results) > 0L) {
  # Index sets by (gmt, gene_set) for lookup in follow-ups
  set_lookup <- setNames(sets, vapply(sets, function(s) paste(s$gmt, s$name, sep = "|"),
                                      FUN.VALUE = character(1)))
  # Fire follow-up diagnostics on any row with p_interaction below
  # `followup_p_threshold` (default 1e-4). This is broader than primary FDR:
  # under set-level bias many nominally-small p-values can appear without any
  # row passing BH-FDR, and the top-25% / outlier tests are exactly the
  # discriminator between real tissue-specific effects and set-level leakage.
  fdr_hits <- which(!is.na(results$p_interaction) &
                    results$p_interaction < followup_p_thr)
  cat(sprintf("magma_interaction: %d row(s) with p_interaction < %g; running follow-up tests\n",
              length(fdr_hits), followup_p_thr))
  run_followup <- function(k) {
    row_gmt <- results$gmt[k]
    row_set <- results$gene_set[k]
    row_tis <- results$tissue[k]
    s <- set_lookup[[paste(row_gmt, row_set, sep = "|")]]
    if (is.null(s)) return(NULL)
    top <- tryCatch(run_top25(s, row_tis),
                    error = function(e) list(p = NA_real_,
                                             cmd = list(stage = "top25", set = s$name,
                                                        tissue = row_tis,
                                                        error = conditionMessage(e))))
    out <- tryCatch(run_outlier_check(s, row_tis),
                    error = function(e) list(p = NA_real_, n_out = NA_integer_,
                                             cmd = list(stage = "outlier", set = s$name,
                                                        tissue = row_tis,
                                                        error = conditionMessage(e))))
    list(k = k,
         p_top25 = top$p,
         p_outliers_removed = out$p,
         n_out = out$n_out,
         cmds = list(top$cmd, out$cmd))
  }
  followups <- par_lapply(fdr_hits, run_followup)
  for (fu in followups) {
    if (is.null(fu)) next
    k <- fu$k
    results$p_top25[k]                        <- fu$p_top25
    results$p_interaction_outliers_removed[k] <- fu$p_outliers_removed
    results$n_outliers[k]                     <- fu$n_out
    for (ce in fu$cmds) if (!is.null(ce)) commands[[length(commands) + 1L]] <- ce
  }
}

results[, retained := !is.na(fdr) & fdr < fdr_thr &
        ((!is.na(p_top25)                        & p_top25 < 0.05) |
         (!is.na(p_interaction_outliers_removed) & p_interaction_outliers_removed < 0.05))]

# Set-level bias flag: count the number of tissues (out of n_tissues) at which
# this set fires below `bias_flag_p_threshold`. Real tissue-specific effects
# concentrate in a small number of biologically related tissues; set-level
# leakage (see plan + de Leeuw 2018 top-25% discussion) spreads across many
# unrelated tissues. `set_bias_suspect` is a coarse but informative flag.
if (nrow(results) > 0L) {
  results[, set_n_tissues_p_low :=
            sum(!is.na(p_interaction) & p_interaction < bias_p_thr),
          by = .(gmt, gene_set)]
  results[, set_bias_suspect := set_n_tissues_p_low >= bias_min_hits]
} else {
  results[, set_n_tissues_p_low := integer()]
  results[, set_bias_suspect    := logical()]
}

# ---------------------------------------------------------------------
# Write outputs (fixed column order matching the spec).
# ---------------------------------------------------------------------
cols <- c("gwas", "tissue", "gene_set", "gmt",
          "n_genes_in_set", "n_genes_tested",
          "beta", "se",
          "p_interaction", "p_interaction_outliers_removed",
          "n_outliers", "p_top25",
          "fdr", "p_bonferroni_tissue", "retained",
          "set_n_tissues_p_low", "set_bias_suspect")
fwrite(results[, ..cols], out_results, sep = "\t")
fwrite(calib, out_calib, sep = "\t")
writeLines(jsonlite::toJSON(commands, auto_unbox = TRUE, pretty = FALSE), out_cmds)

cat(sprintf("magma_interaction: wrote %d rows -> %s\n", nrow(results), out_results))
cat(sprintf("magma_interaction: wrote %d tissue-calibration rows -> %s\n", nrow(calib), out_calib))
cat("magma_interaction: per-tissue GC-lambda summary (all should be near 1.0 under the null):\n")
print(calib[order(-abs(lambda_gc - 1))], nrows = 60)

cat(sprintf("magma_interaction: done in %.1f s\n",
            as.numeric(difftime(Sys.time(), start.time, units = "secs"))))
