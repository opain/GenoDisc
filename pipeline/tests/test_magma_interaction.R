#!/usr/bin/Rscript
# Post-run assertions for the MAGMA gene-set x tissue interaction module.
# Invoked with a single positional arg: the outdir root
# (e.g. "tests/output"). Fails loudly on the first violation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: test_magma_interaction.R <outdir>")
outdir <- args[1]

suppressMessages(library(data.table))
suppressMessages(library(yaml))
suppressMessages(library(jsonlite))

cfg <- yaml::read_yaml("tests/config_test.yaml")
if (!isTRUE(cfg$magma_interaction$enabled)) {
  stop("magma_interaction.enabled is not true in tests/config_test.yaml")
}
min_set_size <- as.integer(cfg$magma_interaction$min_set_size)
gmt_dir      <- cfg$pathway_gmt_dir

gwas_names <- readLines("tests/gwas_list_test.txt")[-1]
gwas_names <- vapply(strsplit(gwas_names, "\\s+"), `[`, character(1), 1)
gwas_names <- gwas_names[nzchar(gwas_names)]

expected_cols <- c(
  "gwas", "tissue", "gene_set", "gmt",
  "n_genes_in_set", "n_genes_tested",
  "beta", "se",
  "p_interaction", "p_interaction_outliers_removed",
  "n_outliers", "p_top25",
  "fdr", "p_bonferroni_tissue", "retained"
)

check_gwas <- function(g) {
  cat(sprintf("[assert] gwas=%s\n", g))
  res_path   <- file.path(outdir, "results", g, "magma", "interaction", "interaction_results.tsv")
  calib_path <- file.path(outdir, "results", g, "magma", "interaction", "interaction_calibration.tsv")
  cmds_path  <- file.path(outdir, "results", g, "magma", "interaction", "interaction_commands.json")
  stopifnot(file.exists(res_path), file.exists(calib_path), file.exists(cmds_path))

  # --- 1. Schema ---------------------------------------------------
  res <- fread(res_path)
  if (!identical(names(res), expected_cols)) {
    stop("schema mismatch:\n  expected: ", paste(expected_cols, collapse = ","),
         "\n  got:      ", paste(names(res), collapse = ","))
  }
  cat(sprintf("  schema OK; %d rows\n", nrow(res)))

  # --- 2. Row count = n_tissues x n_sets_passing --------------------
  # .genes.raw has variable-width rows (correlations appended), so read
  # line-by-line and take the first whitespace-delimited token (gene id).
  genes_raw <- file.path(outdir, "results", g, "magma", "magma_gene_level.genes.raw")
  raw_lines <- readLines(genes_raw)
  raw_lines <- raw_lines[!startsWith(raw_lines, "#")]
  gr_genes <- vapply(strsplit(raw_lines, "\\s+"), `[`, character(1), 1)
  gr_genes <- gr_genes[nzchar(gr_genes)]
  tv <- fread("resources/data/gtex/GTEx_v8_tissue.tsv")
  setnames(tv, 1, "entrez")
  tv[, entrez := as.character(entrez)]
  tv <- tv[complete.cases(tv)]
  gene_universe <- intersect(gr_genes, tv$entrez)
  tissues <- setdiff(names(tv), c("entrez", "Average"))

  n_pass <- 0L
  for (gf in list.files(gmt_dir, "\\.gmt$", full.names = TRUE)) {
    for (ln in readLines(gf)) {
      fs <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(fs) < 3L) next
      kept <- intersect(fs[-c(1, 2)], gene_universe)
      if (length(kept) >= min_set_size) n_pass <- n_pass + 1L
    }
  }
  expected_rows <- length(tissues) * n_pass
  if (nrow(res) != expected_rows) {
    stop("row-count invariant broken: got ", nrow(res),
         "  expected ", expected_rows,
         " (", length(tissues), " tissues x ", n_pass, " sets)")
  }
  cat(sprintf("  row-count invariant OK (%d tissues x %d sets = %d)\n",
              length(tissues), n_pass, expected_rows))

  # --- 3. Every primary MAGMA command must include the Average main-effect
  # AND the Average x S interaction covariate. MAGMA v1.10 rejects the built-in
  # `condition-interaction=Average,S` combined with `interaction-each=S`
  # ("conditioned-on variables cannot be used for interaction analysis"), so
  # the wrapper injects the interaction via an explicit gene-covar column
  # named Avg_x_S and passes `condition-hide=Average,Avg_x_S`. The invariant
  # therefore checks that fixed string.
  cmds <- jsonlite::fromJSON(cmds_path, simplifyDataFrame = FALSE)
  if (n_pass > 0L && length(cmds) == 0L) stop("commands JSON is empty despite ", n_pass, " sets")
  prim <- Filter(function(x) identical(x$stage, "primary"), cmds)
  if (length(prim) != n_pass) {
    stop("primary command count ", length(prim), " != expected sets ", n_pass)
  }
  bad <- Filter(function(x) !grepl("condition-hide=Average,Avg_x_S", x$cmd, fixed = TRUE),
                prim)
  if (length(bad) > 0L) {
    stop("primary MAGMA command missing condition-hide=Average,Avg_x_S: ",
         bad[[1]]$cmd)
  }
  cat(sprintf("  MAGMA command-string check OK (%d primary calls all include condition-hide=Average,Avg_x_S)\n",
              length(prim)))

  # --- 4. Calibration TSV has one row per tested tissue. lambda_gc may be
  # NA for tissues where MAGMA could not compute any p-value (e.g. all
  # models collinear on a chr-only test GWAS) - that is a diagnostic
  # signal, not a failure. Assert at least one tissue produced a finite
  # lambda so the module isn't silently producing empty calibration.
  calib <- fread(calib_path)
  if (nrow(res) > 0L) {
    stopifnot("tissue"    %in% names(calib),
              "lambda_gc" %in% names(calib))
    missing_t <- setdiff(unique(res$tissue), calib$tissue)
    if (length(missing_t) > 0L) {
      stop("tissues missing from calibration: ", paste(missing_t, collapse = ","))
    }
    finite_lam <- calib$lambda_gc[is.finite(calib$lambda_gc)]
    if (length(finite_lam) == 0L) {
      stop("no finite lambda_gc in calibration - MAGMA produced no usable interaction p-values")
    }
    cat(sprintf("  calibration OK; %d/%d tissues with finite lambda_gc, range = [%.3f, %.3f]\n",
                length(finite_lam), nrow(calib),
                min(finite_lam), max(finite_lam)))
  } else {
    cat("  calibration TSV present (results empty)\n")
  }
}

for (g in gwas_names) check_gwas(g)
cat("test_magma_interaction: all assertions passed.\n")
