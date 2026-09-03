# Data-level smoke tests for the cross-GWAS comparison layer.
#
# Verifies that build_comparison_long() + helpers reproduce the numeric
# assertions from the internal 9-GWAS neuropsych review. This does NOT
# exercise the Shiny UI — that is checked manually by running the app.
#
# Run from the app directory:
#   Rscript tests/test_comparison_long.R
# or under testthat:
#   Rscript -e "testthat::test_file('tests/test_comparison_long.R')"

# ---------------------------------------------------------------------------
# Locate the app source and the reference bundle.

here <- tryCatch(dirname(sys.frame(1L)$ofile), error = function(e) getwd())
app_dir <- normalizePath(file.path(here, ".."), mustWork = FALSE)
if (!file.exists(file.path(app_dir, "reader.R"))) {
  # Fallback: assume cwd is the app dir when running with `Rscript tests/...`
  app_dir <- getwd()
}

BUNDLE <- "/users/k1806347/oliverpainfel/Analyses/genodisc_validation/output_neuropsych/results/results_package.rds"

skip_all <- !file.exists(BUNDLE)

if (skip_all) {
  message("Reference bundle not found at:\n  ", BUNDLE,
          "\nSkipping tests.")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Source the app helpers.

source(file.path(app_dir, "reader.R"))
source(file.path(app_dir, "R", "utils_data.R"))
source(file.path(app_dir, "R", "utils_comparison.R"))

suppressPackageStartupMessages({
  library(data.table)
  if (!requireNamespace("testthat", quietly = TRUE)) {
    message("testthat not installed; running as plain assertions.")
    testthat_available <- FALSE
  } else {
    testthat_available <- TRUE
  }
})

# ---------------------------------------------------------------------------
# Build the long tibble once.

gd <- gd_open(BUNDLE)
all9 <- gd_gwas(gd)
long <- build_comparison_long(gd, all9, c("tissue", "atc"))

expect <- if (testthat_available) testthat::expect_true else function(x, info = "") {
  if (!isTRUE(x)) stop("ASSERT FAILED: ", info, call. = FALSE)
  invisible(TRUE)
}
expect_equal <- if (testthat_available) testthat::expect_equal else function(a, b, tol = 1e-8, info = "") {
  if (isTRUE(all.equal(a, b, tolerance = tol))) invisible(TRUE)
  else stop("EQUAL FAILED: ", info, "\n  actual   = ",
            paste(a, collapse=","), "\n  expected = ", paste(b, collapse=","),
            call. = FALSE)
}

pass <- function(msg) message(sprintf("  ✓ %s", msg))

# ---------------------------------------------------------------------------
# Test 1: Tissues — Cerebellum & Cerebellar Hemisphere recurrence.

message("Test 1 — Tissues (Cerebellum / Cerebellar Hemisphere)")

cere <- c("Brain - Cerebellum", "Brain - Cerebellar Hemisphere")
n_fdr <- long[entity_id %in% cere & fdr < 0.05, .N]
expect(n_fdr == 14, sprintf("Cere/CereH FDR<0.05 rows = %d (expected 14)", n_fdr))
pass(sprintf("Cere/CereH FDR<0.05 rows = %d", n_fdr))

n_ret <- long[entity_id %in% cere & evidence == TRUE, .N]
expect(n_ret == 12, sprintf("Cere/CereH retained rows = %d (expected 12)", n_ret))
pass(sprintf("Cere/CereH retained rows = %d", n_ret))

mig_brain <- long[gwas == "MIG" & grepl("Brain|Cerebell", entity_id) & evidence == TRUE, .N]
expect(mig_brain == 0, sprintf("MIG retained brain tissues = %d (expected 0)", mig_brain))
pass(sprintf("MIG retained brain tissues = %d", mig_brain))

# ---------------------------------------------------------------------------
# Test 2: ATC TWAS-GSEA N05A (antipsychotics) — best-per-cell across panels.

message("Test 2 — ATC TWAS-GSEA N05A best-per-cell")

n05a <- pick_best_per_cell(
  long[method == "TWAS-GSEA-ATC" & entity_id == "N05A"],
  c("gwas", "entity_id")
)
expect(nrow(n05a) == 9, sprintf("N05A rows = %d (expected 9)", nrow(n05a)))
pass(sprintf("Best-per-cell rows = %d", nrow(n05a)))

n_nom <- sum(n05a$p < 0.05)
expect(n_nom == 9, sprintf("N05A nominal-sig = %d (expected 9)", n_nom))
pass(sprintf("Nominal-sig = %d", n_nom))

n_fdr_sig <- sum(n05a$fdr < 0.05)
expect(n_fdr_sig == 5, sprintf("N05A FDR-sig = %d (expected 5)", n_fdr_sig))
pass(sprintf("FDR-sig = %d", n_fdr_sig))

fdr_sig_gwas <- sort(n05a$gwas[n05a$fdr < 0.05])
expect(setequal(fdr_sig_gwas, c("ADHD", "ALS", "ALZ", "MDD", "PRK")),
       sprintf("FDR-sig GWAS = %s (expected ADHD,ALS,ALZ,MDD,PRK)",
               paste(fdr_sig_gwas, collapse = ",")))
pass(sprintf("FDR-sig GWAS = %s", paste(fdr_sig_gwas, collapse = ",")))

n_opposes <- sum(n05a$direction == "Opposes disease", na.rm = TRUE)
n_matches <- sum(n05a$direction == "Matches disease", na.rm = TRUE)
expect(n_opposes == 8 && n_matches == 1,
       sprintf("Direction: opposes=%d matches=%d (expected 8,1)", n_opposes, n_matches))
pass(sprintf("opposes=%d matches=%d", n_opposes, n_matches))

matches_gwas <- n05a$gwas[n05a$direction == "Matches disease"]
expect(identical(matches_gwas, "PRK"),
       sprintf("matches GWAS = %s (expected PRK)", paste(matches_gwas, collapse = ",")))
pass(sprintf("matches GWAS = %s", matches_gwas))

expect_equal(n05a$fdr[n05a$gwas == "PRK"], 4.2e-4, tol = 5e-5,
             info = "PRK FDR")
pass(sprintf("PRK FDR = %.3e (expected ~4.2e-4)", n05a$fdr[n05a$gwas == "PRK"]))

expect_equal(n05a$fdr[n05a$gwas == "MDD"], 0.034, tol = 5e-4,
             info = "MDD FDR")
pass(sprintf("MDD FDR = %.3e (expected ~0.034)", n05a$fdr[n05a$gwas == "MDD"]))

expect_equal(n05a$fdr[n05a$gwas == "ADHD"], 2.9e-4, tol = 5e-5,
             info = "ADHD FDR")
pass(sprintf("ADHD FDR = %.3e (expected ~2.9e-4)", n05a$fdr[n05a$gwas == "ADHD"]))

# ---------------------------------------------------------------------------
# Test 3: ATC MAGMA — SCZ/BIP hallmark classes.

message("Test 3 — ATC MAGMA hallmark cells")

magma <- long[method == "MAGMA-ATC"]

check <- function(gwas_, code, expected, tol_rel = 0.05) {
  actual <- magma[gwas == gwas_ & entity_id == code, fdr]
  expect(length(actual) == 1 && !is.na(actual),
         sprintf("Missing cell: %s x %s", gwas_, code))
  expect(abs(actual - expected) < abs(expected) * tol_rel + 1e-15,
         sprintf("%s x %s FDR = %.3e (expected ~%.3e)", gwas_, code, actual, expected))
  pass(sprintf("%s x %s FDR = %.3e", gwas_, code, actual))
}

check("SCZ", "N05A", 5.4e-9)
check("BIP", "N05C", 4.9e-12)
check("BIP", "N03A", 9.0e-8)
check("BIP", "C08C", 1.7e-4)

# ---------------------------------------------------------------------------
# Test 4: Column-preservation invariant on pivot_matrix().

message("Test 4 — pivot_matrix() column preservation")

# B05C is present in only 5 of the 9 GWAS (ALZ, BIP, MIG, PRK, SCZ).
slice_b05c <- long[method == "TWAS-GSEA-ATC" & entity_id == "B05C"]
best_b05c <- pick_best_per_cell(slice_b05c, c("gwas", "entity_id"))
m <- pivot_matrix(best_b05c, "fdr", all9)

expect(ncol(m) == 9, sprintf("B05C matrix ncol = %d (expected 9)", ncol(m)))
pass(sprintf("B05C matrix ncol = %d (all 9 selected GWAS preserved)", ncol(m)))

not_tested <- colnames(m)[is.na(m)]
expect(setequal(not_tested, c("ADHD", "ALS", "ASD", "MDD")),
       sprintf("B05C not-tested GWAS = %s (expected ADHD,ALS,ASD,MDD)",
               paste(sort(not_tested), collapse = ",")))
pass(sprintf("B05C not-tested GWAS = %s", paste(sort(not_tested), collapse = ",")))

# Empty slice invariant: pivoting an empty long tibble still yields the right
# number of columns.
empty <- long[FALSE]
m_empty <- pivot_matrix(empty, "fdr", all9)
expect(ncol(m_empty) == 9,
       sprintf("Empty pivot ncol = %d (expected 9)", ncol(m_empty)))
pass(sprintf("Empty pivot ncol = %d", ncol(m_empty)))

# ---------------------------------------------------------------------------
# Test 5: apply_comparison_filters() with k_min recurrence filter.

message("Test 5 — apply_comparison_filters recurrence k=2")

# Take the union of ATC classes across the 9 GWAS after best-per-cell
# reduction, then apply k_min=2 on FDR<0.05. Expect only classes that are
# FDR-sig in >=2 GWAS to survive; N05A survives (5 FDR-sig).
best_all <- pick_best_per_cell(
  long[method == "TWAS-GSEA-ATC"], c("gwas", "entity_id"))
filt <- apply_comparison_filters(best_all,
  list(sig_basis = "fdr", sig_threshold = 0.05, k_min = 2L))
expect("N05A" %in% filt$entity_id,
       "N05A should survive k_min=2 filter (FDR-sig in 5 GWAS)")
pass("N05A survives k_min=2")

# Sanity: nothing that is FDR-sig in <2 GWAS should survive.
per_class <- best_all[, .(k = sum(fdr < 0.05, na.rm = TRUE)), by = entity_id]
survivors <- unique(filt$entity_id)
non_survivors <- setdiff(per_class$entity_id, survivors)
non_survivor_max_k <- if (length(non_survivors) > 0)
  max(per_class$k[per_class$entity_id %in% non_survivors]) else 0
expect(non_survivor_max_k < 2,
       sprintf("Some non-surviving class has k=%d (expected all <2)",
               non_survivor_max_k))
pass(sprintf("All non-survivors have k<2 (max non-survivor k = %d)",
             non_survivor_max_k))

# ---------------------------------------------------------------------------
# Test 6: Overview builders don't blow up and produce one row per GWAS.

message("Test 6 — Overview QC and Yield builders")

qc <- build_overview_qc(gd, all9)
expect(nrow(qc) == 9, sprintf("Overview QC rows = %d (expected 9)", nrow(qc)))
expect(all(!is.na(qc$lambda_gc)),
       "Some lambda_GC values are NA (should be populated for all 9 traits)")
pass(sprintf("Overview QC: %d rows, all lambda_GC populated", nrow(qc)))

yield <- build_overview_yield(long, gd, all9, sig_basis = "fdr", sig_threshold = 0.05)
expect(nrow(yield) == 9, sprintf("Overview yield rows = %d (expected 9)", nrow(yield)))
pass(sprintf("Overview yield: %d rows", nrow(yield)))

# ---------------------------------------------------------------------------

message("\nAll assertions passed.")
