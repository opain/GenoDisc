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
long <- build_comparison_long(gd, all9, c("tissue", "atc", "gene", "locus", "drug"))

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
# Test 7: gene entity type — PWAS-FUSION empty-ASD column invariant.
#
# ASD has zero FDR-significant PWAS-FUSION rows in this bundle, but pivoting
# a PWAS-FUSION slice across all 9 GWAS must still produce a 9-column matrix
# with ASD present as an all-NA column. This is the invariant the plan
# specifically calls out.

message("Test 7 — Gene layer PWAS-FUSION empty-column invariant")

pwas <- long[method == "PWAS-FUSION"]
expect(nrow(pwas) > 0, "PWAS-FUSION rows must exist in the long tibble")

# Sig PWAS-FUSION per GWAS - ASD should be zero
sig_by_gwas <- pwas[fdr < 0.05, .N, by = gwas]
asd_sig <- sig_by_gwas[gwas == "ASD", N]
expect(length(asd_sig) == 0 || asd_sig == 0,
       sprintf("ASD FDR-sig PWAS-FUSION rows = %s (expected 0)",
               if (length(asd_sig) == 0) "0" else as.character(asd_sig)))
pass("ASD has 0 FDR-sig PWAS-FUSION genes")

# Pivot the sig-only slice and confirm ASD column is preserved (all NA)
sig_pwas <- pick_best_per_cell(pwas[fdr < 0.05], c("gwas", "entity_id"))
m <- pivot_matrix(sig_pwas, "fdr", all9)
expect(ncol(m) == 9,
       sprintf("Sig-PWAS pivot ncol = %d (expected 9)", ncol(m)))
expect(all(is.na(m[, "ASD"])),
       "ASD column should be all-NA in the sig-only pivoted matrix")
pass("Sig-PWAS pivot preserves 9 columns; ASD column is all-NA")

# ---------------------------------------------------------------------------
# Test 8: gene entity type — chr17q21.31 recurrence at k=4 for TWAS-FUSION.
#
# MAPT and the surrounding chr17q21.31 genes (ARHGAP27, DND1, KANSL1,
# LRRC37A4P) are FDR-significant in ALZ, ASD, MDD, PRK, SCZ under
# TWAS-FUSION. Plan text listed 4 traits; actual bundle shows 5 (MDD
# additionally). All 5 genes should therefore survive a k >= 4 recurrence
# filter.

message("Test 8 — Gene layer TWAS-FUSION chr17q21.31 recurrence (k >= 4)")

chr17 <- c("MAPT", "ARHGAP27", "DND1", "KANSL1", "LRRC37A4P")
twas_chr17 <- long[method == "TWAS-FUSION" & entity_id %in% chr17]
best <- pick_best_per_cell(twas_chr17, c("gwas", "entity_id"))
recurrence <- best[, .(k = sum(fdr < 0.05, na.rm = TRUE)), by = entity_id]
for (g in chr17) {
  actual_k <- recurrence[entity_id == g, k]
  expect(length(actual_k) == 1 && actual_k >= 4,
         sprintf("%s recurrence k = %s (expected >= 4)",
                 g, if (length(actual_k) == 0) "0" else as.character(actual_k)))
  pass(sprintf("%s FDR-sig in %d GWAS", g, actual_k))
}

# Apply the recurrence filter at k=4 - all 5 chr17q21.31 genes must survive.
filt <- apply_comparison_filters(best,
  list(sig_basis = "fdr", sig_threshold = 0.05, k_min = 4L))
survivors <- unique(filt$entity_id)
expect(setequal(intersect(chr17, survivors), chr17),
       sprintf("Survivors of k=4 filter = %s (expected superset of %s)",
               paste(sort(survivors), collapse = ","),
               paste(sort(chr17), collapse = ",")))
pass("All 5 chr17q21.31 genes survive k=4 filter")

# ---------------------------------------------------------------------------
# Test 9: locus entity type — clump and COJO builders.

message("Test 9 — Locus layer (clump / COJO)")

for (m in c("clump", "COJO")) {
  n_rows <- long[method == m & entity_type == "locus", .N]
  expect(n_rows > 0, sprintf("%s locus rows > 0 (got %d)", m, n_rows))
  pass(sprintf("%s locus rows = %d", m, n_rows))
}

# The literal string "None" in the NearestGene column must not be used as a
# shared entity_id (that would collapse unrelated SNPs across GWAS). The
# builder falls back to the SNP identifier for those rows.
none_rows <- long[entity_id == "None" & entity_type == "locus", .N]
expect(none_rows == 0, sprintf("Locus rows with entity_id == 'None' = %d (expected 0)", none_rows))
pass("No locus rows use 'None' as entity_id")

# Cross-GWAS recurrence at k=4 on clump surfaces at least KANSL1 (chr17q21.31)
# which is FDR-sig in the gene layer for the same 5 traits.
clump <- long[method == "clump" & entity_type == "locus"]
best <- pick_best_per_cell(clump, c("gwas", "entity_id"))
recur <- best[, .(k = .N), by = entity_id][k >= 4]
expect("KANSL1" %in% recur$entity_id,
       "KANSL1 should appear as a clumped-locus hit in >= 4 GWAS")
pass(sprintf("KANSL1 recurrence >= 4 GWAS (actual k = %d)",
             recur$k[recur$entity_id == "KANSL1"]))

# Column-preservation invariant for the locus layer.
m <- pivot_matrix(best, "p", all9)
expect(ncol(m) == 9,
       sprintf("Locus pivot ncol = %d (expected 9)", ncol(m)))
pass(sprintf("Locus pivot preserves %d columns", ncol(m)))

# ---------------------------------------------------------------------------
# Test 10: drug entity type — MAGMA and TWAS-GSEA builders.

message("Test 10 — Drug layer (MAGMA / TWAS-GSEA)")

for (m in c("MAGMA-drug", "TWAS-GSEA-drug")) {
  n_rows <- long[method == m & entity_type == "drug", .N]
  expect(n_rows > 0, sprintf("%s drug rows > 0 (got %d)", m, n_rows))
  pass(sprintf("%s rows = %d", m, n_rows))
}

# Direction is carried through verbatim for TWAS-GSEA-drug and only there.
gsea_dir <- long[method == "TWAS-GSEA-drug", unique(direction)]
expect(all(c("Matches disease", "Opposes disease") %in% gsea_dir),
       sprintf("TWAS-GSEA-drug Direction values = %s (expected includes 'Matches disease' + 'Opposes disease')",
               paste(sort(gsea_dir), collapse = ",")))
pass("TWAS-GSEA-drug Direction populated with expected labels")

magma_dir <- long[method == "MAGMA-drug", unique(direction)]
expect(length(magma_dir) == 1 && is.na(magma_dir),
       sprintf("MAGMA-drug Direction should be all-NA (got %s)",
               paste(magma_dir, collapse = ",")))
pass("MAGMA-drug Direction is NA (undirected)")

# Column-preservation invariant on the drug layer.
best <- pick_best_per_cell(long[method == "TWAS-GSEA-drug"], c("gwas", "entity_id"))
m <- pivot_matrix(best, "fdr", all9)
expect(ncol(m) == 9,
       sprintf("Drug pivot ncol = %d (expected 9)", ncol(m)))
pass(sprintf("Drug pivot preserves %d columns", ncol(m)))

# ---------------------------------------------------------------------------

message("\nAll assertions passed.")
