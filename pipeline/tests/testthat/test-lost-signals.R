# Tests for the lost-signal locus logic (scripts/functions/lost_signals_functions.R).
# Run with:  Rscript -e "testthat::test_dir('tests/testthat')"  from pipeline/.
#
# Regression focus: the lost-signal check must compare removed vs retained variants on the SAME
# genome build. The bug this guards against fed removed positions in the input build (e.g. GRCh38)
# against retained positions remapped to GRCh37, so nearby retained GW-sig variants fell outside the
# ±500 kb window and loci were falsely flagged lost.

suppressMessages(library(data.table))

# Locate the pipeline dir and source the function under test. wu_pipeline_dir() comes from
# helper-writeup.R (auto-loaded by testthat::test_dir) and resolves the same pipeline directory.
.ls_pd <- tryCatch(wu_pipeline_dir(), error = function(e) {
  d <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) {
    if (file.exists(file.path(d, "scripts", "functions", "lost_signals_functions.R"))) return(d)
    d <- dirname(d)
  }
  stop("could not locate pipeline dir (scripts/functions/lost_signals_functions.R)")
})
source(file.path(.ls_pd, "scripts", "functions", "lost_signals_functions.R"))

# A GW-significant removed lead, recorded in the input build (GRCh38-like coordinate). Mirrors the
# real chr21 example: removed SNP P=2.6e-35 at GRCh38 34,233,898 (GRCh37 35,606,199).
removed_one <- function() {
  data.table(SNP = "rs_removed", CHR = "21", BP = 34233898, P = 2.6e-35, reason = "maf_discordance")
}

test_that("retained GW-sig variant in the window (same build) is NOT a lost signal", {
  removed <- removed_one()
  # Retained GW-sig SNP 200 kb away, in the SAME (input/GRCh38) build as the removed file.
  retained <- data.table(SNP = "rs_kept", CHR = "21", BP = 34233898 + 2e5, P = 1e-20)

  loci <- compute_lost_signal_loci(removed, retained)
  expect_equal(nrow(loci), 1L)
  expect_false(loci$lost_signal[1])
  expect_equal(loci$n_retained_sig_in_window[1], 1L)
})

test_that("the cross-build bug: GRCh37 retained vs GRCh38 removed falsely reads as lost", {
  removed <- removed_one()
  # Same retained variant but given in GRCh37 coordinates (the buggy mix): ~1.37 Mb away from the
  # GRCh38 removed lead, so it falls outside the ±500 kb window.
  retained_grch37 <- data.table(SNP = "rs_kept", CHR = "21", BP = 35606199 + 2e5, P = 1e-20)

  loci <- compute_lost_signal_loci(removed, retained_grch37)
  expect_equal(nrow(loci), 1L)
  expect_true(loci$lost_signal[1])   # demonstrates what the fix prevents
})

test_that("retained variant beyond the window is a lost signal", {
  removed <- removed_one()
  retained <- data.table(SNP = "rs_far", CHR = "21", BP = 34233898 + 7e5, P = 1e-20)  # 700 kb away
  loci <- compute_lost_signal_loci(removed, retained)
  expect_true(loci$lost_signal[1])
})

test_that("retained variant in window but not GW-significant is a lost signal", {
  removed <- removed_one()
  retained <- data.table(SNP = "rs_weak", CHR = "21", BP = 34233898 + 1e5, P = 1e-3)
  loci <- compute_lost_signal_loci(removed, retained)
  expect_equal(loci$n_retained_in_window[1], 1L)
  expect_equal(loci$n_retained_sig_in_window[1], 0L)
  expect_true(loci$lost_signal[1])
})

test_that("removed variant on a different chromosome is not matched to retained signal", {
  removed <- removed_one()
  retained <- data.table(SNP = "rs_otherchr", CHR = "20", BP = 34233898 + 1e5, P = 1e-20)
  loci <- compute_lost_signal_loci(removed, retained)
  expect_true(loci$lost_signal[1])
})

test_that("multiple removed variants within a window collapse into one lead locus", {
  removed <- data.table(
    SNP    = c("rs_lead", "rs_near", "rs_other"),
    CHR    = c("21", "21", "21"),
    BP     = c(34233898, 34233898 + 3e5, 34233898 + 2e6),  # first two share a window, third separate
    P      = c(2.6e-35, 1e-12, 5e-9),
    reason = c("maf_discordance", "maf_discordance", "info")
  )
  retained <- data.table(SNP = "rs_kept", CHR = "21", BP = 34233898 + 1e5, P = 1e-20)
  loci <- compute_lost_signal_loci(removed, retained)
  expect_equal(nrow(loci), 2L)
  lead_locus <- loci[lead_SNP == "rs_lead"]
  expect_equal(lead_locus$n_removed, 2L)   # rs_lead + rs_near
  expect_false(lead_locus$lost_signal)      # retained GW-sig sits in the lead window
})

test_that("no removed GW-sig variants yields an empty loci table (no error)", {
  empty_removed <- data.table(SNP = character(0), CHR = character(0), BP = numeric(0),
                              P = numeric(0), reason = character(0))
  retained <- data.table(SNP = "rs_kept", CHR = "21", BP = 34233898, P = 1e-20)
  expect_equal(nrow(compute_lost_signal_loci(empty_removed, retained)), 0L)
  expect_equal(nrow(compute_lost_signal_loci(NULL, NULL)), 0L)
})

test_that("empty retained set flags every removed locus as lost", {
  removed <- removed_one()
  empty_retained <- data.table(CHR = character(0), BP = numeric(0), P = numeric(0))
  loci <- compute_lost_signal_loci(removed, empty_retained)
  expect_equal(nrow(loci), 1L)
  expect_equal(loci$n_retained_in_window[1], 0L)
  expect_true(loci$lost_signal[1])
})
