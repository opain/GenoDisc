# Per-block sentence-generation tests, including null-result and block-not-run cases.

wu_source()

new_ctx <- function() writeup_new_ctx(wu_pipeline_dir(), fx_config(), schema_version = 2L)

test_that("preprocessing paragraph reports QC thresholds and variant counts", {
  ctx <- new_ctx(); g <- fx_gwas_full()
  p <- methods_preprocessing(ctx, g, NULL)
  expect_true(grepl("quality-controlled", p))
  expect_true(grepl("0.9", p))      # INFO threshold read from config
  expect_true(grepl("9,000,000", p))
  expect_true(grepl("7,000,000", p))
  expect_true(grepl("GRCh37", p))
  expect_true("preprocessing" %in% ctx$ran_blocks)
})

test_that("LDSC methods + results emitted when the block ran, NULL otherwise", {
  ctx <- new_ctx()
  expect_true(grepl("LD Score Regression", methods_ldsc_h2(ctx, fx_gwas_full(), NULL)))
  expect_true(grepl("0.066", results_ldsc_h2(ctx, fx_gwas_full(), NULL)))
  # cross-ancestry GWAS has no LDSC output -> NULL (not "not performed")
  expect_null(methods_ldsc_h2(ctx, fx_gwas_null(), NULL))
  expect_null(results_ldsc_h2(ctx, fx_gwas_null(), NULL))
})

test_that("a not-run block (COJO requested but absent) yields NULL paragraphs", {
  ctx <- new_ctx(); g <- fx_gwas_full()   # cojo = NULL in fixture
  expect_null(methods_cojo(ctx, g, NULL))
  expect_null(results_cojo(ctx, g, NULL))
  expect_false("cojo" %in% ctx$ran_blocks)
})

test_that("a switched-off block (SumHer) yields NULL paragraphs", {
  ctx <- new_ctx()
  expect_null(methods_sumher(ctx, fx_gwas_full(), NULL))
  expect_null(results_sumher(ctx, fx_gwas_full(), NULL))
})

test_that("TWAS results report test count, FDR count and coloc subset", {
  ctx <- new_ctx()
  r <- results_twas(ctx, fx_gwas_full(), NULL)
  expect_true(grepl("tested 3 associations", r))
  expect_true(grepl("2 were significant", r))
  expect_true(grepl("1 also showed colocalisation", r))
  expect_true(grepl("P.FDR = 1.0e-07", r))
})

test_that("TWAS null-result produces an explicit none-significant sentence", {
  ctx <- new_ctx()
  r <- results_twas(ctx, fx_gwas_null(), NULL)
  expect_true(grepl("none were significant after FDR correction", r))
  expect_false(grepl("strongest association", r))
})

test_that("MAGMA results count FDR-significant genes", {
  ctx <- new_ctx()
  r <- results_magma_gene(ctx, fx_gwas_full(), NULL)
  expect_true(grepl("tested 2 associations", r))
  expect_true(grepl("1 was significant", r))
  expect_true(grepl("GENE1", r))
})

test_that("clump results report locus count and top variant", {
  ctx <- new_ctx()
  r <- results_clump(ctx, fx_gwas_full(), NULL)
  expect_true(grepl("2 independent genome-wide significant loci", r))
  expect_true(grepl("rs1", r))
  expect_true(grepl("GENEA", r))
})

test_that("present-but-empty SMR block is not described (results NULL)", {
  ctx <- new_ctx()
  expect_null(methods_smr_exp(ctx, fx_gwas_full(), NULL))   # results = NULL
  expect_null(results_smr_exp(ctx, fx_gwas_full(), NULL))
})

test_that("MAGMA version is flagged unpinned rather than fabricated", {
  ctx <- new_ctx()
  m <- methods_magma_gene(ctx, fx_gwas_full(), NULL)
  expect_true(grepl("version not pinned upstream", m))
  expect_true(any(grepl("MAGMA version is not recorded", ctx$warnings)))
})
