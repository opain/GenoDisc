# Unit tests for the write-up number / string formatters.

wu_source()

test_that("fmt_num uses significant figures and is NA-safe", {
  expect_equal(fmt_num(0.066), "0.066")
  expect_equal(fmt_num(0.9), "0.9")
  expect_equal(fmt_num(1.2345, sig = 3), "1.23")
  expect_equal(fmt_num(1.01), "1.01")
  expect_equal(fmt_num(NA), "NA")
  expect_equal(fmt_num(NULL), "NA")
})

test_that("fmt_int adds thousands separators and is NA-safe", {
  expect_equal(fmt_int(1234), "1,234")
  expect_equal(fmt_int(9000000), "9,000,000")
  expect_equal(fmt_int(NA), "NA")
})

test_that("fmt_p renders scientific for small P and fixed otherwise", {
  expect_equal(fmt_p(1e-9), "1.0e-09")
  expect_equal(fmt_p(1e-12), "1.0e-12")
  expect_equal(fmt_p(0.03), "0.030")
  expect_equal(fmt_p(NA), "NA")
  # plain style never collapses a tiny P to 0
  expect_true(grepl("e", fmt_p(1e-20, writeup_fmt_opts(p_style = "plain"))))
})

test_that("pluralise / was_were / oxford agree grammatically", {
  expect_equal(pluralise("gene", 1), "gene")
  expect_equal(pluralise("gene", 2), "genes")
  expect_equal(pluralise("locus", 2, "loci"), "loci")
  expect_equal(was_were(1), "was")
  expect_equal(was_were(2), "were")
  expect_equal(oxford(c("a")), "a")
  expect_equal(oxford(c("a", "b")), "a and b")
  expect_equal(oxford(c("a", "b", "c")), "a, b and c")
  expect_equal(oxford(character(0)), "")
})

test_that("count_below counts FDR-significant rows, NA-safe", {
  df <- data.frame(P.FDR = c(0.01, 0.2, NA, 0.049))
  expect_equal(count_below(df, "P.FDR", 0.05), 2L)
  expect_equal(count_below(NULL, "P.FDR"), 0L)
  expect_equal(count_below(data.frame(x = 1), "P.FDR"), 0L)
  expect_equal(count_below(df[0, , drop = FALSE], "P.FDR"), 0L)
})

test_that("pluck safely traverses absent branches", {
  x <- list(a = list(b = 1))
  expect_equal(pluck(x, "a", "b"), 1)
  expect_null(pluck(x, "a", "z"))
  expect_null(pluck(x, "nope", "b"))
  expect_null(pluck(NULL, "a"))
})
