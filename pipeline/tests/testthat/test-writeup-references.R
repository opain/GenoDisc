# Coverage tests: every tool/dataset that can appear in a run has a references.bib entry.

wu_source()
pd <- wu_pipeline_dir()

test_that("every registry bibkey exists in references.bib", {
  missing <- writeup_missing_bibkeys(pd)
  expect_equal(length(missing), 0,
               info = paste("bibkeys missing from references.bib:", paste(missing, collapse = ", ")))
})

test_that("every block-enabling config flag is covered by the registry", {
  uncovered <- writeup_uncovered_flags(pd)
  expect_equal(length(uncovered), 0,
               info = paste("config flags with no registry entry:", paste(uncovered, collapse = ", ")))
})

test_that("a generated write-up cites nothing that is missing from references.bib", {
  res <- generate_writeup(fx_output(), fx_config(), fx_manifest(), pd)
  # no CITATION MISSING markers leaked into any output
  all_text <- paste(unlist(res$files), collapse = "\n")
  expect_false(grepl("CITATION MISSING", all_text, fixed = TRUE))
  # every used bibkey resolves to a parsed entry
  entries <- names(writeup_parse_bib(pd))
  expect_true(all(res$meta$bibkeys %in% entries))
})

test_that("the generated references.bib subset contains exactly the cited keys", {
  res <- generate_writeup(fx_output(), fx_config(), fx_manifest(), pd)
  bib <- res$files[["references.bib"]]
  for (k in res$meta$bibkeys) expect_true(grepl(paste0("{", k, ","), bib, fixed = TRUE))
})
