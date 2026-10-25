# Snapshot test on a full fixture run, so template drift is visible in diffs.
# Golden files live in tests/testthat/_snaps/writeup/. To update after an
# intentional change, delete the relevant golden file and re-run (it is
# recreated), then review the diff before committing.

wu_source()
pd <- wu_pipeline_dir()
snap_dir <- file.path(pd, "tests", "testthat", "_snaps", "writeup")

.norm <- function(s) sub("\n+$", "", paste(s, collapse = "\n"))

test_that("full-run write-up matches the committed snapshot", {
  res <- generate_writeup(fx_output(), fx_config(), fx_manifest(), pd)
  dir.create(snap_dir, recursive = TRUE, showWarnings = FALSE)
  rels <- c("methods.md", "TRAIT1/results.md", "TRAIT2/results.md",
            "references.bib", "references.ris")
  golden_of <- function(rel) file.path(snap_dir, gsub("/", "__", rel))
  missing <- rels[!file.exists(vapply(rels, golden_of, character(1)))]
  if (length(missing) > 0) {
    for (rel in missing) writeLines(res$files[[rel]], golden_of(rel))
    skip(paste("created", length(missing), "snapshot(s) - re-run to compare"))
  }
  for (rel in rels) {
    expect_identical(.norm(res$files[[rel]]), .norm(readLines(golden_of(rel), warn = FALSE)),
                     info = paste("snapshot drift in", rel,
                                  "- review and update", golden_of(rel), "if intended"))
  }
})
