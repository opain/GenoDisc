#!/usr/bin/env Rscript
# Build targets.csv from the existing 63 shipped codes intersected with
# sumstat_repo.csv EUR rows, resolving source paths in
# /scratch/prj/gwas_sumstats/cleaned/. Rebuild-only for the currently
# shipping panel; growing the panel is out of scope.
suppressMessages(library(data.table))

SGWAS <- "/scratch/prj/neurohackpain/GenoDisc/pipeline_resources/data/secondary_gwas"
CLEANED <- "/scratch/prj/gwas_sumstats/cleaned"
MISC <- "/scratch/prj/neurohackpain/GenoDisc/repo/current/pipeline/misc/secondary_gwas"

shipped <- sub("\\.sumstats\\.gz$", "",
               list.files(SGWAS, pattern = "\\.sumstats\\.gz$"))
stopifnot(length(shipped) > 0)

repo <- fread(file.path(SGWAS, "sumstat_repo.csv"))
# Handle BOM-mangled first column name if present.
setnames(repo, names(repo)[1], "code")

meta <- repo[code %in% shipped]
missing_from_repo <- setdiff(shipped, meta$code)

# Strict EUR only. Strip trailing metadata parentheticals like
# "European (N=123,665)" or "EUR (UKB)" before matching, so annotated
# EUR-only entries aren't spuriously excluded. Mixed strings like
# "European and Asian ..." still fail the exact match.
meta[, ancestry_norm := tolower(trimws(gsub("\\s*\\([^)]*\\)", "", ancestry)))]
excl <- meta[!ancestry_norm %in% c("eur", "european"),
             .(code, ancestry, status = "excluded_non_eur")]
keep <- meta[ancestry_norm %in% c("eur", "european")]

# Resolve source path in cleaned/.
keep[, source_path := file.path(CLEANED, paste0(code, ".gz"))]
no_src <- keep[!file.exists(source_path),
               .(code, ancestry, source_path,
                 status = "no_source_file")]
keep <- keep[file.exists(source_path)]

# N fallback for --n. Prefer sample_size_discovery; else n_cases+n_controls.
num <- function(x) suppressWarnings(as.numeric(gsub("[,\\s]", "", x, perl = TRUE)))
keep[, n_disc := num(sample_size_discovery)]
keep[, n_cas  := num(n_cases)]
keep[, n_con  := num(n_controls)]
keep[, n := fifelse(!is.na(n_disc), n_disc,
             fifelse(!is.na(n_cas) & !is.na(n_con), n_cas + n_con, NA_real_))]

targets <- keep[, .(
  code,
  source_path,
  n,
  ancestry,
  phenotype_type,
  category,
  trait_label,
  short_label,
  n_cases   = n_cas,
  n_controls = n_con,
  sample_size_discovery = n_disc,
  pmid
)]
setorder(targets, code)

# Sanity: at least one shipped code should survive.
stopifnot(nrow(targets) > 0)

dir.create(MISC, showWarnings = FALSE, recursive = TRUE)
fwrite(targets, file.path(MISC, "targets.csv"))

# Also emit an exclusion prelog so downstream reconciliation can include
# codes that never entered the SLURM array.
parts <- list()
if (length(missing_from_repo) > 0)
  parts[[length(parts) + 1]] <- data.table(
    code = missing_from_repo, ancestry = NA_character_,
    source_path = NA_character_, status = "missing_from_repo")
if (nrow(excl) > 0)
  parts[[length(parts) + 1]] <- excl[, .(code, ancestry,
                                         source_path = NA_character_, status)]
if (nrow(no_src) > 0)
  parts[[length(parts) + 1]] <- no_src[, .(code, ancestry, source_path, status)]
prelog <- if (length(parts)) rbindlist(parts, fill = TRUE) else
  data.table(code = character(), ancestry = character(),
             source_path = character(), status = character())
fwrite(prelog, file.path(MISC, "targets_excluded.csv"))

cat(sprintf("shipped_currently=%d, in_repo=%d, kept=%d, excluded=%d\n",
            length(shipped), nrow(meta), nrow(targets), nrow(prelog)))
