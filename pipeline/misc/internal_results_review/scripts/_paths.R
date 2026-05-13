# Shared path resolution for internal_results_review scripts.
# Sourced at the top of every other script in this folder.
#
# Resolution rules:
#   rds_path     <- env GENODISC_RDS    (default /work/results/results_package.rds)
#   out_base     <- env GENODISC_OUTDIR (default /work/results/internal_results_review)
#   scripts_dir  <- directory of the *calling* script (auto-detected)
#
# Override by exporting env vars before running, e.g.
#   GENODISC_RDS=/path/to/other.rds GENODISC_OUTDIR=/tmp/run1 \
#     Rscript scripts/extract_twas_gsea_results.R

get_caller_script_dir <- function() {
  # When run via Rscript: --file=/absolute/path/to/script.R is in commandArgs.
  # When sourced interactively: fall back to getwd().
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa) > 0) return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  # When sourced from another script: sys.frame(1)$ofile
  for(i in seq_along(sys.frames())){
    fname <- sys.frames()[[i]]$ofile
    if(!is.null(fname)) return(dirname(normalizePath(fname)))
  }
  getwd()
}

rds_path    <- Sys.getenv("GENODISC_RDS",    unset = "/work/results/results_package.rds")
out_base    <- Sys.getenv("GENODISC_OUTDIR", unset = "/work/results/internal_results_review")
scripts_dir <- get_caller_script_dir()

tables_dir  <- file.path(out_base, "tables")
figures_dir <- file.path(out_base, "figures")
report_dir  <- file.path(out_base, "report")

dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir,  showWarnings = FALSE, recursive = TRUE)

if (interactive() || Sys.getenv("GENODISC_VERBOSE", unset = "0") == "1") {
  message("[_paths.R] rds_path = ", rds_path)
  message("[_paths.R] out_base = ", out_base)
}
