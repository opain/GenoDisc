#!/usr/bin/Rscript
# Clean the MAGMA .gsa.out for one pathway gmt into a tidy CSV consumed by
# the bundle reader. FDR is NOT computed here — read_pathway_magma() in
# package_results_functions.R pools across every gmt for the same primary
# GWAS and applies BH once at bundle-read time (per the user's spec).

suppressMessages(library("optparse"))
suppressMessages(library("data.table"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gsa_out", action = "store", type = "character",
              help = "Path to magma_pathway_<gmt>.gsa.out"),
  make_option("--out_csv", action = "store", type = "character",
              help = "Path to write magma_pathway_<gmt>.clean.csv")
)))

# .gsa.out is whitespace-separated with a `#` comment header block; keep
# only the data rows via `grep -v '^#'` — same idiom used by
# format_magma_gsea_results.R.
d <- fread(cmd = paste0("grep -v '^#' ", opt$gsa_out))

# MAGMA emits VARIABLE (or FULL_NAME) + TYPE + NGENES + BETA + BETA_STD +
# SE + P. For pathway sets, the set name is the column of interest. Fall
# back gracefully across MAGMA column-name variants.
name_col <- if ("FULL_NAME" %in% names(d)) {
  "FULL_NAME"
} else if ("VARIABLE" %in% names(d)) {
  "VARIABLE"
} else {
  names(d)[1]
}
d[, Name := as.character(get(name_col))]

keep <- c("Name", "NGENES", "BETA", "SE", "P")
missing <- setdiff(keep, names(d))
if (length(missing) > 0) stop("Missing MAGMA output columns: ",
                              paste(missing, collapse = ","))

fwrite(d[, ..keep], opt$out_csv)
cat(sprintf("format_magma_pathway_results: wrote %d rows -> %s\n",
            nrow(d), opt$out_csv))
