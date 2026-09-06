#!/usr/bin/Rscript
# Clean one TWAS-GSEA-fast .competitive.txt for one (weight × gmt) pair
# into a tidy CSV. Non-directional only — TWAS-GSEA emits P as a
# one-sided right-tail p already, so we just carry it through.
#
# FDR is not applied here: read_pathway_twas_gsea() in
# package_results_functions.R pools across every weight × every gmt for the
# same primary GWAS and applies BH once at bundle-read time (per user spec).

suppressMessages(library("optparse"))
suppressMessages(library("data.table"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--competitive", action = "store", type = "character",
              help = "Path to TWAS-GSEA .competitive.txt"),
  make_option("--out_csv", action = "store", type = "character",
              help = "Path to write cleaned CSV")
)))

d <- fread(opt$competitive, sep = " ")
if (nrow(d) == 0) {
  # Emit a header-only CSV so downstream fan-in doesn't crash on an empty
  # panel × gmt pair (e.g. panel has 0 usable genes in that gmt).
  fwrite(data.table(Name = character(), N = integer(),
                    Estimate = numeric(), SE = numeric(),
                    Z = numeric(), P = numeric()),
         opt$out_csv)
  cat("format_twas_gsea_pathway_results: empty input -> header-only CSV\n")
  quit(status = 0)
}

# Column name defensiveness: TWAS-GSEA-fast has had two schemas (with and
# without N_Mem). Prefer N_Mem_Avail (the count of genes actually usable in
# this panel).
N_col <- if ("N_Mem_Avail" %in% names(d)) "N_Mem_Avail" else "N"
out <- data.table(
  Name     = as.character(d$GeneSet),
  N        = suppressWarnings(as.integer(d[[N_col]])),
  Estimate = suppressWarnings(as.numeric(d$Estimate)),
  SE       = suppressWarnings(as.numeric(d$SE)),
  Z        = suppressWarnings(as.numeric(d$T)),
  P        = suppressWarnings(as.numeric(d$P))
)
fwrite(out, opt$out_csv)
cat(sprintf("format_twas_gsea_pathway_results: %d rows -> %s\n",
            nrow(out), opt$out_csv))
