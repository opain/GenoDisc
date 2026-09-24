#!/usr/bin/Rscript

# Build gene x ATC-class gene sets for TWAS-GSEA, so ATC-class enrichment can be
# tested directly at the gene level ("option C") instead of the two-stage
# per-drug + Wilcoxon aggregation. Each ATC class becomes one gene set;
# TWAS-GSEA-fast.R then tests it as a gene-level GLS/REML mixed model, with drugs
# no longer treated as independent observations (removing within-class
# pseudo-replication).
#
# Emits, per ATC level:
#   * a SIGNED .prop  (directional run, --prop_file --directional T)
#   * an unsigned .gmt (non-directional run, --gmt_file --directional F)
#
# Built from the RAW DrugTargetor database (not the per-drug .prop, which
# collapses each multi-code drug into a single concatenated column). A single
# drug/substance can carry several comma-separated ATC codes; those are exploded
# so the drug contributes to every ATC class it belongs to.

suppressMessages(library(optparse))

option_list <- list(
  make_option("--resdir", type = "character", default = "resources",
              help = "Resources dir (reads data/drug_targetor/wholedatabase_for_targetor and data/magma/NCBI37.3.gene.loc)"),
  make_option("--outdir", type = "character", default = NA,
              help = "Output dir for the .prop/.gmt files [default: {resdir}/data/drug_targetor]"),
  make_option("--pipeline_dir", type = "character", default = NA,
              help = "Path to the pipeline directory"),
  make_option("--agg", type = "character", default = "mean",
              help = "Cross-drug within-class sign aggregation for the .prop: 'mean' (continuous consensus in [-1,1], default) or 'netsign' (discrete sign(sum))"),
  make_option("--levels", type = "character", default = "l2,l3,l4",
              help = "Comma-separated ATC levels to emit: l2 (substr 1-3), l3 (1-4), l4 (1-5)")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!(opt$agg %in% c("mean", "netsign"))) stop("--agg must be 'mean' or 'netsign'")

suppressMessages(library(data.table))

outdir <- if (is.na(opt$outdir)) file.path(opt$resdir, "data", "drug_targetor") else opt$outdir
level_width <- c(l2 = 3L, l3 = 4L, l4 = 5L)   # ATC hierarchy: L2=3 chars, L3=4, L4=5
levels <- trimws(strsplit(opt$levels, ",", fixed = TRUE)[[1]])
if (!all(levels %in% names(level_width))) stop("--levels must be a subset of l2,l3,l4")

# Gene universe = HGNC symbols (V6) of the MAGMA gene.loc, matching the existing
# .prop / .gmt / TWAS ID space (TWAS-GSEA matches on symbol via --use_alt_id ID).
gene_loc <- fread(file.path(opt$resdir, "data", "magma", "NCBI37.3.gene.loc"), header = FALSE)
universe <- unique(gene_loc$V6)

# Raw assertions: atc, original_name, gene, activity_type, source
db <- fread(file.path(opt$resdir, "data", "drug_targetor", "wholedatabase_for_targetor"),
            sep = "\t", header = TRUE)
db <- db[gene %in% universe]

# Signed membership (for .prop), same mapping as format_drug_targetor_for_twas_gsea.R.
db[, sgn := fifelse(activity_type %in% c("DECREASED_EXPRESSION", "NEGATIVE_RESPONSE", "OPPOSITE_RESPONSE"), -1L,
             fifelse(activity_type %in% c("INCREASED_EXPRESSION", "POSITIVE_RESPONSE"), 1L, 0L))]
# One sign per (drug, gene) with -1 precedence, mirroring the existing within-drug tie-break.
dg <- db[sgn != 0L, .(sgn = if (any(sgn == -1L)) -1L else 1L), by = .(drug = atc, gene)]
# Unsigned membership (for .gmt): any targeting assertion, mirroring format_drug_targetor.R.
dg_all <- unique(db[, .(drug = atc, gene)])

# Explode the composite ATC field into individual 7-char codes per drug.
udrugs <- unique(db$atc)
codes_list <- strsplit(sub("^ATC:([^|]+)\\|.*", "\\1", udrugs), ",", fixed = TRUE)
drug_codes <- data.table(drug = rep(udrugs, lengths(codes_list)), code = unlist(codes_list))

for (lv in levels) {
  k <- level_width[[lv]]
  # drug -> distinct classes at this level (dedupe so a drug counts once per class)
  dc <- unique(drug_codes[, .(drug, class = substr(code, 1, k))])

  # --- signed .prop ---
  m <- merge(dc, dg, by = "drug", allow.cartesian = TRUE)   # drug, class, gene, sgn
  agg <- if (opt$agg == "mean") {
    m[, .(val = mean(sgn)), by = .(class, gene)]             # consensus in [-1,1]
  } else {
    m[, .(val = sign(sum(sgn))), by = .(class, gene)]        # discrete -1/0/1
  }
  agg[, setname := paste0("ATC_", class)]
  wide <- dcast(agg, gene ~ setname, value.var = "val", fill = 0)
  out <- merge(data.table(ID = universe), wide, by.x = "ID", by.y = "gene", all.x = TRUE)
  for (j in setdiff(names(out), "ID")) set(out, which(is.na(out[[j]])), j, 0)
  prop_f <- file.path(outdir, paste0("wholedatabase_for_targetor_atc_", lv, ".prop"))
  fwrite(out, prop_f, sep = ",", quote = TRUE, na = "NA")

  # --- unsigned .gmt (setname \t <desc> \t gene1 \t gene2 ...) ---
  mg <- unique(merge(dc, dg_all, by = "drug", allow.cartesian = TRUE)[, .(class, gene)])
  gmt_f <- file.path(outdir, paste0("wholedatabase_for_targetor_atc_", lv, ".gmt"))
  con <- file(gmt_f, "w")
  for (cl in sort(unique(mg$class))) {
    genes <- mg[class == cl, gene]
    writeLines(paste(c(paste0("ATC_", cl), " ", genes), collapse = "\t"), con)
  }
  close(con)

  cat(sprintf("%s: %d genes x %d prop-classes / %d gmt-classes -> %s , %s\n",
              lv, nrow(out), ncol(out) - 1L, length(unique(mg$class)), prop_f, gmt_f))
}
