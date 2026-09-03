#!/usr/bin/Rscript
# Pre-download BioMart gene annotation data for use by downstream scripts.
# Saves:
#   - biomart_genes_grch37.tsv (ensembl_gene_id, external_gene_name, external_synonym, chromosome_name, start_position, end_position)
#   - gene_locations.tsv (chromosome_name, start_position, end_position, strand, external_gene_name, gene_biotype; chrs 1-22, X, Y)

library(biomaRt)
library(optparse)

option_list = list(
  make_option("--resdir", type="character", default="resources")
)
option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

outdir <- paste0(opt$resdir, "/data/biomart")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# GRCh37
ensembl37 <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
biomartCacheClear()
genes37 <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'external_synonym',
                  'chromosome_name', 'start_position', 'end_position'),
  mart = ensembl37
)
write.table(genes37, file.path(outdir, "biomart_genes_grch37.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Gene locations (separate query so the existing output above is untouched)
gene_locations <- getBM(
  attributes = c('chromosome_name', 'start_position', 'end_position', 'strand',
                  'external_gene_name', 'gene_biotype'),
  mart = ensembl37
)
gene_locations <- gene_locations[gene_locations$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
write.table(gene_locations, file.path(outdir, "gene_locations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
