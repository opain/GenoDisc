#!/usr/bin/Rscript
# Pre-download BioMart gene annotation data for use by downstream scripts.
# Saves two files:
#   - biomart_genes_grch37.tsv (ensembl_gene_id, external_gene_name, external_synonym, chromosome_name, start_position, end_position)
#   - biomart_genes_grch38.tsv (ensembl_gene_id, external_gene_name)

library(biomaRt)
library(optparse)

option_list = list(
  make_option("--resdir", type="character", default="resources")
)
opt = parse_args(OptionParser(option_list=option_list))

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
