#!/usr/bin/Rscript
# Convert a MSigDB-style .gmt from Entrez gene IDs to HGNC symbols.
#
# TWAS-GSEA-fast.R matches gene set members to the TWAS gene column by
# symbol, so any Entrez-only .gmt has to be translated first. MAGMA takes
# Entrez directly and doesn't need this step.
#
# Mapping comes from NCBI37.3.gene.loc — cols: entrez chr start end strand
# symbol — the same file MAGMA uses. Genes whose Entrez ID has no matching
# row are dropped; drop counts are logged per set.
#
# Usage:
#   Rscript convert_gmt_entrez_to_symbol.R <in.gmt> <out.gmt> <gene.loc>

suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
  stop("usage: convert_gmt_entrez_to_symbol.R <in.gmt> <out.gmt> <gene.loc>")
in_gmt  <- args[1]
out_gmt <- args[2]
loc_fp  <- args[3]

loc <- fread(loc_fp, header = FALSE,
             col.names = c("entrez", "chr", "start", "end", "strand", "symbol"))
loc[, entrez := as.character(entrez)]
map <- setNames(loc$symbol, loc$entrez)

# .gmt is tab-separated but rows have variable numbers of fields, so read
# with readLines and split. Standard schema: <set_name> <url> <gene1>
# <gene2> ...
lines <- readLines(in_gmt, warn = FALSE)
n_in  <- length(lines)

total_in  <- 0L
total_out <- 0L
out_lines <- character(n_in)
for (i in seq_along(lines)) {
  parts <- strsplit(lines[i], "\t", fixed = TRUE)[[1L]]
  if (length(parts) < 3L) {
    out_lines[i] <- lines[i]
    next
  }
  head  <- parts[1:2]
  genes <- parts[-(1:2)]
  total_in <- total_in + length(genes)
  syms <- map[genes]
  syms <- syms[!is.na(syms) & nzchar(syms)]
  total_out <- total_out + length(syms)
  out_lines[i] <- paste(c(head, syms), collapse = "\t")
}

writeLines(out_lines, out_gmt)

cat(sprintf(
  "convert_gmt_entrez_to_symbol: %s -> %s\n  %d sets, %d/%d genes retained (%.1f%% mapped).\n",
  in_gmt, out_gmt, n_in, total_out, total_in,
  100 * total_out / max(total_in, 1L)))
