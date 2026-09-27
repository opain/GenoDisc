#!/usr/bin/env Rscript
# compute_magma_drugcorr.R
#
# Build the drug-drug correlation matrix Sigma for the MAGMA DRUGSETS-style GLS
# ATC enrichment, faithfully reproducing DRUGSETS' compute_corrs.r
# (Bell, de Leeuw, Posthuma; github.com/nybell/drugsets) from the MAGMA
# gene-level results (.genes.raw), which store the gene-gene correlations.
#
# Deviation from DRUGSETS (per project decision): the per-drug MAGMA run is the
# existing competitive test WITHOUT druggable-background conditioning, so Sigma
# residualises drug membership on the gene-level covariates + intercept only
# (no druggable indicator).
#
# Usage:
#   compute_magma_drugcorr.R <genes.raw> <drug.gmt> <gsa.out> <out.rds>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: compute_magma_drugcorr.R <genes.raw> <drug.gmt> <gsa.out> <out.rds>")
raw.file <- args[1]; set.file <- args[2]; res.file <- args[3]; out.file <- args[4]
prune.thresh <- 0.1

## --- gene sets present in the MAGMA results (drugs to include) ---------------
gsa.res  <- read.table(res.file, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
sets     <- strsplit(scan(set.file, what = "", sep = "\n", quiet = TRUE), "\t")
set.names.all <- vapply(sets, `[`, character(1), 1)
sets.use <- gsa.res$FULL_NAME                       # no "druggable" set (see header)
keep     <- set.names.all %in% sets.use
set.genes <- lapply(sets[keep], function(x) x[-1][nzchar(x[-1])])  # drop name (+ any blank field)
set.names <- set.names.all[keep]

## --- MAGMA .genes.raw: info fields + trailing gene-gene correlations ---------
# fields: gene chr start stop nsnps nparam nsamp mac zstat  <corrs with preceding genes...>
raw.data    <- strsplit(scan(raw.file, what = "", comment.char = "#", sep = "\n", quiet = TRUE), " ")
info.length <- length(raw.data[[1]])                # 9 (first gene has no correlations)
raw.info    <- data.frame(matrix(vapply(raw.data, function(x) x[seq_len(info.length)], character(info.length)),
                                 ncol = info.length, byrow = TRUE),
                          stringsAsFactors = FALSE)[, -(3:4)]
names(raw.info) <- c("gene", "chr", "nsnps", "nparam", "nsamp", "mac", "zstat")
for (cc in c("nsnps", "nparam", "nsamp", "mac", "zstat")) raw.info[[cc]] <- as.numeric(raw.info[[cc]])
raw.corrs <- lapply(raw.data, function(x) if (length(x) > info.length) x[-seq_len(info.length)] else character(0))

## --- membership matrix (genes x drugs), Entrez match -------------------------
sets <- matrix(0, nrow = nrow(raw.info), ncol = length(set.names))
for (i in seq_along(set.names)) sets[raw.info$gene %in% set.genes[[i]], i] <- 1
keep.drug <- colSums(sets) > 0                      # drugs with >=1 gene in the universe
sets <- sets[, keep.drug, drop = FALSE]
set.names <- set.names[keep.drug]

## --- per-chromosome whitening projection from gene-gene correlations ---------
chr.names  <- unique(raw.info$chr)
projection <- list(); proj.index <- list(); no.proj <- 0L
for (chr in chr.names) {
  idx      <- which(raw.info$chr == chr)
  curr.raw <- raw.corrs[idx]
  curr.size <- length(curr.raw)
  curr.corr <- matrix(0, curr.size, curr.size)
  if (curr.size >= 2) for (i in 2:curr.size) {
    len <- length(curr.raw[[i]])
    if (len > 0) curr.corr[i, seq_len(len) + (i - len - 1)] <- as.numeric(curr.raw[[i]])
  }
  curr.corr <- curr.corr + t(curr.corr); diag(curr.corr) <- 1
  eig <- eigen(curr.corr, symmetric = TRUE)
  use <- which(eig$values >= prune.thresh)
  projection[[chr]] <- eig$vectors[, use, drop = FALSE] %*% diag(1 / sqrt(eig$values[use]), length(use))
  proj.index[[chr]] <- seq_along(use) + no.proj
  no.proj <- no.proj + length(use)
}
project.data <- function(data) {
  out <- matrix(NA_real_, nrow = no.proj, ncol = ncol(data))
  for (chr in chr.names) out[proj.index[[chr]], ] <- t(projection[[chr]]) %*% data[raw.info$chr == chr, , drop = FALSE]
  out
}

## --- residualise membership on gene covariates (+ intercept), whitened -------
residualize <- cbind(as.matrix(raw.info[, c("nsnps", "nparam", "nsamp")]), 1 / raw.info$mac)
residualize <- residualize[, apply(residualize, 2, var) > 0, drop = FALSE]
residualize <- cbind(1, scale(cbind(residualize, log(residualize))))
residualize <- project.data(residualize)
sets        <- project.data(sets)

covar   <- residualize[, 1, drop = FALSE]           # intercept only (no druggable)
ctc.inv <- solve(t(covar) %*% covar)
cts     <- t(covar) %*% sets
det     <- colSums(sets^2) - colSums(cts * (ctc.inv %*% cts))
V       <- sweep(sets - covar %*% ctc.inv %*% cts, 2, det, FUN = "/")

set.corrs <- cov2cor(t(V) %*% V)                    # drug x drug
dimnames(set.corrs) <- list(set.names, set.names)

saveRDS(set.corrs, out.file)
cat("Wrote", out.file, ":", nrow(set.corrs), "x", ncol(set.corrs), "drug-drug correlation matrix",
    "(", no.proj, "whitened dims from", length(chr.names), "chromosomes )\n")
