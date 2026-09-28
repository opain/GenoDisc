#!/usr/bin/env Rscript
# format_magma_drugtargetor_gls_results.R
#
# DRUGSETS-style GLS ATC drug-class enrichment on MAGMA per-drug results.
# Regresses the per-drug competitive statistic (BETA/SE from MAGMA's gene-set
# analysis) on an ATC-class indicator + size covariates, using the drug-drug
# correlation Sigma (from compute_magma_drugcorr.R) as the GLS weight. This is
# the MAGMA analogue of format_twas_gsea_drugtargetor_gls_results.R and follows
# DRUGSETS (Bell, de Leeuw, Posthuma; github.com/nybell/drugsets/compute_lnreg.R).
#
# The test is competitive and NON-DIRECTIONAL (MAGMA magnitude): it asks whether
# in-class drugs are, on average, more trait-enriched than out-of-class drugs.
#
# Usage:
#   format_magma_drugtargetor_gls_results.R <gsa.out> <drugcorr.rds> <atc_labels.txt> <out_dir> [floor] [nsize]

suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 4) stop("usage: <gsa.out> <drugcorr.rds> <atc_labels.txt> <out_dir> [floor] [nsize]")
gsa.file <- a[1]; corr.file <- a[2]; atc.file <- a[3]; out.dir <- a[4]
floor <- ifelse(length(a) >= 5, as.numeric(a[5]), 0.1)
nsize <- ifelse(length(a) >= 6, as.integer(a[6]), 5L)
level_width <- c(l2 = 3L, l3 = 4L, l4 = 5L)

# Regularised (pseudo)inverse: drop eigenvalues < floor (DRUGSETS prune 0.1).
reg_inv <- function(S, floor) {
  e <- eigen(S, symmetric = TRUE); k <- e$values >= floor
  list(inv = e$vectors[, k, drop = FALSE] %*% (t(e$vectors[, k, drop = FALSE]) / e$values[k]), rank = sum(k))
}
# GLS of y on [1, s, size, log(size)]; one-sided upper-tail test on the class coef.
# df-corrected: Sinv is a rank-m pseudo-inverse, so the residual quadratic form
# r'Sinv r ~ sigma^2 * chi^2(m - K). Divide by (m - K) and use df = m - K (NOT N - K),
# else sigma^2 is under-estimated and the test is anti-conservative.
gls_group <- function(y, s, size, Sinv, m) {
  X <- cbind(1, s, size, log(size)); K <- ncol(X); df <- m - K
  W <- solve(t(X) %*% Sinv %*% X)
  B <- W %*% t(X) %*% Sinv %*% y
  resid <- y - X %*% B
  sigma <- as.numeric((t(resid) %*% Sinv %*% resid) / df)
  se <- sqrt(sigma * W[2, 2]); tval <- B[2] / se
  list(B = B[2], SE = se, T = tval, P = pt(tval, df, lower.tail = FALSE))  # one-sided (enrichment)
}

## --- per-drug MAGMA statistic ------------------------------------------------
g <- as.data.table(read.table(gsa.file, header = TRUE, comment.char = "#", stringsAsFactors = FALSE))
g <- g[is.finite(BETA) & is.finite(SE) & SE > 0 & NGENES > 0]
g[, stat := BETA / SE]

## --- align to Sigma ----------------------------------------------------------
Sig <- readRDS(corr.file)
common <- intersect(g$FULL_NAME, rownames(Sig))
g <- g[match(common, FULL_NAME)]
Sig <- Sig[common, common]
ri <- reg_inv(Sig, floor); Sinv <- ri$inv; m_rank <- ri$rank   # m = retained-eigenvalue rank (df = m - K)
y <- g$stat; size <- as.numeric(g$NGENES)

## --- drug -> ATC 7-char codes (multi-code exploded) --------------------------
codes7 <- strsplit(sub("^ATC:([^|]+)\\|.*", "\\1", g$FULL_NAME), ",")

## --- ATC labels --------------------------------------------------------------
atc <- fread(atc.file, sep = "!", header = FALSE, col.names = c("Code", "Name"))
atc[, Name := tolower(Name)]

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
for (lv in names(level_width)) {
  k <- level_width[[lv]]
  drug_cls <- lapply(codes7, function(cc) unique(substr(cc, 1, k)))
  classes  <- sort(unique(unlist(drug_cls)))
  classes  <- classes[nchar(classes) == k]
  rows <- list()
  for (cl in classes) {
    s <- as.numeric(vapply(drug_cls, function(x) cl %in% x, logical(1)))
    n1 <- sum(s == 1)
    if (n1 < nsize) next
    r <- gls_group(y, s, size, Sinv, m_rank)
    rows[[length(rows) + 1]] <- data.table(Code = cl, N_Drugs = n1, Estimate = r$B,
                                            SE = r$SE, T = r$T, P = r$P)
  }
  res <- rbindlist(rows)
  if (nrow(res)) {
    res[, FDR := p.adjust(P, "fdr")]
    res[, Direction := ifelse(Estimate > 0, "Enriched", "Depleted")]  # competitive, non-directional
    res <- merge(res, atc[nchar(Code) == k], by = "Code", all.x = TRUE)
    setcolorder(res, c("Code", "Name", "N_Drugs", "Estimate", "SE", "T", "P", "FDR", "Direction"))
    res <- res[order(P)]
  }
  fwrite(res, file.path(out.dir, paste0("magma_drug_targetor_gls_", lv, "_res.csv")))
  cat("Wrote", lv, ":", nrow(res), "ATC classes\n")
}
