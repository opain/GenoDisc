#!/usr/bin/env Rscript
# format_magma_drugtargetor_vif_results.R
#
# CAMERA-style VIF-OLS ATC drug-class enrichment on MAGMA per-drug results -- the
# recommended drug-level ATC test (see docs/atc_enrichment_estimator_validation.Rmd).
# Regresses the per-drug competitive statistic (BETA/SE) on an ATC-class indicator +
# size covariates by OLS, inflating the class-coefficient SE for within-class drug-drug
# correlation: SE = SE_OLS * sqrt(1 + (n1-1) * rho_bar), rho_bar = mean in-class entry of
# the drug-drug correlation Sigma (from compute_magma_drugcorr.R). Unlike the GLS this
# does not invert Sigma, avoiding the pseudo-inverse out-of-class leakage.
#
# The test is competitive and NON-DIRECTIONAL (MAGMA magnitude), one-sided (enrichment):
# are in-class drugs, on average, more trait-enriched than out-of-class drugs?
#
# Usage:
#   format_magma_drugtargetor_vif_results.R <gsa.out> <drugcorr.rds> <atc_labels.txt> <out_dir> [nsize]

suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 4) stop("usage: <gsa.out> <drugcorr.rds> <atc_labels.txt> <out_dir> [nsize]")
gsa.file <- a[1]; corr.file <- a[2]; atc.file <- a[3]; out.dir <- a[4]
nsize <- ifelse(length(a) >= 5, as.integer(a[5]), 5L)
level_width <- c(l2 = 3L, l3 = 4L, l4 = 5L)

# VIF-OLS of y on [1, s, size, log(size)]; one-sided upper-tail test on the class coef.
vif_group <- function(y, s, size, rho_bar, N) {
  X <- cbind(1, s, size, log(size)); K <- ncol(X)
  XtXi <- solve(crossprod(X)); B <- XtXi %*% crossprod(X, y)
  resid <- y - X %*% B
  sigma <- as.numeric(crossprod(resid)) / (N - K)
  vif <- 1 + (sum(s) - 1) * rho_bar
  se <- sqrt(sigma * XtXi[2, 2] * vif); tval <- B[2] / se; df <- N - K
  list(B = B[2], SE = se, T = tval, P = pt(tval, df, lower.tail = FALSE))  # one-sided (enrichment)
}

## --- per-drug MAGMA statistic ------------------------------------------------
g <- as.data.table(read.table(gsa.file, header = TRUE, comment.char = "#", stringsAsFactors = FALSE))
g <- g[is.finite(BETA) & is.finite(SE) & SE > 0 & NGENES > 0]
g[, stat := BETA / SE]

## --- align to Sigma by full drug name ----------------------------------------
Sig <- readRDS(corr.file)
common <- intersect(g$FULL_NAME, rownames(Sig))
g <- g[match(common, FULL_NAME)]
Sig <- Sig[common, common, drop = FALSE]
y <- g$stat; size <- as.numeric(g$NGENES); N <- length(y)

## --- drug -> ATC codes (multi-code exploded) ---------------------------------
codes_all <- strsplit(sub("^ATC:([^|]+)\\|.*", "\\1", g$FULL_NAME), ",")

## --- ATC labels --------------------------------------------------------------
atc <- fread(atc.file, sep = "!", header = FALSE, col.names = c("Code", "Name"))
atc[, Name := tolower(Name)]

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
for (lv in names(level_width)) {
  k <- level_width[[lv]]
  drug_cls <- lapply(codes_all, function(cc) unique(substr(cc, 1, k)))
  classes  <- sort(unique(unlist(drug_cls)))
  classes  <- classes[nchar(classes) == k]
  rows <- list()
  for (cl in classes) {
    s <- as.numeric(vapply(drug_cls, function(x) cl %in% x, logical(1)))
    n1 <- sum(s == 1)
    if (n1 < nsize) next
    ins <- s == 1; Rin <- Sig[ins, ins, drop = FALSE]; rho_bar <- max(0, mean(Rin[upper.tri(Rin)]))
    r <- vif_group(y, s, size, rho_bar, N)
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
  fwrite(res, file.path(out.dir, paste0("magma_drug_targetor_vif_", lv, "_res.csv")))
  cat("Wrote", lv, ":", nrow(res), "ATC classes (VIF-OLS,", N, "drugs)\n")
}
