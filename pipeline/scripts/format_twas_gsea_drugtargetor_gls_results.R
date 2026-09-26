#!/usr/bin/Rscript

# Drug-level GLS ATC-class enrichment (DRUGSETS-style; Bell et al. 2022,
# github.com/nybell/drugsets). Asks the SAME question as the legacy per-drug
# Wilcoxon -- are the drugs in an ATC class more associated than other drugs? --
# but validly, by modelling the correlation between drugs that share target genes.
#
# For each ATC level (L2/L3/L4): GLS regression of the per-drug signed TWAS-GSEA T
# on an ATC-class indicator (+ gene-set-size covariates), weighted by the exact
# drug-drug correlation Sigma of the per-drug statistics (cov2cor(Z_wh'Z_wh),
# dumped by TWAS-GSEA's --drug_corr_out on the directional per-drug run). Sigma is
# inverted keeping eigenvalues >= --floor. The class coefficient is DIRECTIONAL
# (B>0 = in-class drugs match the disease signature more than average) and is
# reported with a 95% CI (read the CI, not the sign, near the null).

suppressMessages(library("optparse"))
option_list = list(
  make_option("--twas", action="store", default=NA, type='character', help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character', help="PANEL [required]"),
  make_option("--config_file", action="store", default=NA, type='character', help="Config file [required]"),
  make_option("--pipeline_dir", action="store", default=NA, type="character", help="Pipeline dir [required]"),
  make_option("--floor", action="store", default=0.1, type='numeric', help="Eigenvalue floor for the Sigma pseudo-inverse"),
  make_option("--levels", action="store", default="l2,l3,l4", type='character', help="ATC levels to emit")
)
opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

suppressMessages(library(data.table))
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

config <- readLines(opt$config_file)
outdir <- gsub('outdir: ', '', config[grepl('outdir: ', config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)
ddir   <- paste0(outdir, '/results/', opt$twas, '/twas/drugtargetor')

level_width <- c(l2 = 3L, l3 = 4L, l4 = 5L)
levels <- trimws(strsplit(opt$levels, ",", fixed = TRUE)[[1]])

## --- GLS group test (DRUGSETS compute_lnreg.R algebra) + 95% CI --------------
reg_inv <- function(S, floor) {
  e <- eigen(S, symmetric = TRUE); keep <- e$values >= floor
  list(inv = e$vectors[, keep, drop=FALSE] %*% (t(e$vectors[, keep, drop=FALSE]) / e$values[keep]),
       rank = sum(keep))
}
gls_group <- function(y, s, size, Sinv) {
  X <- cbind(1, s, size, log(size)); N <- length(y); K <- ncol(X)
  Xt <- t(X); W <- solve(Xt %*% Sinv %*% X)
  B <- W %*% Xt %*% Sinv %*% y
  resid <- y - X %*% B
  sigma <- as.numeric(t(resid) %*% Sinv %*% resid) / (N - K)
  se <- sqrt(sigma * W[2, 2]); tval <- B[2] / se; df <- N - K; tc <- qt(0.975, df)
  list(B = B[2], SE = se, T = tval, P = 2 * pt(-abs(tval), df),
       CI_lo = B[2] - tc * se, CI_hi = B[2] + tc * se)
}

## --- inputs -----------------------------------------------------------------
res <- fread(paste0(ddir, '/twas_gsea_drugtargetor_', opt$panel, '.competitive.clean.csv'))
res <- res[!is.na(T) & N_Mem_Avail >= 2]
res$Code7 <- res$ATC

DC  <- readRDS(paste0(ddir, '/twas_gsea_drugtargetor_', opt$panel, '.drugcorr.rds'))
dc7 <- sub("^ATC:([^|]+)\\|.*", "\\1", rownames(DC))
mi  <- match(res$Code7, dc7)
res <- res[!is.na(mi)]; mi <- mi[!is.na(mi)]
Sigma <- DC[mi, mi, drop = FALSE]
ri <- reg_inv(Sigma, opt$floor); Sinv <- ri$inv

atc <- fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep = '!')
names(atc) <- c('Code', 'Name'); atc$Name <- tolower(atc$Name)

y <- res$T; size <- res$N_Mem_Avail

for (lv in levels) {
  k <- level_width[[lv]]
  res[[paste0('cls', lv)]] <- substr(res$Code7, 1, k)
  classes <- sort(unique(res[[paste0('cls', lv)]]))
  out <- rbindlist(lapply(classes, function(cl) {
    s <- as.numeric(res[[paste0('cls', lv)]] == cl); n1 <- sum(s)
    if (n1 < 2) return(NULL)
    g <- gls_group(y, s, size, Sinv)
    data.table(Code = cl, N_Drugs = n1,
               Estimate = g$B, SE = g$SE, CI_lo = g$CI_lo, CI_hi = g$CI_hi,
               T = g$T, P = g$P,
               Direction = ifelse(g$B > 0, 'Matches disease', ifelse(g$B < 0, 'Opposes disease', NA_character_)),
               Reversal_Z = -g$T)   # heatmap convention: positive = opposes disease
  }))
  out[, FDR := p.adjust(P, method = 'fdr')]
  out <- merge(out, atc[nchar(Code) == k, ], by = 'Code', all.x = TRUE)
  out <- out[order(P)]
  keep <- c('Code','Name','N_Drugs','Estimate','SE','CI_lo','CI_hi','T','P','FDR','Direction','Reversal_Z')
  fwrite(out[, ..keep], paste0(ddir, '/twas_gsea_drugtargetor_gls_', lv, '_', opt$panel, '_res.csv'))
  cat(sprintf("%s %s: %d classes tested (Sigma eff.rank %d/%d @floor %.2f)\n",
              opt$panel, lv, nrow(out), ri$rank, nrow(Sigma), opt$floor))
}
