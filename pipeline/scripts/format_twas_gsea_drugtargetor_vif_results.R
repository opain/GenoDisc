#!/usr/bin/Rscript

# Drug-level VIF-OLS ATC-class enrichment (CAMERA-style; Wu & Smyth 2012).
# The recommended drug-level ATC test (see docs/atc_enrichment_estimator_validation.Rmd).
# Asks the SAME question as the legacy Wilcoxon and the DRUGSETS GLS -- are the drugs in
# an ATC class more directionally concordant with the disease TWAS signature than other
# drugs? -- but via an OLS class contrast whose SE is inflated for within-class correlation
# (SE = SE_OLS * sqrt(1 + (n1-1) * rho_bar), rho_bar = mean in-class drug-drug correlation).
# Unlike the GLS this does not invert Sigma, so it does not suffer the pseudo-inverse
# out-of-class leakage that inflated small mixed classes and diluted large coherent ones.
#
# Directional (two-sided): the class coefficient B on the per-drug signed T. B>0 = in-class
# drugs match the disease signature more than average ("Matches disease"); B<0 = oppose it
# ("Opposes disease"). Reported with a 95% CI. Sigma = the exact drug-drug correlation of the
# per-drug statistics (twas_gsea's --drug_corr_out). Drugs are aligned to Sigma by CID and
# multi-ATC drugs are exploded into every class they belong to.

suppressMessages(library("optparse"))
option_list = list(
  make_option("--twas", action="store", default=NA, type='character', help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character', help="PANEL [required]"),
  make_option("--config_file", action="store", default=NA, type='character', help="Config file [required]"),
  make_option("--pipeline_dir", action="store", default=NA, type="character", help="Pipeline dir [required]"),
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

## --- VIF-OLS group test (two-sided) -----------------------------------------
vif_group <- function(y, s, size, rho_bar, N) {
  X <- cbind(1, s, size, log(size)); K <- ncol(X)
  XtXi <- solve(crossprod(X)); B <- XtXi %*% crossprod(X, y)
  resid <- y - X %*% B
  sigma <- as.numeric(crossprod(resid)) / (N - K)
  vif <- 1 + (sum(s) - 1) * rho_bar
  se <- sqrt(sigma * XtXi[2, 2] * vif); tval <- B[2] / se; df <- N - K; tc <- qt(0.975, df)
  list(B = B[2], SE = se, T = tval, P = 2 * pt(-abs(tval), df),
       CI_lo = B[2] - tc * se, CI_hi = B[2] + tc * se)
}

## --- inputs -----------------------------------------------------------------
res <- fread(paste0(ddir, '/twas_gsea_drugtargetor_', opt$panel, '.competitive.clean.csv'))
res <- res[!is.na(T) & N_Mem_Avail >= 2]
# Align by CID. clean.csv punctuation-normalises the whole GeneSet (multi-CID drugs become
# "CID.111.222") while the correlation rownames keep original punctuation ("CID:111,222"),
# so normalise BOTH to digits-only tokens; otherwise multi-CID drugs (e.g. flupentixol,
# CID 5281881,5281878) are silently dropped -- this recovers the full drug set (alignment bug #1).
norm_cid <- function(x) gsub("[^0-9]+", "_", x)
res[, cid := norm_cid(sub(".*CID\\.", "", GeneSet))]

DC   <- readRDS(paste0(ddir, '/twas_gsea_drugtargetor_', opt$panel, '.drugcorr.rds'))
rcid <- norm_cid(sub(".*CID:", "", rownames(DC)))
common <- intersect(res$cid, rcid)
res <- res[match(common, cid)]; ri <- match(common, rcid); DC <- DC[ri, ri, drop = FALSE]

y <- res$T; size <- as.numeric(res$N_Mem_Avail); N <- length(y)

atc <- fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep = '!')
names(atc) <- c('Code', 'Name'); atc$Name <- tolower(atc$Name)

for (lv in levels) {
  k <- level_width[[lv]]
  # Canonical membership: drug's own ATC codes (clean.csv 'ATC', dot-joined), exploded, unique per drug.
  drug_cls <- lapply(strsplit(res$ATC, '.', fixed = TRUE), function(x) unique(substr(x, 1, k)))
  classes <- sort(unique(unlist(drug_cls))); classes <- classes[nchar(classes) == k]
  out <- rbindlist(lapply(classes, function(cl) {
    s <- as.numeric(vapply(drug_cls, function(x) cl %in% x, logical(1))); n1 <- sum(s)
    if (n1 < 5) return(NULL)
    ins <- s == 1; Rin <- DC[ins, ins, drop = FALSE]; rho_bar <- max(0, mean(Rin[upper.tri(Rin)]))
    g <- vif_group(y, s, size, rho_bar, N)
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
  fwrite(out[, ..keep], paste0(ddir, '/twas_gsea_drugtargetor_vif_', lv, '_', opt$panel, '_res.csv'))
  cat(sprintf("%s %s: %d classes tested (VIF-OLS, %d drugs aligned by CID)\n",
              opt$panel, lv, nrow(out), N))
}
