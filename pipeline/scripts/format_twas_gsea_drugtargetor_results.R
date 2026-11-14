#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--twas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character',
              help="PANEL [required]"),
  make_option("--mode", action="store", default='directional', type='character',
              help="Either 'directional' (default) or 'nondirectional'"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

if(!(opt$mode %in% c('directional','nondirectional'))) stop("--mode must be 'directional' or 'nondirectional'")
suffix <- if(opt$mode == 'nondirectional') '_nondir' else ''

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

# Read in TWAS-GSEA results.
# The engine writes a space-delimited table; the drug GeneSet names contain
# spaces (e.g. "...NAME:LIOTHYRONINE SODIUM|...") and are not always quoted, so
# a plain space-delimited read mis-splits those rows. Read robustly: if the
# parsed column count doesn't match the header, peel the fixed number of
# trailing numeric columns off each line and treat the remainder as GeneSet.
res_file <- paste0(outdir,'/results/',opt$twas,'/twas/drugtargetor/twas_gsea_drugtargetor',suffix,'_',opt$panel,'.competitive.txt')
res_hdr  <- strsplit(trimws(readLines(res_file, n = 1L)), '\\s+')[[1]]
res<-fread(res_file, sep = ' ', quote = '"', fill = TRUE)
if(ncol(res) != length(res_hdr)){
  res_raw <- readLines(res_file)[-1]
  res_raw <- res_raw[nzchar(res_raw)]
  res_k   <- length(res_hdr) - 1L
  res_pat <- paste0('^(.*?)\\s+((?:\\S+\\s+){', res_k - 1L, '}\\S+)\\s*$')
  res_gs  <- sub(res_pat, '\\1', res_raw)
  res_num <- do.call(rbind, strsplit(trimws(sub(res_pat, '\\2', res_raw)), '\\s+'))
  res     <- as.data.table(matrix(as.numeric(res_num), nrow = nrow(res_num)))
  res     <- cbind(GeneSet = res_gs, res)
  setnames(res, res_hdr)
}

if(opt$mode == 'directional'){
  # Both signs of T are interpretable (drug mimics vs reverses disease signature),
  # so use a two-sided p-value.
  res$P <- 2*pnorm(-abs(res$T))
  wilcox_alt <- 'two.sided'
} else {
  # TWAS-GSEA-fast already emits a one-sided right-tail P (probit(1-P) outcome),
  # which is what we want: only positive enrichment of drug-target genes is meaningful.
  wilcox_alt <- 'greater'
}

res$P.CORR<-p.adjust(res$P, method='fdr')

res$GeneSet <- gsub('[[:punct:]]', '.', res$GeneSet)

res$ATC<-gsub('ATC.','',gsub('\\.NAME.*','',res$GeneSet))
res$NAME<-tolower(gsub('\\.CID\\..*','',gsub('.*NAME\\.','',res$GeneSet)))
res$NAME<-gsub('_', ' ', res$NAME)

# Direction-of-effect columns.
# Directional: T is signed (T = Estimate/SE). T > 0 = drug-target gene-set
# direction matches disease TWAS-Z = mimics. T < 0 = opposes (reversal,
# repurposing direction). Reversal_Z is positive when the drug opposes disease,
# so any downstream consumer can use it directly without sign flips.
# Non-directional: P is one-sided right-tail enrichment, sign of T is not
# interpretable as mimics/opposes; Direction is NA and Reversal_Z = -qnorm(P)
# (positive = enriched, no direction).
if(opt$mode == 'directional'){
  res$Direction <- ifelse(res$T < 0, 'Opposes disease',
                   ifelse(res$T > 0, 'Matches disease', NA_character_))
  res$Reversal_Z <- -res$T
} else {
  res$Direction  <- NA_character_
  res$Reversal_Z <- -qnorm(res$P)
}

# Sort by p-value
res<-res[order(res$P),]

write.csv(res, paste0(outdir,'/results/',opt$twas,'/twas/drugtargetor/twas_gsea_drugtargetor',suffix,'_',opt$panel,'.competitive.clean.csv'), row.names=F)

# Insert ATC codes
atc<-fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep='!')
names(atc)<-c('Code','Name')
atc$Name<-tolower(atc$Name)

# Canonical membership: each drug's own ATC codes (clean.csv 'ATC', dot-joined), exploded to every
# class, unique per drug. The rank-sum runs over UNIQUE drugs (no drug*code pseudo-replication).
# Test for enrichment at each ATC level (L2/L3/L4 = 1st 3/4/5 chars), threshold >= 5 drugs, writing
# one file per level (mirrors format_twas_gsea_drugtargetor_vif_results.R). The packaging reader
# stacks the per-level files and adds a Level column + per-level FDR.
#
# Direction (directional mode): Wilcoxon HL on rank(T) ~ class_bin takes out-class as x and in-class
# as y, giving HL(out - in); so Estimate > 0 => in-class drugs have LOWER T = 'Opposes disease',
# Estimate < 0 => 'Matches disease'; Reversal_Z = qnorm(1 - P/2) * sign(Estimate). Non-directional:
# P is one-sided right-tail (in > out), Direction NA, Reversal_Z = qnorm(1 - P).
level_width <- c(l2 = 3L, l3 = 4L, l4 = 5L)
for (lv in names(level_width)) {
  k <- level_width[[lv]]
  drug_cls <- lapply(strsplit(res$ATC, '.', fixed = TRUE), function(x) unique(substr(x, 1, k)))
  cls      <- sort(unique(unlist(drug_cls))); cls <- cls[nchar(cls) == k]
  # Empty-panel guard: a panel with few testable drug sets may clear no class; start from an empty
  # typed frame so merge()/direction columns operate on a valid 0-row table and emit a well-formed CSV.
  atc_enrich <- data.frame(ATC=character(), Estimate=numeric(), Class_Median=numeric(),
                           Non_Class_Median=numeric(), P=numeric(), N=numeric())
  for (cat in cls) {
    class_bin <- as.numeric(vapply(drug_cls, function(x) cat %in% x, logical(1)))
    if (sum(class_bin == 1) >= 5) {
      if (opt$mode == 'directional') {
        wil_cox_res <- wilcox.test(rank(res$T) ~ class_bin, conf.int = T)
      } else {
        in_T  <- rank(res$T)[class_bin == 1]
        out_T <- rank(res$T)[class_bin == 0]
        wil_cox_res <- wilcox.test(in_T, out_T, conf.int = T, alternative = 'greater')
      }
      atc_enrich <- rbind(atc_enrich, data.frame(ATC = cat,
                                                 Estimate = as.numeric(wil_cox_res$estimate),
                                                 Class_Median = median(res$Estimate[class_bin == 1]),
                                                 Non_Class_Median = median(res$Estimate[class_bin == 0]),
                                                 P = wil_cox_res$p.value,
                                                 N = sum(class_bin)))
    }
  }

  atc_labels <- atc[nchar(atc$Code) == k, ]
  atc_enrich <- merge(atc_enrich, atc_labels, by.x = 'ATC', by.y = 'Code')
  atc_enrich <- atc_enrich[order(atc_enrich$P), ]

  if (opt$mode == 'directional') {
    atc_enrich$Direction  <- ifelse(atc_enrich$Estimate > 0, 'Opposes disease',
                             ifelse(atc_enrich$Estimate < 0, 'Matches disease', NA_character_))
    atc_enrich$Reversal_Z <- qnorm(1 - atc_enrich$P/2) * sign(atc_enrich$Estimate)
  } else {
    atc_enrich$Direction  <- NA_character_
    atc_enrich$Reversal_Z <- qnorm(1 - atc_enrich$P)
  }

  write.csv(atc_enrich, paste0(outdir, '/results/', opt$twas, '/twas/drugtargetor/twas_gsea', suffix, '_', opt$panel, '_res_atc_', lv, '_res.csv'), row.names = F)
}

