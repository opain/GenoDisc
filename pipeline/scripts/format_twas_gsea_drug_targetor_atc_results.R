#!/usr/bin/Rscript

# Clean the ATC-class-level TWAS-GSEA output ("option C": each ATC class is one
# gene set, tested at the gene level by TWAS-GSEA-fast.R). Unlike
# format_twas_gsea_drugtargetor_results.R this does NOT aggregate per-drug T with
# a Wilcoxon -- the class competitive.txt IS the class test. Here we only parse
# the ATC code, attach the label, set the direction/FDR, and write a tidy CSV.

suppressMessages(library("optparse"))

option_list = list(
  make_option("--twas", action="store", default=NA, type='character', help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character', help="PANEL [required]"),
  make_option("--level", action="store", default=NA, type='character', help="ATC level: l2, l3 or l4 [required]"),
  make_option("--mode", action="store", default='directional', type='character', help="'directional' (default) or 'nondirectional'"),
  make_option("--config_file", action="store", default=NA, type='character', help="Path to config file [required]"),
  make_option("--pipeline_dir", action="store", default=NA, type="character", help="Path to the pipeline directory [required]")
)
opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

if(!(opt$mode %in% c('directional','nondirectional'))) stop("--mode must be 'directional' or 'nondirectional'")
suffix <- if(opt$mode == 'nondirectional') '_nondir' else ''

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

config <- readLines(opt$config_file)
outdir <- gsub('outdir: ', '', config[grepl('outdir: ', config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

base <- paste0(outdir, '/results/', opt$twas, '/twas/drugtargetor/twas_gsea_drugtargetor_atc_', opt$level, suffix, '_', opt$panel)
res <- fread(paste0(base, '.competitive.txt'), sep = ' ')

# Strip the "ATC_" prefix to recover the bare class code. The non-directional
# .gmt path comes back from TWAS-GSEA with the separator munged to a dot
# ("ATC.C03C"), so accept either separator or the label merge fails (Name = NA).
res$Code <- sub('^ATC[._]', '', res$GeneSet)

# Direction-of-effect columns, same convention as the per-drug clean file.
# Directional: TWAS-GSEA emits a one-sided P; recompute two-sided (both mimics
# and reversal are interpretable). T>0 = class up-regulated targets track disease
# TWAS-Z = Matches; T<0 = Opposes (reversal). Reversal_Z positive = opposes.
# Non-directional: keep TWAS-GSEA's one-sided right-tail P; sign not interpretable.
if(opt$mode == 'directional'){
  res$P <- 2 * pnorm(-abs(res$T))
  res$Direction <- ifelse(res$T < 0, 'Opposes disease',
                   ifelse(res$T > 0, 'Matches disease', NA_character_))
  res$Reversal_Z <- -res$T
} else {
  res$Direction  <- NA_character_
  res$Reversal_Z <- -qnorm(res$P)
}
res$FDR <- p.adjust(res$P, method = 'fdr')
res$P.CORR <- res$FDR

# Attach ATC class label
atc <- fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep = '!')
names(atc) <- c('Code', 'Name')
atc$Name <- tolower(atc$Name)
res <- merge(res, atc, by = 'Code', all.x = TRUE)

res <- res[order(res$P), ]

keep <- intersect(c('GeneSet','Code','Name','Estimate','SE','T','N_Mem_Avail','N_Mem',
                    'P','FDR','P.CORR','Direction','Reversal_Z'), names(res))
res <- res[, ..keep]

write.csv(res, paste0(base, '_res.csv'), row.names = FALSE)
