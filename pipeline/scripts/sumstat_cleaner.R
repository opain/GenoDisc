#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics file [required]"),
  make_option("--ref_chr", action="store", default=NA, type='character',
              help="Path to per chromosome reference .rds files [required]"),
  make_option("--population", action="store", default=NA, type='character',
              help="Reference population code [required]"),
  make_option("--sampling", action="store", default=NA, type='numeric',
              help="Sampling fraction of cases vs. controls [optional]"),
  make_option("--n", action="store", default=NA, type='numeric',
              help="Total sample size [optional]"),
  make_option("--info", action="store", default=0.9, type='numeric',
              help="INFO threshold [optional]"),
  make_option("--maf", action="store", default=0.01, type='numeric',
              help="MAF threshold [optional]"),
  make_option("--maf_diff", action="store", default=0.2, type='numeric',
              help="Difference between reference and reported MAF threshold [optional]"),
  make_option("--output", action="store", default='./Output', type='character',
              help="Path for output files [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'sumstat_cleaner.R', start.time = start.time)

#####
# Read in sumstats
#####

log_add(log_file = log_file, message = 'Reading in sumstats.')

GWAS <- fread(opt$sumstats)

log_add(log_file = log_file, message = paste0('GWAS contains ',nrow(GWAS),' variants.'))

# Insert N using opt$N if N data isn't already present
if(!(all(c('N_CAS','N_CON') %in% names(GWAS)) | 'N' %in% names(GWAS))){
  GWAS$N <- opt$n
}

#####
# Interpret and update header
#####

# Update header and drop columns that aren't recognised
GWAS <- format_header(sumstats = GWAS, log_file = log_file)

#####
# Format FREQ and N columns
#####

if(!('FREQ' %in% names(GWAS)) & all(c('FRQ_A','FRQ_U') %in% names(GWAS))){
  # Calculate average allele frequency across cases and controls
  GWAS$FREQ <- calc_mean_freq(sumstats = GWAS, sampling = opt$sampling, log_file = log_file)
}

# If N_CAS and N_CON are provided, but N is not, calculate total N
if(all(c('N_CAS','N_CON') %in% names(GWAS)) & !('N' %in% names(GWAS))){
  GWAS$N <- GWAS$N_CAS + GWAS$N_CON
}

# Retain only expected columns
GWAS <- GWAS[, names(GWAS) %in% c('CHR','BP','SNP','A1','A2','BETA','SE','OR','Z','FREQ','N','INFO','P'), with=F]

#####
# Insert IUPAC codes into target
#####

# Insert IUPAC codes into target
GWAS$IUPAC<-snp_iupac(GWAS$A1, GWAS$A2)

# Retain only non-ambiguous SNPs
GWAS <- remove_ambig(GWAS)

log_add(log_file = log_file, message = paste0('After removal of variants that are not SNPs or are ambiguous, ',nrow(GWAS),' variants remain.'))

#####
# Harmonise per chromosome with reference
#####
# We identify variants present in GWAS and reference, insert missing data (SNP, CHR, BP), and REF.FREQ

GWAS <- ref_harmonise(targ = GWAS, ref_rds = opt$ref_chr, population = opt$population, log_file = log_file)

#####
# Remove SNPs with INFO < opt$info
#####

GWAS <- filter_info(targ = GWAS, thresh = opt$info, log_file = log_file)

#####
# Remove SNPs with reported MAF < opt$maf
#####

GWAS <- filter_maf(targ = GWAS, thresh = opt$maf, ref = F, log_file = log_file)

#####
# Remove SNPs with reference MAF < opt$maf
#####

GWAS <- filter_maf(targ = GWAS, thresh = opt$maf, ref = T, log_file = log_file)

#####
# Remove SNPs with discordant MAF
#####

GWAS <- discord_maf(targ = GWAS, thresh = opt$maf_diff, log_file = log_file, plot_file = paste0(opt$output,'.MAF_plot.png'))

#####
# Remove SNPs with out-of-bounds p-values
#####

GWAS<-GWAS[GWAS$P <= 1 & GWAS$P > 0,]
log_add(log_file = log_file, message = paste0('After removal of SNPs with out-of-bound P values, ',nrow(GWAS),' variants remain.'))

#####
# Remove SNPs with duplicated rs numbers
#####

dups<-GWAS$SNP[duplicated(GWAS$SNP)]
GWAS<-GWAS[!(GWAS$SNP %in% dups),]
log_add(log_file = log_file, message = paste0('After removal of SNPs with duplicate IDs, ',nrow(GWAS),' variants remain.'))

#####
# Remove SNPs with N < 3SD from median N
#####

GWAS <- filter_n(targ = GWAS, log_file = log_file)

#####
# Create BETA column if not present
#####

GWAS <- insert_beta(targ = GWAS, log_file = log_file)

# Remove OR and Z columns if present
GWAS<-GWAS[, names(GWAS) %in% c('CHR','BP','SNP','A1','A2','BETA','SE','FREQ','REF.FREQ','N','INFO','P'), with=F]

#####
# Check for genomic control
#####

GWAS<-avoid_gc(targ = GWAS, log_file = log_file)

#####
# Insert SE column
#####

GWAS <- insert_se(targ = GWAS, log_file = log_file)

#####
# Remove SNPs with SE == 0
#####

GWAS<-GWAS[GWAS$SE != 0,]
log_add(log_file = log_file, message = paste0('After removal of SNPs with SE == 0, ',nrow(GWAS),' variants remain.'))

#####
# Write out results
#####

if(file.exists(paste0(opt$output,'.gz'))){
  system(paste0(paste0('rm ',opt$output,'.gz')))
}

fwrite(GWAS, paste0(opt$output,'.gz'), sep='\t')

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat0('Analysis finished at ',as.character(end.time),'\n')
cat0('Analysis duration was ',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
