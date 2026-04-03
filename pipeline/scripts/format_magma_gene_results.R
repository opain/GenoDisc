#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
source('scripts/functions/utils_functions.R')

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

# Read in MAGMA gene results
res<-fread(paste0(outdir,'/results/',opt$gwas,'/magma/magma_gene_level.genes.out'))

# Insert gene ID
loc<-fread(paste0(resdir, '/data/magma/NCBI37.3.gene.loc'))
loc<-loc[,c('V1','V6'),with=F]
names(loc)<-c('GENE','ID')

res<-merge(res, loc, by='GENE')

res<-res[,c('CHR','START','STOP','ID','P'),with=F]

write.csv(res, paste0(outdir,'/results/',opt$gwas,'/magma/magma_gene_level.clean.csv'), row.names=F, quote=F)

