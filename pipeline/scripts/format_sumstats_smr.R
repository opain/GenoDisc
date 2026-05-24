#!/usr/bin/Rscript

suppressMessages(library("optparse"))
option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")

)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

ss<-fread(paste0(outdir,'/results/',opt$gwas,'/gwas_sumstat/',opt$gwas,'.cleaned.gz'))

if(all(names(ss) != 'BETA')){
  ss$BETA<-log(ss$OR)
}

if(any(names(ss) != 'FREQ')){
  ss$FREQ<-ss$REF.FREQ
}

ss<-ss[,c('SNP','A1','A2','FREQ','BETA','SE','P','N'),with=F]
names(ss)<-c('SNP','A1','A2','freq','b','se','p','N')

fwrite(ss, paste0(outdir,'/results/',opt$gwas,'/gwas_sumstat/',opt$gwas,'.cleaned.cojo'), sep=' ', na='NA', quote=F)

