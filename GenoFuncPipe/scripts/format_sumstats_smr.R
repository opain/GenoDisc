#!/usr/bin/Rscript

suppressMessages(library("optparse"))
option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

ss<-fread(paste0('resources/data/gwas_sumstat/',opt$gwas,'/',opt$gwas,'.cleaned.gz'))

if(any(names(ss) != 'BETA')){
  ss$BETA<-log(ss$OR)
}

if(any(names(ss) != 'FREQ')){
  ss$FREQ<-ss$REF.FREQ
}

ss<-ss[,c('SNP','A1','A2','FREQ','BETA','SE','P','N'),with=F]
names(ss)<-c('SNP','A1','A2','freq','b','se','p','N')

fwrite(ss, paste0('resources/data/gwas_sumstat/',opt$gwas,'/',opt$gwas,'.cleaned.cojo'), sep=' ', na='NA', quote=F)

