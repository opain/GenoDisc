#!/usr/bin/Rscript
library(data.table)
library(optparse)

option_list = list(
  make_option("--resdir", type="character", default="resources")
)
option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

tissue<-c('Basalganglia','Cerebellum','Cortex','Hippocampus','Spinalcord')
chr<-1:22

for(tissue_i in tissue){
  for(chr_i in chr){
    esi<-fread(paste0(opt$resdir, '/data/MetaBrain/',tissue_i,'/2020-05-26-',tissue_i,'-EUR-',chr_i,'-SMR-besd.esi'))
    esi$V2_new<-esi$V2
    esi$V2_new<-gsub(':.*','',gsub('.*rs','rs',esi$V2_new))
    esi$V2_new[esi$V2_new == 'rs']<-'nors'
    esi$V2<-esi$V2_new
    esi$V2_new<-NULL
    fwrite(esi, paste0(opt$resdir, '/data/MetaBrain/',tissue_i,'/2020-05-26-',tissue_i,'-EUR-',chr_i,'-SMR-besd.esi'), col.names=F, row.names=F, quote=F, sep='\t', na='NA')
  }
}
