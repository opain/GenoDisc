#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--panel", action="store", default=NA, type='character',
              help="Panel ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

panels<-c("CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","YFS.BLOOD.RNAARR")

panels<-data.frame(panel=panels,
                   N=c(452,452,1247,1264))

library(data.table)

pos<-fread(paste0('resources/data/fusion_data/',opt$panel,'/',opt$panel,'.pos'))
pos$PANEL<-opt$panel
pos$N<-panels$N[panels$panel == opt$panel]
pos<-pos[,c('PANEL','WGT','ID','CHR','P0','P1','N'), with=F]
write.table(pos, paste0('resources/data/fusion_data/',opt$panel,'/',opt$panel,'.pos'), quote=F, col.names=T, row.names=F)
