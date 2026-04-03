#!/usr/bin/Rscript
library(data.table)
library(optparse)

option_list = list(
  make_option("--resdir", type="character", default="resources")
)
opt = parse_args(OptionParser(option_list=option_list))

pathways<-fread(paste0(opt$resdir, '/data/drug_targetor/wholedatabase_for_targetor'))
gene_id<-fread(paste0(opt$resdir, '/data/magma/NCBI37.3.gene.loc'))
gene_id<-gene_id[,c('V1','V6'),with=F]
names(gene_id)<-c('GENE','ID')
pathways<-merge(pathways,gene_id, by.x='gene',by.y='ID')
drugs<-unique(pathways$atc)
for(i in 1:length(drugs)){
  drug_tmp<-pathways[pathways$atc == drugs[i],]
  drug_tmp<-drug_tmp[!duplicated(drug_tmp$gene),]
  drug_row<-t(data.frame(c(drug_tmp$atc[1], ' ', drug_tmp$GENE)))
  write.table(drug_row, paste0(opt$resdir, '/data/drug_targetor/wholedatabase_for_targetor.gmt'), append=T, col.names=F, row.names=F, quote=F, sep='\t')
}
