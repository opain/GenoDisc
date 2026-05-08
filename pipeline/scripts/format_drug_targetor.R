#!/usr/bin/Rscript
library(data.table)
library(optparse)

option_list = list(
  make_option("--resdir", type="character", default="resources")
)
opt = parse_args(OptionParser(option_list=option_list))

pathways<-fread(paste0(opt$resdir, '/data/drug_targetor/wholedatabase_for_targetor'))
pathways$atc<-gsub(' ', '_', pathways$atc)
gene_id<-fread(paste0(opt$resdir, '/data/magma/NCBI37.3.gene.loc'))
gene_id<-gene_id[,c('V1','V6'),with=F]
names(gene_id)<-c('GENE','ID')
pathways<-merge(pathways,gene_id, by.x='gene',by.y='ID')
drugs<-unique(pathways$atc)
drugs<-gsub(' ', '_', drugs)
gmt_entrez<-paste0(opt$resdir, '/data/drug_targetor/wholedatabase_for_targetor.gmt')
gmt_symbol<-paste0(opt$resdir, '/data/drug_targetor/wholedatabase_for_targetor_symbols.gmt')
unlink(gmt_entrez)
unlink(gmt_symbol)
for(i in 1:length(drugs)){
  drug_tmp<-pathways[pathways$atc == drugs[i],]
  drug_tmp<-drug_tmp[!duplicated(drug_tmp$gene),]
  # entrez-id .gmt (used by MAGMA)
  drug_row_entrez<-t(data.frame(c(drug_tmp$atc[1], ' ', drug_tmp$GENE)))
  write.table(drug_row_entrez, gmt_entrez, append=T, col.names=F, row.names=F, quote=F, sep='\t')
  # gene-symbol .gmt (used by TWAS-GSEA, which matches on symbols).
  # Sanitise the gene-set ID with make.names() so it (a) survives the
  # whitespace-delimited TWAS-GSEA output and (b) matches the dot-separated
  # column-name format produced by the directional .prop pipeline. This lets
  # downstream parsing be identical for both modes.
  drug_row_symbol<-t(data.frame(c(make.names(drug_tmp$atc[1]), ' ', drug_tmp$gene)))
  write.table(drug_row_symbol, gmt_symbol, append=T, col.names=F, row.names=F, quote=F, sep='\t')
}
