#!/usr/bin/Rscript
library(data.table)

# Read in database
pathways<-fread('resources/data/drug_targetor/wholedatabase_for_targetor')

# Read in gene list
gene_id<-fread('resources/data/magma/NCBI37.3.gene.loc')
gene_id<-gene_id[,'V6', with=F]
names(gene_id)<-'ID'

drugs<-unique(pathways$atc)
for(i in 1:length(drugs)){
  drug_tmp<-pathways[pathways$atc == drugs[i],]
  
  gene_id[[drugs[i]]]<-ifelse(gene_id$ID %in% drug_tmp$gene[drug_tmp$activity_type %in% c('DECREASED_EXPRESSION','NEGATIVE_RESPONSE','OPPOSITE_RESPONSE')], -1, ifelse(gene_id$ID %in% drug_tmp$gene[drug_tmp$activity_type %in% c('INCREASED_EXPRESSION','POSITIVE_RESPONSE')], 1, 0))
}

fwrite(gene_id, 'resources/data/drug_targetor/wholedatabase_for_targetor_directional.prop', sep=',', quote=T, na='NA')
