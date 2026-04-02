#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--panel", action="store", default=NA, type='character',
              help="Panel ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

panels<-c("CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","YFS.BLOOD.RNAARR","METSIM.ADIPOSE.RNASEQ")

panels<-data.frame(panel=panels,
                   N=c(452,452,1247,1264,563))

library(data.table)

pos<-fread(paste0('resources/data/fusion_snp_weights/',opt$panel,'/',opt$panel,'.pos'))
pos$PANEL<-opt$panel
pos$N<-panels$N[panels$panel == opt$panel]

# Convert IDs to ENSEMBL IDs
biomart<-read.delim('resources/data/biomart/biomart_genes_grch37.tsv', stringsAsFactors=FALSE)
Genes<-biomart[,c('ensembl_gene_id','external_gene_name','external_synonym')]

pos_1<-merge(pos,Genes[,c('ensembl_gene_id','external_gene_name')], by.x='ID', by.y='external_gene_name')
pos_2<-merge(pos,Genes[,c('ensembl_gene_id','external_synonym')], by.x='ID', by.y='external_synonym')

pos_new<-do.call(rbind, list(pos_1, pos_2))
pos_new<-pos_new[!duplicated(pos_new$WGT),]

pos_new$ID<-pos_new$ensembl_gene_id
pos_new$ensembl_gene_id<-NULL

pos_new<-pos_new[,c('PANEL','WGT','ID','CHR','P0','P1','N'), with=F]
write.table(pos_new, paste0('resources/data/fusion_snp_weights/',opt$panel,'/',opt$panel,'.pos'), quote=F, col.names=T, row.names=F)
