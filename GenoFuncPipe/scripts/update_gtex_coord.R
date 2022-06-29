#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--panel", action="store", default=NA, type='character',
              help="Panel ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Read in gene locations from build GRCh37
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart = ensembl)

library(data.table)

pos<-fread(paste0('resources/data/fusion_data/',opt$panel,'/',opt$panel,'.pos'))
pos$ensembl_gene_id<-gsub('\\..*','',pos$ID)
pos<-merge(pos, Genes, by='ensembl_gene_id')
pos<-pos[!duplicated(pos$ID),]
pos<-pos[,c('PANEL','WGT','ID','CHR','start_position','end_position','N'), with=F]
names(pos)<-c('PANEL','WGT','ID','CHR','P0','P1','N')
fwrite(pos, paste0('resources/data/fusion_data/',opt$panel,'/',opt$panel,'.pos'), quote=F, sep=' ', na='NA')
