#!/usr/bin/Rscript

suppressMessages(library("optparse"))
option_list = list(
  make_option("--twas", action="store", default=NA, type='character',
              help="GWAS ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

twas_psychencode_files<-list.files(path=paste0('results/',opt$twas,'/twas/psychencode/'), pattern=paste0(opt$twas,'_twas_psychencode_chr'))

twas_psychencode<-NULL
for(i in twas_psychencode_files){
  twas_psychencode<-rbind(twas_psychencode, fread(paste0('results/',opt$twas,'/twas/psychencode/',i)))
}

# Insert gene IDs
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)

Genes<-Genes[!duplicated(Genes$ensembl_gene_id),]

twas_psychencode<-twas_psychencode[!duplicated(twas_psychencode$ID),]

twas_psychencode<-merge(twas_psychencode,Genes, by.x='ID', by.y='ensembl_gene_id', all.x=T)
names(twas_psychencode)[names(twas_psychencode) == 'ID']<-'ensembl_gene_id'

twas_psychencode<-twas_psychencode[,c('ensembl_gene_id','external_gene_name','CHR','P0','P1','TWAS.Z','TWAS.P','COLOC.PP3','COLOC.PP4'), with=F]

fwrite(twas_psychencode, paste0('results/',opt$twas,'/twas/psychencode/',opt$twas,'_twas_psychencode_GW_clean.txt.gz'), sep=' ', na='NA', quote=F)

