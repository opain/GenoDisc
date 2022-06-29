#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(biomaRt)

# Read in snakefile to find cofgi file
snakefile<-readLines('Snakefile')

# Read in config file
configfile<-snakefile[grepl('config', snakefile)]
configfile<-gsub('\"','',gsub('.* \\"','',configfile))
config<-readLines(configfile)

# Determine which panels were requested
non_gtex_weights<-gsub('non_gtex_weights: ','',config[grepl('^non_gtex_weights:',config)])
non_gtex_weights<-gsub('\\[','', non_gtex_weights)
non_gtex_weights<-gsub('\\]','', non_gtex_weights)
non_gtex_weights<-gsub("'",'', non_gtex_weights)
non_gtex_weights<-gsub('"','', non_gtex_weights)
non_gtex_weights<-unlist(strsplit(non_gtex_weights, ','))

gtex_weights<-gsub('gtex_weights: ','',config[grepl('^gtex_weights:',config)])
gtex_weights<-gsub('\\[','', gtex_weights)
gtex_weights<-gsub('\\]','', gtex_weights)
gtex_weights<-gsub("'",'', gtex_weights)
gtex_weights<-gsub('"','', gtex_weights)
gtex_weights<-unlist(strsplit(gtex_weights, ','))

weights<-c(non_gtex_weights, gtex_weights)

if(config[grepl('^twas_panel_psychencode:',config)] == "twas_panel_psychencode: T"){
  weights<-c(weights, 'psychencode')
}

# Write out this list of SNP-weights as this might be useful elsewhere
write.table(weights, paste0('results/',opt$gwas,'/twas/list_of_weights.txt'), col.names=F, row.names=F, quote=F) 

all<-NULL

# Read in gene names from biomart
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)

# Read in and format PsychENCODE TWAS results
if(any(weights == 'psychencode')){
  twas_psychencode_files<-list.files(path=paste0('results/',opt$gwas,'/twas/psychencode/'), pattern=paste0(opt$twas,'_twas_psychencode_chr'))
  
  twas_psychencode<-NULL
  for(i in twas_psychencode_files){
    twas_psychencode<-rbind(twas_psychencode, fread(paste0('results/',opt$gwas,'/twas/psychencode/',i)))
  }

  twas_psychencode<-merge(twas_psychencode,Genes, by.x='ID', by.y='ensembl_gene_id', all.x=T)
  twas_psychencode<-twas_psychencode[!duplicated(twas_psychencode$ID),]
  names(twas_psychencode)[names(twas_psychencode) == 'ID']<-'ensembl_gene_id'
  
  twas_psychencode<-twas_psychencode[,c('FILE','PANEL','MODEL','BEST.GWAS.Z','ensembl_gene_id','external_gene_name','CHR','P0','P1','TWAS.Z','TWAS.P','COLOC.PP3','COLOC.PP4'), with=F]
  
  all<-rbind(all, twas_psychencode)

}

# Read in FUSION non-gtex weights
for(weight in weights[weights %in% c("CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","YFS.BLOOD.RNAARR")]){
  twas_fusion_files<-list.files(path=paste0('results/',opt$gwas,'/twas/fusion_',weight,'/'), pattern=paste0(opt$twas,'_twas_',weight,'_chr'))
  
  twas_fusion<-NULL
  for(i in twas_fusion_files){
    twas_fusion<-rbind(twas_fusion, fread(paste0('results/',opt$gwas,'/twas/fusion_',weight,'/',i)))
  }
  
  twas_fusion<-merge(twas_fusion,Genes, by.x='ID', by.y='external_gene_name', all.x=T)
  twas_fusion<-twas_fusion[!duplicated(twas_fusion$ensembl_gene_id),]
  names(twas_fusion)[names(twas_fusion) == 'ID']<-'external_gene_name'
  
  twas_fusion<-twas_fusion[,c('FILE','PANEL','MODEL','BEST.GWAS.Z','ensembl_gene_id','external_gene_name','CHR','P0','P1','TWAS.Z','TWAS.P','COLOC.PP3','COLOC.PP4'), with=F]
  
  all<-rbind(all, twas_fusion)
  
}

# Read in FUSION gtex weights
for(weight in weights[weights %in% c("Adrenal_Gland","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Substantia_nigra","Pituitary","Thyroid","Whole_Blood")]){
  twas_fusion_files<-list.files(path=paste0('results/',opt$gwas,'/twas/fusion_',weight,'/'), pattern=paste0(opt$twas,'_twas_',weight,'_chr'))
  
  twas_fusion<-NULL
  for(i in twas_fusion_files){
    twas_fusion<-rbind(twas_fusion, fread(paste0('results/',opt$gwas,'/twas/fusion_',weight,'/',i)))
  }
  
  twas_fusion$ID<-gsub('\\..*','',twas_fusion$ID)
  
  twas_fusion<-merge(twas_fusion,Genes, by.x='ID', by.y='ensembl_gene_id', all.x=T)
  twas_fusion<-twas_fusion[!duplicated(twas_fusion$ID),]
  names(twas_fusion)[names(twas_fusion) == 'ID']<-'ensembl_gene_id'
  
  twas_fusion<-twas_fusion[,c('FILE','PANEL','MODEL','BEST.GWAS.Z','ensembl_gene_id','external_gene_name','CHR','P0','P1','TWAS.Z','TWAS.P','COLOC.PP3','COLOC.PP4'), with=F]
  
  all<-rbind(all, twas_fusion)
}

all<-all[!duplicated(all$FILE),]

###
# Write out combined results
###

fwrite(all, paste0('results/',opt$gwas,'/twas/',opt$gwas,'_twas_GW_clean.txt.gz'), quote=F, sep=' ', na='NA')

# Write file listing chromosomes with transcriptome-wide significant results
all$TWAS.P.FDR<-p.adjust(all$TWAS.P, method='BY')
write.table(unique(all$CHR[which(all$TWAS.P.FDR < 0.05)]), paste0('results/',opt$gwas,'/twas/',opt$gwas,'_twas_sig_chr.txt'), row.names=F, col.names=T, quote=F)

# Write out transcriptome-wide significant results
names(all)[names(all) == 'external_gene_name']<-'ID'
write.table(all[which(all$TWAS.P.FDR < 0.05),], paste0('results/',opt$gwas,'/twas/',opt$gwas,'_twas_GW_clean_sig.txt'), row.names=F, col.names=T, quote=F)
