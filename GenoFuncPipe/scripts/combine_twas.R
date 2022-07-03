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

# Read in gene names from biomart
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)

# Read in TWAS resukts and insert external_gene_name
all<-NULL
for(weight_i in weights){
  twas_files<-list.files(path=paste0('results/',opt$gwas,'/twas/',weight_i,'/'), pattern=paste0(opt$twas,'_twas_',weight_i,'_chr'))
  
  twas<-NULL
  for(i in twas_files){
    twas<-rbind(twas, fread(paste0('results/',opt$gwas,'/twas/',weight_i,'/',i)))
  }
  
  twas<-merge(twas,Genes, by.x='ID', by.y='ensembl_gene_id', all.x=T)
  twas<-twas[!duplicated(twas$ID),]
  names(twas)[names(twas) == 'ID']<-'ensembl_gene_id'
  
  twas<-twas[,c('FILE','PANEL','MODEL','BEST.GWAS.Z','NSNP','NWGT','MODELCV.R2','ensembl_gene_id','external_gene_name','CHR','P0','P1','TWAS.Z','TWAS.P','COLOC.PP3','COLOC.PP4'), with=F]
  
  # Create ID column for compatability with post TWAS analyses
  twas$ID<-twas$external_gene_name
  
  fwrite(twas, paste0('results/',opt$gwas,'/twas/',opt$gwas,'_twas_',weight_i,'_GW_clean.txt.gz'), quote=F, sep=' ', na='NA')
  
  all<-rbind(all, twas)
  
}

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
