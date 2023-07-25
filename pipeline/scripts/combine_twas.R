#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(biomaRt)

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

gtex_weights<-config[grepl('^gtex_weights', config)]
gtex_weights<-unlist(strsplit(gsub('"','',gsub('\\]','',gsub('.*\\[','',gtex_weights))),','))

non_gtex_weights<-config[grepl('^non_gtex_weights', config)]
non_gtex_weights<-unlist(strsplit(gsub('"','',gsub('\\]','',gsub('.*\\[','',non_gtex_weights))),','))

weights<-c(gtex_weights, non_gtex_weights)

psychencode_weights_log<-config[grepl('^twas_panel_psychencode:', config)]
if(psychencode_weights_log == "twas_panel_psychencode: T"){
  weights<-c(weights, 'psychencode')
}

external_weights_log<-config[grepl('^external_weights:', config)]
if(external_weights_log == "external_weights: T"){
  external_weights<-config[grepl('^external_weights_pos_path', config)]
  external_weights<-gsub('.pos','',basename(unlist(strsplit(gsub('"','',gsub('\\]','',gsub('.*\\[','',external_weights))),','))))
  weights<-c(weights, external_weights)
}

# Write out this list of SNP-weights as this might be useful elsewhere
write.table(weights, paste0(outdir,'/results/',opt$gwas,'/twas/list_of_weights.txt'), col.names=F, row.names=F, quote=F) 

# Read in gene names from biomart
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
biomartCacheClear()
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)

# Read in TWAS resukts and insert external_gene_name
all<-NULL
for(weight_i in weights){
  twas_files<-list.files(path=paste0(outdir,'/results/',opt$gwas,'/twas/',weight_i,'/'), pattern=paste0(opt$twas,'_twas_',weight_i,'_chr'))
  
  twas<-NULL
  for(i in twas_files){
    twas<-rbind(twas, fread(paste0(outdir,'/results/',opt$gwas,'/twas/',weight_i,'/',i)))
  }
  
  twas<-merge(twas,Genes, by.x='ID', by.y='ensembl_gene_id', all.x=T)
  twas<-twas[!duplicated(twas$ID),]
  names(twas)[names(twas) == 'ID']<-'ensembl_gene_id'
  
  twas<-twas[,c('FILE','PANEL','MODEL','BEST.GWAS.Z','NSNP','NWGT','MODELCV.R2','ensembl_gene_id','external_gene_name','CHR','P0','P1','TWAS.Z','TWAS.P','COLOC.PP3','COLOC.PP4'), with=F]
  
  # Create ID column for compatability with post TWAS analyses
  twas$ID<-twas$external_gene_name
  
  fwrite(twas, paste0(outdir,'/results/',opt$gwas,'/twas/',opt$gwas,'_twas_',weight_i,'_GW_clean.txt.gz'), quote=F, sep=' ', na='NA')
  
  all<-rbind(all, twas)
  
}

###
# Write out combined results
###

fwrite(all, paste0(outdir,'/results/',opt$gwas,'/twas/',opt$gwas,'_twas_GW_clean.txt.gz'), quote=F, sep=' ', na='NA')

# Write file listing chromosomes with transcriptome-wide significant results
all$TWAS.P.FDR<-p.adjust(all$TWAS.P, method='fdr')
write.table(unique(all$CHR[which(all$TWAS.P.FDR < 0.05)]), paste0(outdir,'/results/',opt$gwas,'/twas/',opt$gwas,'_twas_sig_chr.txt'), row.names=F, col.names=T, quote=F)

# Write out transcriptome-wide significant results
names(all)[names(all) == 'external_gene_name']<-'ID'
write.table(all[which(all$TWAS.P.FDR < 0.05),], paste0(outdir,'/results/',opt$gwas,'/twas/',opt$gwas,'_twas_GW_clean_sig.txt'), row.names=F, col.names=T, quote=F)
