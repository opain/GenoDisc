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
source('scripts/functions/utils_functions.R')
source_all('scripts/functions')

# Read in config parameters
outdir <- read_param(config = opt$config_file, param = 'outdir', return_obj = F)
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

gtex_weights <- read_param(config = opt$config_file, param = 'gtex_weights', return_obj = F)
non_gtex_weights <- read_param(config = opt$config_file, param = 'non_gtex_weights', return_obj = F)

weights<-c(gtex_weights, non_gtex_weights)

twas_panel_psychencode <- read_param(config = opt$config_file, param = 'twas_panel_psychencode', return_obj = F)
if(twas_panel_psychencode == "T"){
  weights<-c(weights, 'psychencode')
}

external_weights_flag <- read_param(config = opt$config_file, param = 'external_weights', return_obj = F)
if(external_weights_flag == "T"){
  external_weights_pos_path <- read_param(config = opt$config_file, param = 'external_weights_pos_path', return_obj = F)
  external_weights <- gsub('.pos','',basename(external_weights_pos_path))
  weights<-c(weights, external_weights)
}

# Write out this list of SNP-weights as this might be useful elsewhere.
# Use writeLines so each panel is on its own line; downstream readers
# (e.g., read_twas_gsea_* in package_results_functions.R) parse this file with
# read.table(file)$V1, which only picks up the first space-separated token of
# each row -- so a single-line space-separated file would silently drop all
# but the first panel.
writeLines(as.character(weights), paste0(outdir,'/results/',opt$gwas,'/twas/list_of_weights.txt'))

# Read in gene names from pre-downloaded biomart data
biomart<-read.delim(paste0(resdir, '/data/biomart/biomart_genes_grch37.tsv'), stringsAsFactors=FALSE)
Genes<-biomart[,c('ensembl_gene_id','external_gene_name')]
Genes<-Genes[!duplicated(Genes),]

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
