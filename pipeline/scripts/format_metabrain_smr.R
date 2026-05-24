#!/usr/bin/Rscript

suppressMessages(library("optparse"))
option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'functions', 'utils_functions.R'))

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

biomart<-read.delim(paste0(resdir, '/data/biomart/biomart_genes_grch37.tsv'), stringsAsFactors=FALSE)
Genes<-biomart[,c('ensembl_gene_id','external_gene_name')]
Genes<-Genes[!duplicated(Genes),]

Genes<-Genes[!duplicated(Genes$ensembl_gene_id),]

# Read in config file
smr_expression_panel_metabrain_basalganglia_logical<-config[grepl('smr_expression_panel_metabrain_basalganglia:',config)] == "smr_expression_panel_metabrain_basalganglia: T"
smr_expression_panel_metabrain_cerebellum_logical<-config[grepl('smr_expression_panel_metabrain_cerebellum:',config)] == "smr_expression_panel_metabrain_cerebellum: T"
smr_expression_panel_metabrain_cortex_logical<-config[grepl('smr_expression_panel_metabrain_cortex:',config)] == "smr_expression_panel_metabrain_cortex: T"
smr_expression_panel_metabrain_hippocampus_logical<-config[grepl('smr_expression_panel_metabrain_hippocampus:',config)] == "smr_expression_panel_metabrain_hippocampus: T"
smr_expression_panel_metabrain_spinalcord_logical<-config[grepl('smr_expression_panel_metabrain_spinalcord:',config)] == "smr_expression_panel_metabrain_spinalcord: T"

metabrain_tissues<-c('Basalganglia','Cerebellum','Cortex','Hippocampus','Spinalcord')
  
smr_metabrain<-list()
for(tissues in metabrain_tissues[c(smr_expression_panel_metabrain_basalganglia_logical,
                                   smr_expression_panel_metabrain_cerebellum_logical,
                                   smr_expression_panel_metabrain_cortex_logical,
                                   smr_expression_panel_metabrain_hippocampus_logical,
                                   smr_expression_panel_metabrain_spinalcord_logical)]){
  
  smr_metabrain_files<-list.files(path=paste0(outdir,'/results/',opt$gwas,'/smr/metabrain/',tissues,'/'), pattern=paste0(opt$gwas,'_smr_metabrain_',tissues,'_chr'))
  smr_metabrain_files<-smr_metabrain_files[grepl('.smr$', smr_metabrain_files)]

  smr_metabrain_tissue<-NULL
  for(i in smr_metabrain_files){
    smr_metabrain_tissue<-rbind(smr_metabrain_tissue, fread(paste0(outdir,'/results/',opt$gwas,'/smr/metabrain/',tissues,'/',i)))
  }
  
  smr_metabrain_tissue<-smr_metabrain_tissue[!duplicated(smr_metabrain_tissue$probeID),]
  smr_metabrain_tissue$ensembl_gene_id<-gsub('\\..*','',smr_metabrain_tissue$probeID)
  
  smr_metabrain_tissue<-merge(smr_metabrain_tissue,Genes, by='ensembl_gene_id', all.x=T)
  
  smr_metabrain_tissue$PANEL<-tissues

  smr_metabrain[[tissues]]<-smr_metabrain_tissue

}

smr_metabrain_all<-do.call(rbind, smr_metabrain)
smr_metabrain_all$external_gene_name[smr_metabrain_all$external_gene_name == '']<-NA

fwrite(smr_metabrain_all, paste0(outdir,'/results/',opt$gwas,'/smr/metabrain/',opt$gwas,'_smr_metabrain_GW.txt.gz'), quote=F, sep=' ', na='NA')

