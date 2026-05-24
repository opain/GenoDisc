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
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

biomart<-read.delim(paste0(resdir, '/data/biomart/biomart_genes_grch37.tsv'), stringsAsFactors=FALSE)
Genes<-biomart[,c('ensembl_gene_id','external_gene_name')]
Genes<-Genes[!duplicated(Genes),]

Genes<-Genes[!duplicated(Genes$ensembl_gene_id),]

smr_eqtlgen_files<-list.files(path=paste0(outdir,'/results/',opt$gwas,'/smr/eqtlgen/'), pattern=paste0(opt$gwas,'_smr_eqtlgen_chr'))
smr_eqtlgen_files<-smr_eqtlgen_files[grepl('.smr$', smr_eqtlgen_files)]

smr_eqtlgen<-NULL
for(i in smr_eqtlgen_files){
  smr_eqtlgen<-rbind(smr_eqtlgen, fread(paste0(outdir,'/results/',opt$gwas,'/smr/eqtlgen/',i)))
}

smr_eqtlgen<-smr_eqtlgen[!duplicated(smr_eqtlgen$probeID),]
smr_eqtlgen$ensembl_gene_id<-gsub('\\..*','',smr_eqtlgen$probeID)

smr_eqtlgen<-merge(smr_eqtlgen,Genes, by='ensembl_gene_id', all.x=T)

smr_eqtlgen$external_gene_name[smr_eqtlgen$external_gene_name == '']<-NA

fwrite(smr_eqtlgen, paste0(outdir,'/results/',opt$gwas,'/smr/eqtlgen/',opt$gwas,'_smr_eqtlgen_GW.txt.gz'), quote=F, sep=' ', na='NA')

