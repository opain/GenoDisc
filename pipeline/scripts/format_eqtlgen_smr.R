#!/usr/bin/Rscript

suppressMessages(library("optparse"))
option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

biomart<-read.delim('resources/data/biomart/biomart_genes_grch37.tsv', stringsAsFactors=FALSE)
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

