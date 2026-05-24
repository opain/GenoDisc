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

Genes<-Genes[!duplicated(Genes$external_gene_name),]

# Read in the rosmap smr results
smr_rosmap_files<-list.files(path=paste0(outdir,'/results/',opt$gwas,'/smr/rosmap/'), pattern=paste0(opt$gwas,'_smr_rosmap_chr'))
smr_rosmap_files<-smr_rosmap_files[grepl('.smr$', smr_rosmap_files)]

smr_rosmap<-NULL
for(i in smr_rosmap_files){
  smr_rosmap<-rbind(smr_rosmap, fread(paste0(outdir,'/results/',opt$gwas,'/smr/rosmap/',i)))
}

# Split rows containing a string of gene names into seperate rows
smr_rosmap_one_gene<-smr_rosmap[!grepl(';', smr_rosmap$Gene),]
smr_rosmap_mult_gene<-smr_rosmap[grepl(';', smr_rosmap$Gene),]

smr_rosmap_mult_gene_split<-NULL
for(i in 1:nrow(smr_rosmap_mult_gene)){
  ids<-unlist(strsplit(smr_rosmap_mult_gene$Gene[i], ';'))
  for(j in ids){
    tmp<-smr_rosmap_mult_gene[i,]
    tmp$Gene<-j
    smr_rosmap_mult_gene_split<-rbind(smr_rosmap_mult_gene_split, tmp)
  }
}

smr_rosmap<-rbind(smr_rosmap_one_gene, smr_rosmap_mult_gene_split)

smr_rosmap<-smr_rosmap[!is.na(smr_rosmap$Gene),]
smr_rosmap<-smr_rosmap[!duplicated(smr_rosmap$Gene),]
smr_rosmap$external_gene_name<-smr_rosmap$Gene

smr_rosmap<-merge(smr_rosmap,Genes, by='external_gene_name', all.x=T)
smr_rosmap<-smr_rosmap[!duplicated(smr_rosmap$ensembl_gene_id),]

fwrite(smr_rosmap, paste0(outdir,'/results/',opt$gwas,'/smr/rosmap/',opt$gwas,'_smr_rosmap_GW.txt.gz'), quote=F, sep=' ', na='NA')

