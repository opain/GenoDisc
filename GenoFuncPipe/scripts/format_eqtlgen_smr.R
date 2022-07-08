#!/usr/bin/Rscript

suppressMessages(library("optparse"))
option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)

Genes<-Genes[!duplicated(Genes$ensembl_gene_id),]

smr_eqtlgen_files<-list.files(path=paste0('results/',opt$gwas,'/smr/eqtlgen/'), pattern=paste0(opt$gwas,'_smr_eqtlgen_chr'))
smr_eqtlgen_files<-smr_eqtlgen_files[grepl('.smr$', smr_eqtlgen_files)]

smr_eqtlgen<-NULL
for(i in smr_eqtlgen_files){
  smr_eqtlgen<-rbind(smr_eqtlgen, fread(paste0('results/',opt$gwas,'/smr/eqtlgen/',i)))
}

smr_eqtlgen<-smr_eqtlgen[!duplicated(smr_eqtlgen$probeID),]
smr_eqtlgen$ensembl_gene_id<-gsub('\\..*','',smr_eqtlgen$probeID)

smr_eqtlgen<-merge(smr_eqtlgen,Genes, by='ensembl_gene_id', all.x=T)

smr_eqtlgen$external_gene_name[smr_eqtlgen$external_gene_name == '']<-NA

fwrite(smr_eqtlgen, paste0('results/',opt$gwas,'/smr/eqtlgen/',opt$gwas,'_smr_eqtlgen_GW.txt.gz'), quote=F, sep=' ', na='NA')

