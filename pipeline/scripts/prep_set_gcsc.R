#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--config", action="store", default=NA, type='character',
              help="config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(biomaRt)

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

# Read in database
pathways<-fread('resources/data/drug_targetor/wholedatabase_for_targetor')

# Read in GCSC gene universe
universe<-fread('resources/software/GCSC/gene_universe.txt', header=F)

gene_id<-data.frame(universe)
names(gene_id)<-'ID'

# Use BioMart to convert gene IDs
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
biomartCacheClear()
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)
Genes<-Genes[!duplicated(Genes$ensembl_gene_id),]

pathways_id<-merge(pathways, Genes, by.x='gene', by.y='external_gene_name')

drugs<-unique(pathways$atc)
for(i in 1:length(drugs)){
  drug_tmp<-pathways_id[pathways_id$atc == drugs[i],]
  
  gene_id[[drugs[i]]]<-ifelse(gene_id$ID %in% drug_tmp$ensembl_gene_id, 1, 0)
  
}

# Remove sets with less than 2 genes present in each gene set and in the TWAS
config<-readLines(opt$config)
gcsc_tissues<-config[grepl('^gcsc_tissues', config)]
gcsc_tissues<-unlist(strsplit(gsub('"','',gsub('\\]','',gsub('.*\\[','',gcsc_tissues))),','))

twas<-NULL
for(weight_i in gcsc_tissues){
  for(chr_i in 1:22){
    twas<-rbind(
      twas,
      fread(paste0(outdir,'/results/',opt$gwas,'/gcsc/twas/',weight_i,'/',opt$gwas,'_twas_',weight_i,'_chr',chr_i,'.dat')))
  }
}

twas$ensembl_id<-gsub('\\..*','',gsub('.*ENSG','ENSG',twas$FILE))
twas<-twas[!is.na(twas$TWAS.Z),]

gene_id_all<-gene_id[gene_id$ID %in% twas$ensembl_id,]
gene_id_all_2<-gene_id[,c(1, which(colSums(abs(gene_id_all[,-1])) >= 2)+1)]

gene_id_all_mat<-matrix(as.matrix(gene_id_all_2[,-1]), nrow=nrow(gene_id_all_2), ncol=ncol(gene_id_all_2)-1, dimnames=list(gene_id_all_2$ID,gsub(' ' ,'',gsub('[[:punct:]]','_',names(gene_id_all_2)[-1]))))
gene_id_all_mat<-t(gene_id_all_mat)

chunks<-split(1:nrow(gene_id_all_mat), ceiling(seq_along(1:nrow(gene_id_all_mat))/10))
for(chunks_i in 1:length(chunks)){
  write.csv(gene_id_all_mat[chunks[[chunks_i]],], paste0(outdir,'/results/',opt$gwas,'/gcsc/drugtargetor_gcsc_sets_',chunks_i,'.csv'), quote=F)
}

write.table(data.frame(x=1:length(chunks)), paste0(outdir,'/results/',opt$gwas,'/gcsc/drugtargetor_gcsc_sets.nset.txt'), col.names=T, row.names=F, quote=F)
