#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--panel", action="store", default=NA, type='character',
              help="Panel ID [required]"),
  make_option("--resdir", type="character", default="resources")
)

opt = parse_args(OptionParser(option_list=option_list))

# Read in gene locations from build GRCh37
biomart<-read.delim(paste0(opt$resdir, '/data/biomart/biomart_genes_grch37.tsv'), stringsAsFactors=FALSE)
Genes<-biomart[,c('ensembl_gene_id','chromosome_name','start_position','end_position')]
Genes<-Genes[!duplicated(Genes),]

library(data.table)

pos<-fread(paste0(opt$resdir, '/data/fusion_snp_weights/',opt$panel,'/',opt$panel,'.pos'))
pos$PANEL<-opt$panel
pos$WGT<-gsub('.*/',paste0(opt$panel,'/'), pos$WGT)
pos$ensembl_gene_id<-gsub('\\..*','',pos$ID)
pos<-merge(pos, Genes, by='ensembl_gene_id')
pos<-pos[!duplicated(pos$ID),]
pos<-pos[,c('PANEL','WGT','ID','CHR','start_position','end_position','N'), with=F]
names(pos)<-c('PANEL','WGT','ID','CHR','P0','P1','N')
pos$ID<-gsub('\\..*','',pos$ID)
fwrite(pos, paste0(opt$resdir, '/data/fusion_snp_weights/',opt$panel,'/',opt$panel,'.pos'), quote=F, sep=' ', na='NA')

# Delete WGT files and hsq files for WGTs not in pos, i.e with non-sig h2
dir.create(paste0(opt$resdir, '/data/fusion_snp_weights/',opt$panel,'/',opt$panel,'_new'))
for(i in 1:nrow(pos)){
  system(paste0('mv ', opt$resdir, '/data/fusion_snp_weights/',opt$panel,'/',pos$WGT[i],' ', opt$resdir, '/data/fusion_snp_weights/',opt$panel,'/',opt$panel,'_new/'))
}

system(paste0('rm -r ', opt$resdir, '/data/fusion_snp_weights/',opt$panel,'/',opt$panel))
system(paste0('mv ', opt$resdir, '/data/fusion_snp_weights/',opt$panel,'/',opt$panel,'_new ', opt$resdir, '/data/fusion_snp_weights/',opt$panel,'/',opt$panel))


