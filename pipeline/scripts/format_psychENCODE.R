#!/usr/bin/Rscript
library(optparse)

option_list = list(
  make_option("--resdir", type="character", default="resources")
)
opt = parse_args(OptionParser(option_list=option_list))

# Create a list of ensemble IDs
IDs<-list.files(paste0(opt$resdir, '/data/fusion_snp_weights/psychencode/psychencode'))
IDs<-IDs[grepl('.wgt.RDat', IDs)]
IDs<-gsub('.wgt.RDat','',IDs)

biomart<-read.delim(paste0(opt$resdir, '/data/biomart/biomart_genes_grch37.tsv'), stringsAsFactors=FALSE)
Genes<-biomart[,c('ensembl_gene_id','chromosome_name','start_position','end_position')]
Genes<-Genes[!duplicated(Genes),]
Genes<-Genes[(Genes$ensembl_gene_id %in% IDs),]
Genes$chromosome_name<-as.numeric(Genes$chromosome_name)
Genes<-Genes[,c("chromosome_name","start_position","end_position","ensembl_gene_id")]
Genes<-Genes[order(Genes$chromosome_name, Genes$start_position),]

write.table(Genes, paste0(opt$resdir, '/data/fusion_snp_weights/psychencode.coord'), col.names=T, row.names=F, quote=F)

system(paste0('Rscript ', opt$resdir, '/software/Calculating-FUSION-TWAS-weights-pipeline/OP_packaging_fusion_weights.R --RDat_dir ', opt$resdir, '/data/fusion_snp_weights/psychencode/psychencode --coordinate_file ', opt$resdir, '/data/fusion_snp_weights/psychencode.coord --output_name psychencode --output_dir ', opt$resdir, '/data/fusion_snp_weights/psychencode/psychencode_new'))

# Delete original SNP-weights and rename new folder
system(paste0('rm -r ', opt$resdir, '/data/fusion_snp_weights/psychencode/psychencode'))
system(paste0('rm ', opt$resdir, '/data/fusion_snp_weights/psychencode.coord'))
system(paste0('mv ', opt$resdir, '/data/fusion_snp_weights/psychencode/psychencode_new/* ', opt$resdir, '/data/fusion_snp_weights/psychencode/'))
system(paste0('rm -r ', opt$resdir, '/data/fusion_snp_weights/psychencode/psychencode_new'))

pos<-read.table(paste0(opt$resdir, '/data/fusion_snp_weights/psychencode/psychencode.pos'), header=T)
pos$PANEL<-'psychencode'
pos$N<-1321
pos<-pos[,c('PANEL', 'WGT', 'ID', 'CHR', 'P0', 'P1', 'N')]

write.table(pos, paste0(opt$resdir, '/data/fusion_snp_weights/psychencode/psychencode.pos'), col.names=T, row.names=F, quote=F)

# Update SNP IDs to be RSIDs
library(data.table)

for(i in 1:22){
  print(i)
  pos_i<-pos[pos$CHR == i,]
  ref_i<-fread(paste0(opt$resdir, '/data/1kg/1KG.Phase3.EUR.MAF_001.chr',i,'.bim'))
  ref_i$ID<-paste0(ref_i$V1,':',ref_i$V4)
  names(ref_i)[2]<-'RSID'
  
  print(dim(pos_i)[1])
  
  for(k in 1:dim(pos_i)[1]){
    print(k)
    load(paste0(opt$resdir, '/data/fusion_snp_weights/psychencode/',pos_i$WGT[k]))
    
    ref_i_k<-ref_i[ref_i$V4 > (pos_i$P0[k] - 5e6) & ref_i$V4 < (pos_i$P1[k] + 5e6),]
    
    snps<-data.table(snps)
    snps_2<-merge(snps, ref_i_k[,c('RSID','ID')], by.x='V2', by.y='ID')
    snps_2<-snps_2[match(intersect(snps$V2,snps_2$V2), snps_2$V2),]
    snps_2<-snps_2[,c('V1','RSID','V3','V4','V5','V6')]
    names(snps_2)[2]<-'V2'
    
    snps<-data.frame(snps_2)
    rm(snps_2)
    
    wgt.matrix<-data.frame(wgt.matrix)
    wgt.matrix$ID<-row.names(wgt.matrix)
    wgt.matrix_2<-merge(wgt.matrix, ref_i_k[,c('RSID','ID')], by.x='ID', by.y='ID')
    wgt.matrix_2<-wgt.matrix_2[match(intersect(wgt.matrix$ID,wgt.matrix_2$ID), wgt.matrix_2$ID),]
    wgt.matrix_2$ID<-NULL
    row.names(wgt.matrix_2)<-wgt.matrix_2$RSID
    wgt.matrix_2$RSID<-NULL
    wgt.matrix_2<-as.matrix(wgt.matrix_2)
    
    wgt.matrix<-wgt.matrix_2
    rm(wgt.matrix_2)
    
    save(wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot, file = paste0(opt$resdir, '/data/fusion_snp_weights/psychencode/',pos_i$WGT[k]))
  }
}

file.create(paste0(opt$resdir, '/data/format_psychencode.done'))
