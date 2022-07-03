#!/usr/bin/Rscript

# Create a list of ensemble IDs
IDs<-list.files('resources/data/fusion_snp_weights/psychencode/psychencode')
IDs<-IDs[grepl('.wgt.RDat', IDs)]
IDs<-gsub('.wgt.RDat','',IDs)

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart = ensembl)
Genes<-Genes[(Genes$ensembl_gene_id %in% IDs),]
Genes$chromosome_name<-as.numeric(Genes$chromosome_name)
Genes<-Genes[,c("chromosome_name","start_position","end_position","ensembl_gene_id")]
Genes<-Genes[order(Genes$chromosome_name, Genes$start_position),]

write.table(Genes, 'resources/data/fusion_snp_weights/psychencode.coord', col.names=T, row.names=F, quote=F)

system('Rscript resources/software/Calculating-FUSION-TWAS-weights-pipeline/OP_packaging_fusion_weights.R --RDat_dir resources/data/fusion_snp_weights/psychencode/psychencode --coordinate_file resources/data/fusion_snp_weights/psychencode.coord --output_name psychencode --output_dir resources/data/fusion_snp_weights/psychencode/psychencode_new')

# Delete original SNP-weights and rename new folder
system('rm -r resources/data/fusion_snp_weights/psychencode/psychencode')
system('rm resources/data/fusion_snp_weights/psychencode.coord')
system('mv resources/data/fusion_snp_weights/psychencode/psychencode_new/* resources/data/fusion_snp_weights/psychencode/')
system('rm -r resources/data/fusion_snp_weights/psychencode/psychencode_new')

pos<-read.table('resources/data/fusion_snp_weights/psychencode/psychencode.pos', header=T)
pos$PANEL<-'psychencode'
pos$N<-1321
pos<-pos[,c('PANEL', 'WGT', 'ID', 'CHR', 'P0', 'P1', 'N')]

write.table(pos, 'resources/data/fusion_snp_weights/psychencode/psychencode.pos', col.names=T, row.names=F, quote=F)

# Update SNP IDs to be RSIDs
library(data.table)

for(i in 1:22){
  print(i)
  pos_i<-pos[pos$CHR == i,]
  ref_i<-fread(paste0('resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr',i,'.bim'))
  ref_i$ID<-paste0(ref_i$V1,':',ref_i$V4)
  names(ref_i)[2]<-'RSID'
  
  print(dim(pos_i)[1])
  
  for(k in 1:dim(pos_i)[1]){
    print(k)
    load(paste0('resources/data/fusion_snp_weights/psychencode/',pos_i$WGT[k]))
    
    ref_i_k<-ref_i[ref_i$V4 > (pos_i$P0[k] - 5e6) & ref_i$V4 < (pos_i$P1[k] + 5e6),]
    
    snps<-data.table(snps)
    snps_2<-merge(snps, ref_i_k[,c('RSID','ID')], by.x='V2', by.y='ID')
    snps_2<-snps_2[match(intersect(snps$V2,snps_2$V2), snps_2$V2),]
    snps_2<-snps_2[,c('V1','RSID','V3','V4','V5','V6')]
    names(snps_2)[2]<-'V2'
    
    snps<-snps_2
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
    
    save(wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot, file = paste0('resources/data/fusion_snp_weights/psychencode/',pos_i$WGT[k]))
  }
}

file.create('resources/data/format_psychencode.done')
