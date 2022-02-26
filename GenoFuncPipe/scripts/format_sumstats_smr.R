#!/usr/bin/Rscript

library(data.table)

ss<-fread('/users/k1806347/brc_scratch/Data/GWAS_sumstats/ALS/GWAS/ALS_sumstats_EUR_only.txt.gz')

names(ss)<-c('SNP','CHR','BP','A1','A2','FREQ','BETA','SE','P','N')
ss$A1<-toupper(ss$A1)
ss$A2<-toupper(ss$A2)

ss<-ss[,c('SNP','A1','A2','FREQ','BETA','SE','P','N'),with=F]
names(ss)<-c('SNP','A1','A2','freq','b','se','p','N')

fwrite(ss, '/users/k1806347/brc_scratch/Data/GWAS_sumstats/ALS/munged/ALS_sumstats_EUR_only.txt.cojo', sep=' ', na='NA', quote=F)

