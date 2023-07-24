#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--twas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character',
              help="Panel [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(dplyr)
library(Hmisc)
library(data.table)

# Read in the results of step1 for all batches
batches<-list()
for(i in 1:10){
  batches[[paste0('batch_',i)]]<-readRDS(paste0('results/',opt$twas,'/twas/cmap/res_',opt$panel,'_batch_',i,'.RDS'))
}

# Combine batches for each the real data (obs) and each permutation, and calculate average ranks
average_ranks<-list()
for(i in c('obs',paste0('perm_',1:100))){
  drugs<-NULL
  cor.pearson<-NULL
  cor.spearman<-NULL
  ks.signed<-NULL
  extreme.cor.pearson<-NULL
  extreme.cor.spearman<-NULL

  for(k in names(batches)){
    drugs<-c(drugs,batches[[k]][[i]]$drugs)
    cor.pearson<-c(cor.pearson,batches[[k]][[i]]$cor.pearson)
    cor.spearman<-c(cor.spearman,batches[[k]][[i]]$cor.spearman)
    ks.signed<-rbind(ks.signed,batches[[k]][[i]]$ks.signed)
    extreme.cor.pearson<-rbind(extreme.cor.pearson,batches[[k]][[i]]$extreme.cor.pearson)
    extreme.cor.spearman<-rbind(extreme.cor.spearman,batches[[k]][[i]]$extreme.cor.spearman)
  }
  
  no.thres.N<-ncol(extreme.cor.pearson)
  nodrugs<-length(drugs)
  
  # Now insert average rank calculation
  rank.ks = matrix(nrow=nodrugs, ncol=no.thres.N)
  rank.spearman = matrix(nrow=nodrugs, ncol=no.thres.N)
  rank.pearson = matrix(nrow=nodrugs, ncol=no.thres.N)
  
  for (j in 1:no.thres.N){
    rank.ks[,j]=rank(ks.signed[,j])
    rank.spearman[,j]=rank(extreme.cor.spearman[,j])
    rank.pearson[,j]=rank(extreme.cor.pearson[,j])
  }
  
  mean.rank.ks = rowMeans(rank.ks)
  mean.rank.spearman = rowMeans(rank.spearman)
  mean.rank.pearson = rowMeans(rank.pearson)
  
  five.method.rank = cbind(rank(cor.pearson),
                           rank(cor.spearman), 
                           rank(mean.rank.ks), 
                           rank(mean.rank.spearman), 
                           rank(mean.rank.pearson) )
  
  res.5method.avg = data.frame( drugs, rowMeans(five.method.rank) ) 
  colnames(res.5method.avg)   <- c("Drug",i)    
  
  average_ranks[[i]]<-res.5method.avg
}

# Then insert permutation p.value calculation
# Note. Average ranks across all drugs are used to calculate p
average_ranks_tab<-Reduce(function(...) merge(..., all=T, by='Drug'), average_ranks)
all_perm_ranks<-as.vector(as.matrix(average_ranks_tab[,grepl('perm_',names(average_ranks_tab))]))
  
for(i in 1:nrow(average_ranks_tab)){
  average_ranks_tab$p[i]<-  sum(average_ranks_tab$obs[i] > all_perm_ranks) / length(all_perm_ranks)
}

write.csv(average_ranks_tab,  paste0('results/',opt$twas,'/twas/cmap/So_res_',opt$panel,'.csv'), row.names=F)
