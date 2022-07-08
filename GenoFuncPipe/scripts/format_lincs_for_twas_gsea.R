#!/usr/bin/Rscript
library(data.table)

# Read in cmap data and recombine batches
cmap<-NULL
for(i in 1:10){
  if(i == 1){
    cmap<-readRDS(paste0("resources/data/lincs/lincs_core_subset_batch_",i,".rds"))
  } else {
    tmp<-readRDS(paste0("resources/data/lincs/lincs_core_subset_batch_",i,".rds"))
    cmap<-cbind(cmap, tmp[,-1])
  }
}

names(cmap)[1]<-'ID'

fwrite(cmap, 'resources/data/lincs/lincs_core_subset.txt.gz', na='NA', sep=' ')

