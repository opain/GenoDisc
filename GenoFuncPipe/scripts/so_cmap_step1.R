#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--twas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character',
              help="Panel [required]"),
  make_option("--batch", action="store", default=NA, type='character',
              help="batch number [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(dplyr)
library(Hmisc)
library(data.table)

# Read in the lincs data
cmap<-readRDS(paste0('resources/data/lincs/lincs_core_subset_batch_',opt$batch,'.rds'))

# Read in TWAS results
twas<-fread(paste0('results/',opt$twas,'/twas/',opt$twas,'_twas_',opt$panel,'_GW_clean.txt.gz'))

# Subset relevent data from TWAS
twas<-twas[,c('ensembl_gene_id','TWAS.Z','TWAS.P'),with=T]
names(twas)[names(twas) == 'ensembl_gene_id']<-'ID'
twas<-twas[complete.cases(twas),]
twas<-twas[abs(twas$TWAS.Z) != Inf,]
twas<-twas[order(twas$TWAS.P),]

# Combine the lincs data and twas results
finalres<-merge(twas, cmap, by.x='ID', by.y='ensembl_gene_id')

rm(twas, cmap)

#################
# Run So et al. method
#################

thres.N.vector = c(50,100,250,500)
noperm = 100 # The true number of permuations is nodrugs*noperm

nodrugs = ncol(finalres)-3
no.thres.N = length(thres.N.vector)

cor.spearman=NULL
cor.pearson=NULL
p.spearman = NULL
p.pearson = NULL
ks.ES = NULL
ks.p = NULL 

ks.signed = matrix(ncol=no.thres.N,nrow=nodrugs)
extreme.cor.spearman = matrix(ncol=no.thres.N,nrow=nodrugs)
extreme.cor.pearson = matrix(ncol=no.thres.N,nrow=nodrugs)
connect.score = vector()
finalres = data.frame(finalres)

#******************************************************
#         comparing expression in drugs vs diseases
#******************************************************
for (i in 1:nodrugs) {
  print(i)
  disease.zscore = finalres$TWAS.Z
  spearman.obj = cor.test(finalres[,i+3],disease.zscore, method="spearman",use="na.or.complete")
  cor.spearman[i] = spearman.obj$estimate
  pearson.obj  = cor.test(finalres[,i+3],disease.zscore, method = "pearson",use="na.or.complete")
  cor.pearson[i]= pearson.obj$estimate
  p.spearman[i] = spearman.obj$p.value
  p.pearson[i]= pearson.obj$p.value
  
  nogenes = nrow(finalres)
  pos.zscore = disease.zscore[disease.zscore>=0]
  no.pos.genes = length(pos.zscore)
  neg.zscore = disease.zscore[disease.zscore<0]
  no.neg.genes = length(neg.zscore)
  
  for (j in 1:no.thres.N) {
    thres.N = thres.N.vector[j]
    
    #****************************************
    #   KS method
    #****************************************
    rank.pos = no.pos.genes+1-rank(pos.zscore)  
    up.thres = pos.zscore[rank.pos==thres.N] 
    rank.neg = rank(neg.zscore)
    down.thres = neg.zscore[rank.neg==thres.N]
    
    ind.upreg = order(disease.zscore,decreasing=T)[1:thres.N]
    
    geneset2 = rank(-finalres[,i+3])[ind.upreg] 
    geneset2 = sort(geneset2)
    a.up = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
    b.up = max(   geneset2/nogenes -  (1:thres.N-1)/thres.N )                                           
    ks.test.obj.up = ifelse(a.up>b.up, a.up , -b.up)
    
    ind.downreg= order(disease.zscore,decreasing=F)[1:thres.N]
    
    geneset2 = rank(-finalres[,i+3])[ind.downreg]
    geneset2 = sort(geneset2)
    
    a.down = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
    b.down = max(geneset2/nogenes - (1:thres.N-1)/thres.N )                                           
    ks.test.obj.down = ifelse(a.down>b.down, a.down , -b.down)
    
    
    ks.signed[i,j]=0
    if (  sign(ks.test.obj.up)!=sign(ks.test.obj.down)  ) 
    {
      ks.signed[i,j] = ks.test.obj.up - ks.test.obj.down
    }
    
    #****************************************
    #   spearman and Pearson correlations (either use all observations or only compare the most extreme observations)
    #****************************************
    spearman.obj.extr = cor.test( finalres[,i+3][c(ind.upreg,ind.downreg)], disease.zscore[c(ind.upreg,ind.downreg)] 
                                  ,method="spearman",use="na.or.complete"   )
    extreme.cor.spearman[i,j] =  spearman.obj.extr$estimate 
    pearson.obj.extr  = cor.test(finalres[,i+3][c(ind.upreg,ind.downreg)],disease.zscore[c(ind.upreg,ind.downreg)], method = "pearson",use="na.or.complete")
    extreme.cor.pearson[i,j]= pearson.obj.extr$estimate
    
  }
}

res_list_list<-list()
res_list<-list()
res_list[['drugs']]<-colnames(finalres)[4:ncol(finalres)]
res_list[['cor.pearson']]<-cor.pearson
res_list[['cor.spearman']]<-cor.spearman
res_list[['ks.signed']]<-ks.signed
res_list[['extreme.cor.pearson']]<-extreme.cor.pearson
res_list[['extreme.cor.spearman']]<-extreme.cor.spearman

res_list_list[['obs']]<-res_list

#*********************************************************************
# permutation to determine the significance of the avg. rank obtained (repeats the above code)
#*********************************************************************

for (r in 1:noperm) {
  print(r)
  disease.zscore = sample(finalres$TWAS.Z)
  
  cor.spearman=NULL
  cor.pearson=NULL
  
  ks.signed = matrix(ncol=no.thres.N,nrow=nodrugs)
  extreme.cor.spearman = matrix(ncol=no.thres.N,nrow=nodrugs)
  extreme.cor.pearson = matrix(ncol=no.thres.N,nrow=nodrugs)
  avgrank.perm <- matrix(nrow=nodrugs, ncol=noperm)
  
  for (i in 1:nodrugs) {
    spearman.obj = cor.test(finalres[,i+3],disease.zscore, method="spearman",use="na.or.complete")
    cor.spearman[i] = spearman.obj$estimate
    pearson.obj  = cor.test(finalres[,i+3],disease.zscore, method = "pearson",use="na.or.complete")
    cor.pearson[i]= pearson.obj$estimate
    
    nogenes = nrow(finalres)
    
    pos.zscore = disease.zscore[disease.zscore>=0]
    no.pos.genes = length(pos.zscore)
    neg.zscore = disease.zscore[disease.zscore<0]
    no.neg.genes = length(neg.zscore)
    
    
    for (j in 1:no.thres.N) {
      thres.N = thres.N.vector[j]
      ind.upreg = order(disease.zscore,decreasing=T)[1:thres.N]
      
      #****************************************
      #   KS method
      #******************************************
      geneset2 = rank(-finalres[,i+3])[ind.upreg] 
      geneset2 = sort(geneset2) 
      
      a.up = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
      b.up = max(   geneset2/nogenes -  (1:thres.N-1)/thres.N )                                           
      ks.test.obj.up = ifelse(a.up>b.up, a.up , -b.up)
      
      
      ind.downreg= order(disease.zscore,decreasing=F)[1:thres.N]
      
      geneset2 = rank(-finalres[,i+3])[ind.downreg]
      geneset2 = sort(geneset2)
      a.down = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
      b.down = max(geneset2/nogenes - (1:thres.N-1)/thres.N )                                           
      ks.test.obj.down = ifelse(a.down>b.down, a.down , -b.down)
      
      
      ks.signed[i,j]=0
      if (  sign(ks.test.obj.up)!=sign(ks.test.obj.down)  ) 
      {
        ks.signed[i,j] = ks.test.obj.up - ks.test.obj.down
        
      }
      
      #****************************************
      #   spearman and Pearson correlations (either use all observations or only compare the most extreme observations)
      #****************************************
      spearman.obj.extr = cor.test( finalres[,i+3][c(ind.upreg,ind.downreg)], disease.zscore[c(ind.upreg,ind.downreg)] 
                                    ,method="spearman",use="na.or.complete"   )
      extreme.cor.spearman[i,j] =  spearman.obj.extr$estimate 
      pearson.obj.extr  = cor.test(finalres[,i+3][c(ind.upreg,ind.downreg)],disease.zscore[c(ind.upreg,ind.downreg)], method = "pearson",use="na.or.complete")
      extreme.cor.pearson[i,j]= pearson.obj.extr$estimate
      
    }
    
  } #end of looping over all drugs
  
  perm_res_list<-list()
  perm_res_list[['drugs']]<-colnames(finalres)[4:ncol(finalres)]
  perm_res_list[['cor.pearson']]<-cor.pearson
  perm_res_list[['cor.spearman']]<-cor.spearman
  perm_res_list[['ks.signed']]<-ks.signed
  perm_res_list[['extreme.cor.pearson']]<-extreme.cor.pearson
  perm_res_list[['extreme.cor.spearman']]<-extreme.cor.spearman
  
  res_list_list[[paste0('perm_',r)]]<-perm_res_list

}  # end of permutation loop

dir.create(paste0('results/',opt$twas,'/twas/cmap'))

saveRDS(res_list_list, paste0('results/',opt$twas,'/twas/cmap/res_',opt$panel,'_batch_',opt$batch,'.RDS'))


