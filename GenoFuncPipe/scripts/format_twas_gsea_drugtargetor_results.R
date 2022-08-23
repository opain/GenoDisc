#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--twas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character',
              help="PANEL [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

# Read in TWAS-GSEA results
res<-fread(paste0('results/',opt$twas,'/twas/drugtargetor/twas_gsea_drugtargetor_',opt$panel,'.competitive.txt'))

# Flip one-sided hypothesis to be testing for a negative correlation
res$P<-2*pnorm(-abs(res$T))
res$P.CORR<-p.adjust(res$P, method='fdr')

res$ATC<-gsub('ATC.','',gsub('\\.NAME.*','',res$GeneSet))
res$NAME<-tolower(gsub('\\.CID\\..*','',gsub('.*NAME\\.','',res$GeneSet)))

# Sort by p-value
res<-res[order(res$P),]

write.csv(res, paste0('results/',opt$twas,'/twas/drugtargetor/twas_gsea_drugtargetor_',opt$panel,'.competitive.clean.csv'), row.names=F)

# Insert ATC codes
atc<-fread('resources/data/atc/atc_20220201.txt', sep='!')
names(atc)<-c('Code','Name')
atc$Name<-tolower(atc$Name)

res_atc<-merge(res, atc, by.x='NAME', by.y='Name')
res_atc$atc_cat<-substr(res_atc$Code, 1, 4) 

# Test for enrichment for each ATC catagory
atc_enrich<-NULL
for(cat in unique(res_atc$atc_cat)){
  class_bin<-rep(0, nrow(res_atc))
  class_bin[res_atc$atc_cat == cat]<-1
  
  if(sum(class_bin == 1) > 5){
    
    # Use wilcoxon test
    wil_cox_res<-wilcox.test(rank(res_atc$Estimate) ~ class_bin, conf.int =T, alternative='greater')
    
    atc_enrich<-rbind(atc_enrich, data.frame(ATC=cat,
                                             Estimate=as.numeric(wil_cox_res$estimate),
                                             Class_Median=median(res_atc$Estimate[class_bin == 1]),
                                             Non_Class_Median=median(res_atc$Estimate[class_bin == 0]),
                                             P=wil_cox_res$p.value,
                                             N=sum(class_bin)))
  }
}

atc_labels<-atc[nchar(atc$Code) == 4,]
atc_enrich<-merge(atc_enrich, atc_labels, by.x='ATC', by.y='Code')
atc_enrich<-atc_enrich[order(atc_enrich$P),]

write.csv(atc_enrich, paste0('results/',opt$twas,'/twas/drugtargetor/twas_gsea_',opt$panel,'_res_atc_res.csv'), row.names=F)

