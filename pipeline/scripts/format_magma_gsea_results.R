#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

# Read in MAGMA gene set results
res_gs<-fread(cmd=paste0("grep -v '^#' ",outdir,"/results/",opt$gwas,'/magma/magma_drug_targetor.gsa.out'))

# Remove gene sets with <5 genes present
#res_gs<-res_gs[res_gs$NGENES >= 5,]

res_gs$ATC<-gsub('ATC:','',gsub('\\|.*','',res_gs$FULL_NAME))
res_gs$NAME<-tolower(gsub('\\|.*','',gsub('.*NAME:','',res_gs$FULL_NAME)))

res_gs<-res_gs[,c('NAME','NGENES','BETA','SE','P','ATC'), with=F]
 
write.csv(res_gs, paste0(outdir,'/results/',opt$gwas,'/magma/magma_drug_targetor.clean.csv'), row.names=F, quote=T)

# Perform enrichment analysis of ATC codes
res_gs_enrich<-NULL
for(i in 1:nrow(res_gs)){
  res_gs_enrich<-rbind(res_gs_enrich, data.frame( ATC=unlist(strsplit(res_gs$ATC[i],',')),
                                                  P=res_gs$P[i]))
}

res_gs_enrich$atc_cat<-substr(res_gs_enrich$ATC, 1, 4) 

# Test for enrichment for each ATC category
atc_enrich<-NULL
for(cat in unique(res_gs_enrich$atc_cat)){
  class_bin<-rep(0, nrow(res_gs_enrich))
  class_bin[res_gs_enrich$atc_cat == cat]<-1
  
  if(sum(class_bin == 1) > 5){
    
    wil_cox_res<-wilcox.test(rank(res_gs_enrich$P) ~ class_bin, conf.int =T, alternative='greater')
    
    atc_enrich<-rbind(atc_enrich, data.frame(ATC=cat,
                                             Estimate=as.numeric(wil_cox_res$estimate),
                                             Class_Median=median(res_gs_enrich$P[class_bin == 1]),
                                             Non_Class_Median=median(res_gs_enrich$P[class_bin == 0]),
                                             P=wil_cox_res$p.value,
                                             N=sum(class_bin)))
  }
}

atc<-fread('resources/data/atc/atc_20220201.txt', sep='!')
names(atc)<-c('Code','Name')
atc$Name<-tolower(atc$Name)

atc_labels<-atc[nchar(atc$Code) == 4,]
atc_enrich<-merge(atc_enrich, atc_labels, by.x='ATC', by.y='Code')
atc_enrich<-atc_enrich[order(atc_enrich$P),]

write.csv(atc_enrich, paste0(outdir,'/results/',opt$gwas,'/magma/magma_drug_targetor_atc_res.csv'), row.names=F)

# Test for enrichment for each level 4 ATC category
res_gs_enrich$atc_cat_2<-substr(res_gs_enrich$ATC, 1, 5) 
atc_enrich_2<-NULL
for(cat in unique(res_gs_enrich$atc_cat_2)){
  class_bin<-rep(0, nrow(res_gs_enrich))
  class_bin[res_gs_enrich$atc_cat_2 == cat]<-1
  
  if(sum(class_bin == 1) > 1){
    
    wil_cox_res<-wilcox.test(rank(res_gs_enrich$P) ~ class_bin, conf.int =T, alternative='greater')
    
    atc_enrich_2<-rbind(atc_enrich_2, data.frame(ATC=cat,
                                             Estimate=as.numeric(wil_cox_res$estimate),
                                             Class_Median=median(res_gs_enrich$P[class_bin == 1]),
                                             Non_Class_Median=median(res_gs_enrich$P[class_bin == 0]),
                                             P=wil_cox_res$p.value,
                                             N=sum(class_bin)))
  }
}

names(atc)<-c('Code','Name')
atc$Name<-tolower(atc$Name)

atc_labels<-atc[nchar(atc$Code) == 5,]
atc_enrich_2<-merge(atc_enrich_2, atc_labels, by.x='ATC', by.y='Code')
atc_enrich_2<-atc_enrich_2[order(atc_enrich_2$P),]

write.csv(atc_enrich_2, paste0(outdir,'/results/',opt$gwas,'/magma/magma_drug_targetor_atc_res_level4.csv'), row.names=F)


