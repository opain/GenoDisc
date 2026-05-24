#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'functions', 'utils_functions.R'))

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

chunks<-fread(paste0(outdir,'/results/',opt$gwas,'/gcsc/drugtargetor_gcsc_sets.nset.txt'))$x
      
res<-NULL
for(i in 1:length(chunks)){
  res<-rbind(res, fread(paste0(outdir,'/results/',opt$gwas,'/gcsc/drugtargetor/',i,'/GCSCresults.txt'), fill=T, sep=' '))
}

res<-res[grepl('enrichment:',res$Parameter),]
res<-res[order(res$`P-value`),]

res$Drug<-gsub('_CID.*','',gsub('.*NAME_','',res$Parameter))
res$Drug<-tolower(res$Drug)

res$ATC<-gsub('ATC_','',gsub('_NAME.*','',res$Parameter))

res<-res[,c('Drug','Value','Standard_error','P-value','ATC'), with=F]
names(res)<-c('Drug','Enrichment','SE','P','ATC')

res$P.FDR<-p.adjust(res$P, method='fdr')

fwrite(res, paste0(outdir,'/results/',opt$gwas,'/gcsc/',opt$gwas,'_drugtargetor_gcsc_res.txt'), sep=' ', quote=F, na='NA')

# Perform enrichment analysis of ATC codes
res_enrich<-NULL
for(i in 1:nrow(res)){
  res_enrich<-rbind(res_enrich, data.frame( ATC=unlist(strsplit(res$ATC[i],'_')),
                                                  Enrichment=res$Enrichment[i],
                                                  P=res$P[i]))
}

res_enrich$atc_cat<-substr(res_enrich$ATC, 1, 4) 

# Test for enrichment for each ATC category
atc_enrich<-NULL
for(cat in unique(res_enrich$atc_cat)){
  class_bin<-rep(0, nrow(res_enrich))
  class_bin[res_enrich$atc_cat == cat]<-1
  
  if(sum(class_bin == 1) > 5){
    
    # Rank genes by an enrichment signed -log10(p-value). 
    wil_cox_res<-wilcox.test(rank(sign(res_enrich$Enrichment) * -log10(res_enrich$P)) ~ class_bin, conf.int =T, alternative='less')
    
    atc_enrich<-rbind(atc_enrich, data.frame(ATC=cat,
                                             Estimate=as.numeric(wil_cox_res$estimate),
                                             Class_Median=median(res_enrich$Enrichment[class_bin == 1]),
                                             Non_Class_Median=median(res_enrich$Enrichment[class_bin == 0]),
                                             P=wil_cox_res$p.value,
                                             N=sum(class_bin)))
  }
}

atc<-fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep='!')
names(atc)<-c('Code','Name')
atc$Name<-tolower(atc$Name)

atc_labels<-atc[nchar(atc$Code) == 4,]
atc_enrich<-merge(atc_enrich, atc_labels, by.x='ATC', by.y='Code')
atc_enrich<-atc_enrich[order(atc_enrich$P),]
atc_enrich$P.FDR<-p.adjust(atc_enrich$P, method='fdr')

fwrite(atc_enrich, paste0(outdir,'/results/',opt$gwas,'/gcsc/',opt$gwas,'_drugtargetor_gcsc_res_atc.csv'), sep=',', quote=T, na='NA')

