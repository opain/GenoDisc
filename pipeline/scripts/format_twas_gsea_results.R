#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--twas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character',
              help="PANEL [required]"),
  make_option("--lincs_siginfo_path", action="store", default=NA, type='character',
              help="Path to LINCS sig_info file [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

# Read in TWAS-GSEA results
res<-fread(paste0(outdir,'/results/',opt$twas,'/twas/cmap/twas_gsea_cmap_',opt$panel,'.competitive.txt'))

# Flip one-sided hypothesis to be testing for a negative correlation
res$P<-2*pnorm(-abs(res$T))
res$P.CORR<-p.adjust(res$P, method='fdr')

# Read in LINCS sig info file
sig_info<-fread(opt$lincs_siginfo_path)
sig_info$sig_id_2<-gsub('-','.',sig_info$sig_id)
sig_info$sig_id_2<-gsub(':','.',sig_info$sig_id_2)
sig_info$sig_id_2<-gsub('\\+','.',sig_info$sig_id_2)
sig_info<-sig_info[(sig_info$sig_id_2 %in% res$GeneSet),]

# Merge sig_info file and results
res<-merge(res, sig_info[,c('sig_id_2','sig_id','pert_mfc_id','cmap_name','pert_idose','pert_itime','cell_iname'),with=F], by.x='GeneSet', by.y='sig_id_2', all.x=T)
res<-res[,c('sig_id','pert_mfc_id','cmap_name','pert_idose','pert_itime','cell_iname','Estimate','SE','T','P','P.CORR'), with=F]

# Sort by p-value
res<-res[order(res$P),]

write.csv(res, paste0(outdir,'/results/',opt$twas,'/twas/cmap/twas_gsea_cmap_',opt$panel,'.competitive.clean.csv'), row.names=F)

# Insert ATC codes
atc<-fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep='!')
names(atc)<-c('Code','Name')
atc$Name<-tolower(atc$Name)
res$cmap_name<-tolower(res$cmap_name)

res_atc<-merge(res, atc, by.x='cmap_name', by.y='Name')
res_atc$atc_cat<-substr(res_atc$Code, 1, 4) 

# Test for enrichment for each ATC catagory
atc_enrich<-NULL
for(cat in unique(res_atc$atc_cat)){
  class_bin<-rep(0, nrow(res_atc))
  class_bin[res_atc$atc_cat == cat]<-1
  
  if(sum(class_bin == 1) > 5){
    
    wil_cox_res<-wilcox.test(rank(res_atc$Estimate) ~ class_bin, conf.int =T)
    
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

write.csv(atc_enrich, paste0(outdir,'/results/',opt$twas,'/twas/cmap/twas_gsea_',opt$panel,'_res_atc_res.csv'), row.names=F)

