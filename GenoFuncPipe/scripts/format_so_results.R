#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--twas", action="store", default=NA, type='character',
              help="Name of GWAs/TWAS [required]"),
  make_option("--lincs_siginfo_path", action="store", default=NA, type='character',
              help="Path to LINCS signiture information data [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

# Read in the results
res<-fread(paste0('results/',opt$twas,'/twas/cmap/So_res.csv'))

# Subset set ID, rank and p
res<-res[,!grepl('perm',names(res)), with=F]

# Read in sig_info
sig_info<-fread(opt$lincs_siginfo_path)
sig_info$sig_id_2<-gsub('-','.',sig_info$sig_id)
sig_info$sig_id_2<-gsub(':','.',sig_info$sig_id_2)
sig_info$sig_id_2<-gsub('\\+','.',sig_info$sig_id_2)

sig_info<-sig_info[(sig_info$sig_id_2 %in% res$Drug),]

# Merge
res<-merge(res, sig_info[,c('sig_id_2','sig_id','pert_mfc_id','cmap_name','pert_idose','pert_itime','cell_iname'),with=F], by.x='Drug', by.y='sig_id_2', all.x=T)
res<-res[,c('sig_id','pert_mfc_id','cmap_name','pert_idose','pert_itime','cell_iname','obs','p'), with=F]
names(res)<-c('sig_id','pert_mfc_id','cmap_name','pert_idose','pert_itime','cell_iname','avg.rank','p')
  
# Sort by p-value
res<-res[order(res$p),]

# FDR correct
res$p.fdr<-p.adjust(res$p, method='fdr')

write.csv(res, paste0('results/',opt$twas,'/twas/cmap/So_res_clean.csv'), row.names=F)

# Insert ATC codes
atc<-fread('resources/data/atc/atc_20220201.txt', sep='!')
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
  
    wil_cox_res<-wilcox.test(rank(res_atc$avg.rank) ~ class_bin, conf.int =T, alternative='greater')
    
    atc_enrich<-rbind(atc_enrich, data.frame(ATC=cat,
                                             Estimate=as.numeric(wil_cox_res$estimate),
                                             Class_Median=median(res_atc$avg.rank[class_bin == 1]),
                                             Non_Class_Median=median(res_atc$avg.rank[class_bin == 0]),
                                             P=wil_cox_res$p.value,
                                             N=sum(class_bin)))
  }
}

atc_labels<-atc[nchar(atc$Code) == 4,]
atc_enrich<-merge(atc_enrich, atc_labels, by.x='ATC', by.y='Code')
atc_enrich<-atc_enrich[order(atc_enrich$P),]

write.csv(atc_enrich, paste0('results/',opt$twas,'/twas/cmap/So_res_atc_res.csv'), row.names=F)

