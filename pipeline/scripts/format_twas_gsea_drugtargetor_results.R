#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--twas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--panel", action="store", default=NA, type='character',
              help="PANEL [required]"),
  make_option("--mode", action="store", default='directional', type='character',
              help="Either 'directional' (default) or 'nondirectional'"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

if(!(opt$mode %in% c('directional','nondirectional'))) stop("--mode must be 'directional' or 'nondirectional'")
suffix <- if(opt$mode == 'nondirectional') '_nondir' else ''

library(data.table)
source('scripts/functions/utils_functions.R')

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

# Read in TWAS-GSEA results
res<-fread(paste0(outdir,'/results/',opt$twas,'/twas/drugtargetor/twas_gsea_drugtargetor',suffix,'_',opt$panel,'.competitive.txt'), sep = ' ')

if(opt$mode == 'directional'){
  # Both signs of T are interpretable (drug mimics vs reverses disease signature),
  # so use a two-sided p-value.
  res$P <- 2*pnorm(-abs(res$T))
  wilcox_alt <- 'two.sided'
} else {
  # TWAS-GSEA-fast already emits a one-sided right-tail P (probit(1-P) outcome),
  # which is what we want: only positive enrichment of drug-target genes is meaningful.
  wilcox_alt <- 'greater'
}

res$P.CORR<-p.adjust(res$P, method='fdr')

res$GeneSet <- gsub('[[:punct:]]', '.', res$GeneSet)

res$ATC<-gsub('ATC.','',gsub('\\.NAME.*','',res$GeneSet))
res$NAME<-tolower(gsub('\\.CID\\..*','',gsub('.*NAME\\.','',res$GeneSet)))

# Sort by p-value
res<-res[order(res$P),]

write.csv(res, paste0(outdir,'/results/',opt$twas,'/twas/drugtargetor/twas_gsea_drugtargetor',suffix,'_',opt$panel,'.competitive.clean.csv'), row.names=F)

# Insert ATC codes
atc<-fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep='!')
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
    if(opt$mode == 'directional'){
      # Preserve historical Estimate sign convention (group1 = out-class, group2 = in-class).
      wil_cox_res<-wilcox.test(rank(res_atc$T) ~ class_bin, conf.int =T)
    } else {
      # Two-sample form so the alternative direction is unambiguous: x = in-class,
      # y = out-of-class, so alternative='greater' tests in-class T > out-of-class T.
      in_T <- rank(res_atc$T)[class_bin == 1]
      out_T <- rank(res_atc$T)[class_bin == 0]
      wil_cox_res<-wilcox.test(in_T, out_T, conf.int =T, alternative = 'greater')
    }

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

write.csv(atc_enrich, paste0(outdir,'/results/',opt$twas,'/twas/drugtargetor/twas_gsea',suffix,'_',opt$panel,'_res_atc_res.csv'), row.names=F)

# Test for enrichment for each level 4 ATC category
res_atc$atc_cat_2<-substr(res_atc$Code, 1, 5)
atc_enrich_2<-NULL
for(cat in unique(res_atc$atc_cat_2)){
  class_bin<-rep(0, nrow(res_atc))
  class_bin[res_atc$atc_cat_2 == cat]<-1

  if(sum(class_bin == 1) > 2){

    # Use wilcoxon test
    if(opt$mode == 'directional'){
      # Preserve historical Estimate sign convention (group1 = out-class, group2 = in-class).
      wil_cox_res<-wilcox.test(rank(res_atc$T) ~ class_bin, conf.int =T)
    } else {
      # Two-sample form so the alternative direction is unambiguous: x = in-class,
      # y = out-of-class, so alternative='greater' tests in-class T > out-of-class T.
      in_T <- rank(res_atc$T)[class_bin == 1]
      out_T <- rank(res_atc$T)[class_bin == 0]
      wil_cox_res<-wilcox.test(in_T, out_T, conf.int =T, alternative = 'greater')
    }

    atc_enrich_2<-rbind(atc_enrich_2, data.frame(ATC=cat,
                                             Estimate=as.numeric(wil_cox_res$estimate),
                                             Class_Median=median(res_atc$Estimate[class_bin == 1]),
                                             Non_Class_Median=median(res_atc$Estimate[class_bin == 0]),
                                             P=wil_cox_res$p.value,
                                             N=sum(class_bin)))
  }
}

atc_labels<-atc[nchar(atc$Code) == 5,]
atc_enrich_2<-merge(atc_enrich_2, atc_labels, by.x='ATC', by.y='Code')
atc_enrich_2<-atc_enrich_2[order(atc_enrich_2$P),]

write.csv(atc_enrich_2, paste0(outdir,'/results/',opt$twas,'/twas/drugtargetor/twas_gsea',suffix,'_',opt$panel,'_res_atc_res_level4.csv'), row.names=F)

