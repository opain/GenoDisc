#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
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

# Read in MAGMA gene set results
res_gs<-fread(cmd=paste0("grep -v '^#' ",outdir,"/results/",opt$gwas,'/magma/magma_drug_targetor.gsa.out'))

# Remove gene sets with <5 genes present
#res_gs<-res_gs[res_gs$NGENES >= 5,]

res_gs$ATC<-gsub('ATC:','',gsub('\\|.*','',res_gs$FULL_NAME))
res_gs$NAME<-tolower(gsub('\\|.*','',gsub('.*NAME:','',res_gs$FULL_NAME)))
res_gs$NAME<-gsub('_',' ',res_gs$NAME)

res_gs<-res_gs[,c('NAME','NGENES','BETA','SE','P','ATC'), with=F]

write.csv(res_gs, paste0(outdir,'/results/',opt$gwas,'/magma/magma_drug_targetor.clean.csv'), row.names=F, quote=T)

# Canonical membership: each drug's own ATC codes (comma-joined in res_gs$ATC), exploded to every
# class, unique per drug; the rank-sum runs over UNIQUE drugs (no drug*code pseudo-replication).
# Test for enrichment for each ATC L3 class (threshold >= 5).
drug_cls3 <- lapply(strsplit(res_gs$ATC, ',', fixed = TRUE), function(x) unique(substr(x, 1, 4)))
cls3      <- sort(unique(unlist(drug_cls3))); cls3 <- cls3[nchar(cls3) == 4]
atc_enrich<-NULL
for(cat in cls3){
  class_bin <- as.numeric(vapply(drug_cls3, function(x) cat %in% x, logical(1)))

  if(sum(class_bin == 1) >= 5){

    wil_cox_res<-wilcox.test(rank(res_gs$P) ~ class_bin, conf.int =T, alternative='greater')

    atc_enrich<-rbind(atc_enrich, data.frame(ATC=cat,
                                             Estimate=as.numeric(wil_cox_res$estimate),
                                             Class_Median=median(res_gs$P[class_bin == 1]),
                                             Non_Class_Median=median(res_gs$P[class_bin == 0]),
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

write.csv(atc_enrich, paste0(outdir,'/results/',opt$gwas,'/magma/magma_drug_targetor_atc_res.csv'), row.names=F)

# Test for enrichment for each level 4 ATC class, over UNIQUE drugs (threshold >= 5).
drug_cls4 <- lapply(strsplit(res_gs$ATC, ',', fixed = TRUE), function(x) unique(substr(x, 1, 5)))
cls4      <- sort(unique(unlist(drug_cls4))); cls4 <- cls4[nchar(cls4) == 5]
atc_enrich_2<-NULL
for(cat in cls4){
  class_bin <- as.numeric(vapply(drug_cls4, function(x) cat %in% x, logical(1)))

  if(sum(class_bin == 1) >= 5){

    wil_cox_res<-wilcox.test(rank(res_gs$P) ~ class_bin, conf.int =T, alternative='greater')

    atc_enrich_2<-rbind(atc_enrich_2, data.frame(ATC=cat,
                                             Estimate=as.numeric(wil_cox_res$estimate),
                                             Class_Median=median(res_gs$P[class_bin == 1]),
                                             Non_Class_Median=median(res_gs$P[class_bin == 0]),
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


