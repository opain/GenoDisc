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

# Subset sets with FDR < 0.05 or top 10 sets
res_gs$P.FDR<-p.adjust(res_gs$P, method='fdr')
if(sum(res_gs$P.FDR < 0.05) < 5){
  res_gs_subset<-res_gs[order(res_gs$P),][1:10,]
} else {
  res_gs_subset<-res_gs[res_gs$P.FDR < 0.05,]
}

# Read in drug targetor data
drugtargetor<-fread('resources/data/drug_targetor/wholedatabase_for_targetor')
drugtargetor$activity_score<-NA
drugtargetor$activity_score[drugtargetor$activity_type %in% c('DECREASED_EXPRESSION','NEGATIVE_RESPONSE','OPPOSITE_RESPONSE')]<- -1
drugtargetor$activity_score[drugtargetor$activity_type %in% c('INCREASED_EXPRESSION','POSITIVE_RESPONSE')]<- 1
drugtargetor<-drugtargetor[!is.na(drugtargetor$activity_score),]

# Read in the TWAS results
twas<-fread(paste0(outdir,'/results/',opt$gwas,'/twas/',opt$gwas,'_twas_GW_clean.txt.gz'))
twas<-twas[!grepl('SPLICE', twas$PANEL),]
twas$TWAS.P.FDR<-p.adjust(twas$TWAS.P, method='fdr')

# For each drug, list TWAS associations for genes interacting with the drug
library(ggplot2)
library(cowplot)

twas_drug_cor<-NULL
plots<-list()

for(i in 1:nrow(res_gs_subset)){
  drugtargetor_i<-drugtargetor[grepl(res_gs_subset$FULL_NAME[i], drugtargetor$atc),]
  twas_i<-twas[twas$external_gene_name %in% drugtargetor_i$gene,]
  twas_i<-twas_i[!is.na(twas_i$TWAS.Z),]
  
  # Calculate directional consistency
  if(any(duplicated(twas_i$external_gene_name))){
    twas_i_mean<-aggregate(twas_i$TWAS.Z, list(twas_i$external_gene_name), FUN=mean)
    names(twas_i_mean)<-c('ID','Z')
  } else {
    twas_i_mean<-data.frame(twas_i$external_gene_name, twas_i$TWAS.Z)
    names(twas_i_mean)<-c('ID','Z')
  }
  
  twas_i_mean$Dir<-ifelse(twas_i_mean$Z > 0, 1, -1)

  if(any(duplicated(drugtargetor_i$gene))){
    drugtargetor_i_mean<-aggregate(drugtargetor_i$activity_score, list(drugtargetor_i$gene), FUN=mean)
    names(drugtargetor_i_mean)<-c('ID','Dir')
  } else {
    drugtargetor_i_mean<-data.frame(drugtargetor_i$gene, drugtargetor_i$activity_score)
    names(drugtargetor_i_mean)<-c('ID','Dir')
  }
  
  both_i<-merge(twas_i_mean, drugtargetor_i_mean, by='ID')
  
  twas_drug_cor<-rbind(twas_drug_cor, data.frame(Drug=res_gs_subset$FULL_NAME[i],
                                                 N=nrow(both_i),
                                                 N_con=sum(both_i$Dir.x == both_i$Dir.y),
                                                 N_dis=sum(both_i$Dir.x != both_i$Dir.y)))
  
  twas_i_clean<-twas_i
  twas_i_clean$Sig<-twas_i_clean$TWAS.P.FDR < 0.05
  twas_i_clean$Coloc<-twas_i_clean$TWAS.P.FDR < 0.1 & (twas_i_clean$COLOC.PP4-twas_i_clean$COLOC.PP3)/twas_i_clean$COLOC.PP4 > 0.8
  
  twas_i_clean<-twas_i_clean[,c('external_gene_name','TWAS.Z','Sig','Coloc','PANEL'),with=F]
  names(twas_i_clean)<-c('ID','Z','Sig','Coloc','Panel')

  drugtargetor_i_mean_clean<-drugtargetor_i_mean
  names(drugtargetor_i_mean_clean)<-c('ID','DrugTargetor')

  if(length(unique(twas_i$external_gene_name[twas_i$TWAS.P.FDR < 0.05])) < 20){
    twas_i_lead<-unique(twas_i$external_gene_name[order(twas_i$TWAS.P)])[1:20]
  } else {
    twas_i_lead<-unique(twas_i$external_gene_name[twas_i$TWAS.P.FDR < 0.05])
  }

  both_merge_i<-merge(twas_i_clean,drugtargetor_i_mean_clean, by='ID')
  # Flip TWAS.Z according to DrugTargetor direction to highlight drugs that reverse risk
  both_merge_i$Z<- both_merge_i$Z * both_merge_i$DrugTargetor
  # Negative Z means drug reverse risk
  
  both_merge_i<-both_merge_i[both_merge_i$ID %in% twas_i_lead,]

  plots[[i]]<-ggplot(data = both_merge_i, aes(x = Panel, y = ID)) +
    theme_bw()	+
    geom_tile(data=both_merge_i, aes(fill = Z), width=0.95, height=0.95) +
    geom_tile(data=both_merge_i[both_merge_i$Sig == T,], aes(x = Panel, y = ID), colour='black', fill=NA, size=0.3, width=0.95, height=0.95) +
    geom_tile(data=both_merge_i[both_merge_i$Coloc ==T & both_merge_i$Sig == T,], aes(x = Panel, y = ID), colour='green2', fill=NA, size=0.3, width=0.95, height=0.95) +
    scale_fill_gradientn(colours=c("dodgerblue2","white","red"), na.value = NA,name = "Z-score", limits = c(-max(abs(both_merge_i$Z)),max(abs(both_merge_i$Z)))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
    labs(x='', y='', title=tolower(gsub('\\|.*','',gsub('.*NAME:','',res_gs_subset$FULL_NAME[i]))))

}

bitmap(paste0(outdir,'/results/',opt$gwas,'/magma/magma_drug_targetor_twas_comp.png'), units='px',res=300, height=5000, width=4000)
plot_grid(plotlist=plots, ncol=2)
dev.off()


