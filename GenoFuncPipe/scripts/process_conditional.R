#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Read in the report files
library(data.table)

# Read in all jointly significant associations
temp = list.files(path=paste0('results/',opt$gwas,'/twas/conditional/'),pattern=glob2rx("*chr*.report"))

report<-do.call(rbind, lapply(temp, function(x) read.table(paste0('results/',opt$gwas,'/twas/conditional/',x), header=T,stringsAsFactors=F)))
report$JOINT.ID<-NA
report$MARGIN.ID<-NA
report$JOINT.N<-NA
report$MARGIN.N<-NA
report$loc<-gsub('.*loc_','',report$FILE)
joint_res<-NULL
margin_res<-NULL

# Insert names of jointly significant genes
for(i in unique(report$CHR)){
  joint_i<-read.table(paste0('results/',opt$gwas,'/twas/conditional/',opt$gwas,'_twas_conditional_chr',i,'.joint_included.dat'), header=T,stringsAsFactors=F)
  margin_i<-read.table(paste0('results/',opt$gwas,'/twas/conditional/',opt$gwas,'_twas_conditional_chr',i,'.joint_dropped.dat'), header=T,stringsAsFactors=F)
  
  temp = list.files(path=paste0('results/',opt$gwas,'/twas/conditional/'), pattern=glob2rx(paste0("*chr",i,".loc*.genes")))

  for(k in gsub('.genes','',gsub('.*loc_','', temp))){
    loc_k<-read.table(paste0('results/',opt$gwas,'/twas/conditional/',opt$gwas,'_twas_conditional_chr',i,'.loc_',k,'.genes'), header=T, stringsAsFactors=F)
    
    loc_k_joint<-loc_k[(loc_k$FILE %in% joint_i$FILE),]
    joint_res<-rbind(joint_res,loc_k_joint)
        
    if(dim(margin_i)[1] > 0){
      loc_k_margin<-loc_k[(loc_k$FILE %in% margin_i$FILE),]
      margin_res<-rbind(margin_res,loc_k_margin)
    } else {
      loc_k_margin<-data.frame(ID=NULL)
    }
    
    g_list<-NULL
    for(g in unique(loc_k_joint$ID)){
      g_list<-c(g_list,paste0(g, " (",paste(loc_k_joint$PANEL[loc_k_joint$ID == g], collapse=', '),")"))
    }
    report[report$CHR == i & report$loc == k,]$JOINT.ID<-paste(g_list,collapse=', ')

    if(dim(loc_k_margin)[1] > 0){
      g_list<-NULL
      for(g in unique(loc_k_margin$ID)){
        g_list<-c(g_list,paste0(g, " (",paste(unique(loc_k_margin$PANEL[loc_k_margin$ID == g]), collapse=', '),")"))
      }
      report[report$CHR == i & report$loc == k,]$MARGIN.ID<-paste(g_list,collapse=', ')
    } else {
      report[report$CHR == i & report$loc == k,]$MARGIN.ID<-'-'
    }
    
    report[report$CHR == i & report$loc == k,]$JOINT.N<-dim(loc_k_joint)[1]
    report[report$CHR == i & report$loc == k,]$MARGIN.N<-dim(loc_k_margin)[1]
  }
}

report$LOCUS<-paste0(report$CHR,':',report$P0,':',report$P1)
report$BP<-paste0(report$P0,'-',report$P1)
report$VAR.EXP<-paste0(report$VAR.EXP*100,'%')

report<-report[,c('CHR','P0','P1','BP','LOCUS',"JOINT.N",'MARGIN.N','BEST.TWAS.P','BEST.SNP.P','VAR.EXP','JOINT.ID','MARGIN.ID')]

report<-report[order(report$CHR, report$P0),]

# Save full conditional results table
write.csv(report[,c("CHR","BP","JOINT.ID","MARGIN.ID","BEST.TWAS.P","BEST.SNP.P","VAR.EXP")],paste0('results/',opt$gwas,'/twas/conditional/',opt$gwas,'_twas_conditional_clean_full.csv'), row.names=F, quote=T)

# Save brief conditional results table
write.csv(report[,c('CHR','BP','JOINT.ID','MARGIN.N','BEST.TWAS.P','BEST.SNP.P','VAR.EXP')],paste0('results/',opt$gwas,'/twas/conditional/',opt$gwas,'_twas_conditional_clean_brief.csv'), row.names=F, quote=T)

# Combine gene results for marginal and joint genes
joint_res$Type<-'Joint'

if(!is.null(margin_res)){
  margin_res$Type<-'Marginal'
  gene_res<-rbind(joint_res, margin_res)
} else {
  gene_res<-joint_res
}

gene_res$Novel<-'No'
gene_res$Novel[(2*pnorm(-abs(gene_res$BEST.GWAS.Z)) < 5e-8 & gene_res$TOP.SNP.COR^2 < 0.1) | 2*pnorm(-abs(gene_res$BEST.GWAS.Z)) > 5e-8]<-'Yes'

gene_res$BP<-paste0(gene_res$P0,'-',gene_res$P1)
gene_res$BEST.GWAS.P<-2*pnorm(-abs(gene_res$BEST.GWAS.Z))

gene_res<-gene_res[order(gene_res$CHR, gene_res$P0),]

gene_res$Colocalised<-F
gene_res$Colocalised[gene_res$COLOC.PP4 >0.8]<-T

gene_res<-gene_res[,c('CHR','BP','P0','P1','ID','PANEL','TWAS.Z','TWAS.P','BEST.GWAS.P','TOP.SNP.COR','Type','Novel','COLOC.PP3','COLOC.PP4','Colocalised')]

# Save table showing whether gene associations are novel
write.csv(gene_res,paste0('results/',opt$gwas,'/twas/',opt$gwas,'_twas_novelty.csv'), row.names=F, quote=T)
