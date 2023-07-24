#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),
  make_option("--chr", action="store", default=NA, type='character',
              help="chromosome number [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(susieR)

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

# Read in the sumstats
ss<-fread(paste0(outdir,'/data/gwas_sumstat/',opt$gwas,'/',opt$gwas,'.cleaned.gz'))
ss<-ss[ss$CHR == opt$chr,]

if(!any(names(ss) == 'BETA')){
  ss$BETA<-log(ss$OR)
}

ss$Z<--qnorm(ss$P/2)
ss$Z<-ss$Z*sign(ss$BETA)

# Read in list of lead GW sig variants
lead<-fread(paste0(outdir,'/results/',opt$gwas,'/clump/',opt$gwas,'.GW.clump.clean.csv'))
lead<-lead[lead$P < 5e-8,]
lead<-lead[lead$CHR == opt$chr,]

if(nrow(lead) == 0){
  file.create(paste0(outdir,'/results/',opt$gwas,'/checks/',opt$gwas,'.chr',opt$chr,'.finemap.done'))
  q()
}

lead$NearestGene<-NULL

# Read in reference SNP data to match alleles
bim<-fread(paste0(outdir,'/resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr',opt$chr,'.bim'))

# Flip sumstat Z to match alleles across gwas and reference
ss_bim_match<-merge(ss, bim, by.x=c('SNP','A1','A2'), by.y=c('V2','V5','V6'))
ss_bim_swap<-merge(ss, bim, by.x=c('SNP','A1','A2'), by.y=c('V2','V6','V5'))

ss<-ss[ss$SNP %in% ss_bim_match$SNP | ss$SNP %in% ss_bim_swap$SNP,] 
ss$Z[ss$SNP %in% ss_bim_swap$SNP]<- -ss$Z[ss$SNP %in% ss_bim_swap$SNP]

dir.create(paste0(outdir,'/results/',opt$gwas,'/finemap'), recursive = T)
for(loc in 1:nrow(lead)){
  # Identify variants within 500kb of lead variants
  ss_subset<-ss[ss$BP > lead$BP[loc] - 5e5 & ss$BP < lead$BP[loc] + 5e5,]
  ss_subset<-ss_subset[order(ss_subset$CHR, ss_subset$BP),]
  write.table(ss_subset$SNP, paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc],'.snp_for_ld.txt'), col.names=F, row.names=F, quote=F)
  
  # Calculate LD matrix for variants surrounding lead variants
  system(paste0('plink --bfile resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr',opt$chr,' --extract results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc],'.snp_for_ld.txt --r square --out results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc]))
  
  ld<-as.matrix(fread(paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc],'.ld')))
  dimnames(ld)<-list(ss_subset$SNP, ss_subset$SNP)
  
  ld[lower.tri(ld)] <- t(ld)[lower.tri(ld)]
  ld<-ld[match(ss_subset$SNP, rownames(ld)),match(ss_subset$SNP, colnames(ld))]

  skip_to_next<-F
  tryCatch(fitted_rss <- susie_rss(ss_subset$Z, ld, L = 10), error = function(e){skip_to_next <<- TRUE})
  
  print(skip_to_next)
  
  if(skip_to_next == F){
    saveRDS(fitted_rss, file = paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc],'.rds'))
  }
  
  # Run finemapping with L=1 (single causal variant), thereby being robust to LD misspecification.
  skip_to_next<-F
  tryCatch(fitted_rss_L1 <- susie_rss(ss_subset$Z, ld, L = 1), error = function(e){skip_to_next <<- TRUE})
  
  print(skip_to_next)

  if(skip_to_next == F){
    saveRDS(fitted_rss_L1, file = paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc],'.L1.rds'))
  }
  
  system(paste0('rm ',outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc],'.ld'))
  system(paste0('rm ',outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc],'.log'))
  system(paste0('rm ',outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.chr',opt$chr,'.',lead$SNP[loc],'.snp_for_ld.txt'))
}

file.create(paste0(outdir,'/results/',opt$gwas,'/checks/',opt$gwas,'.chr',opt$chr,'.finemap.done'))
