#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics file [required]"),
  make_option("--ref_chr", action="store", default=NA, type='character',
              help="Path to per chromosome reference .rds files [required]"),
  make_option("--info", action="store", default=0.9, type='numeric',
              help="INFO threshold [optional]"),
  make_option("--maf", action="store", default=0.01, type='numeric',
              help="MAF threshold [optional]"),
  make_option("--maf_diff", action="store", default=0.2, type='numeric',
              help="Difference between reference and reported MAF threshold [optional]"),
  make_option("--gz", action="store", default=T, type='logical',
              help="Set to T to gzip summary statistics [optional]"),
  make_option("--output", action="store", default='./Output', type='character',
              help="Path for output files [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
  '#################################################################
# sumstat_cleaner.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Read in sumstats
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS sumstats.\n')
sink()

GWAS<-fread(opt$sumstats)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

#####
# Insert IUPAC codes into target
#####

# Create IUPAC code function
snp_iupac<-function(x=NA, y=NA){
  if(length(x) != length(y)){
    print('x and y are different lengths')
  } else {
    iupac<-rep(NA, length(x))
    iupac[x == 'A' & y =='T' | x == 'T' & y =='A']<-'W'
    iupac[x == 'C' & y =='G' | x == 'G' & y =='C']<-'S'
    iupac[x == 'A' & y =='G' | x == 'G' & y =='A']<-'R'
    iupac[x == 'C' & y =='T' | x == 'T' & y =='C']<-'Y'
    iupac[x == 'G' & y =='T' | x == 'T' & y =='G']<-'K'
    iupac[x == 'A' & y =='C' | x == 'C' & y =='A']<-'M'
    return(iupac)
  }
}

# Insert IUPAC codes into target
GWAS$IUPAC<-snp_iupac(GWAS$A1, GWAS$A2)

# Retain only non-ambiguous SNPs
GWAS<-GWAS[(GWAS$IUPAC %in% c('R', 'Y', 'K', 'M')),]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of variants that are not SNPs or are ambiguous,',dim(GWAS)[1],'variants remain.\n')
sink()

#####
# Harmonise per chromosome with reference
#####

# Check whether RSIDs are available for majority of SNPs in GWAS
chr_bp_avail<-sum(c('CHR','ORIGBP') %in% names(GWAS)) == 2 
rsid_avail<-(sum(grepl('rs', GWAS$SNP)) > 0.9*length(GWAS$SNP))

# Create function to change allele to complement
snp_allele_comp<-function(x=NA){
  x_new<-x
  x_new[x == 'A']<-'T'
  x_new[x == 'T']<-'A'
  x_new[x == 'G']<-'C'
  x_new[x == 'C']<-'G'
  x_new[!(x %in% c('A','T','G','C'))]<-NA
  return(x_new)
}

target_build<-NA

if(chr_bp_avail){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Using CHR, BP, A1 and A2 to merge with the reference.\n')
  sink()
  
  ###
  # Determine build
  ###
  # Read in chromosome 22 of reference data
  i<-22
  
  tmp<-readRDS(file = paste0(opt$ref_chr,i,'.rds'))
  
  # Check target-ref condordance of BP across builds
  tmp$CHR<-as.character(tmp$CHR)
  GWAS$CHR<-as.character(GWAS$CHR)
  matched<-list()
  matched[['GRCh36']]<-merge(GWAS, tmp, by.x=c('CHR','ORIGBP','IUPAC'), by.y=c('CHR','BP_GRCh36','IUPAC'))
  matched[['GRCh37']]<-merge(GWAS, tmp, by.x=c('CHR','ORIGBP','IUPAC'), by.y=c('CHR','BP_GRCh37','IUPAC'))
  matched[['GRCh38']]<-merge(GWAS, tmp, by.x=c('CHR','ORIGBP','IUPAC'), by.y=c('CHR','BP_GRCh38','IUPAC'))
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('GRCh36 match: ',round(nrow(matched[['GRCh36']])/sum(GWAS$CHR == i)*100, 2),'%\n',sep='')
  cat('GRCh37 match: ',round(nrow(matched[['GRCh37']])/sum(GWAS$CHR == i)*100, 2),'%\n',sep='')
  cat('GRCh38 match: ',round(nrow(matched[['GRCh38']])/sum(GWAS$CHR == i)*100, 2),'%\n',sep='')
  sink()
  
  if((nrow(matched[['GRCh36']])/sum(GWAS$CHR == i)) > 0.7){
    target_build<-'GRCh36'
  }

  if((nrow(matched[['GRCh37']])/sum(GWAS$CHR == i)) > 0.7){
    target_build<-'GRCh37'
  }
  
  if((nrow(matched[['GRCh38']])/sum(GWAS$CHR == i)) > 0.7){
    target_build<-'GRCh38'
  }
  
  rm(matched,tmp)
  
  if(!is.na(target_build)){
    # Build detected, so can continue reference harmonisation
    # Run per chromosome to reduce memory requirements
    chrs<-1:22
    GWAS_matched<-NULL
    for(i in chrs){
      print(i)
      
      # Read reference data
      tmp<-readRDS(file = paste0(opt$ref_chr,i,'.rds'))
      
      # Subset relevent data
      tmp<-tmp[nchar(tmp$A1) == 1 & nchar(tmp$A2) == 1,]
      names(tmp)[names(tmp) == 'CHR']<-'REF.CHR'
      names(tmp)[names(tmp) == 'SNP']<-'REF.SNP'
      names(tmp)[names(tmp) == 'BP_GRCh36']<-'REF.BP_GRCh36'
      names(tmp)[names(tmp) == 'BP_GRCh37']<-'REF.BP_GRCh37'
      names(tmp)[names(tmp) == 'BP_GRCh38']<-'REF.BP_GRCh38'
      names(tmp)[names(tmp) == 'REF.FRQ']<-'REF.FREQ'
      
      # Merge target and reference by BP
      GWAS_chr<-GWAS[GWAS$CHR == i,]
      ref_target<-merge(GWAS_chr, tmp, by.x='ORIGBP', by.y=paste0('REF.BP_',target_build))
      
      # Identify SNPs that are opposite strands
      flipped<-ref_target[(ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'Y') | 
                            (ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'R') | 
                            (ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'M') |
                            (ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'K'),]
      
      # Change target alleles to compliment for flipped variants
      flipped$A1.x<-snp_allele_comp(flipped$A1.x)
      flipped$A2.x<-snp_allele_comp(flipped$A2.x)
      
      # Update IUPAC codes
      flipped$IUPAC.x<-snp_iupac(flipped$A1.x, flipped$A2.x)
      
      # Identify SNPs that have matched alleles
      matched<-ref_target[ref_target$IUPAC.x == ref_target$IUPAC.y,]
      matched<-rbind(matched, flipped)
      
      # Flip REF.FREQ if alleles are swapped
      matched$REF.FREQ[matched$A1.x != matched$A1.y]<-1-matched$REF.FREQ[matched$A1.x != matched$A1.y]
      
      # Retain reference CHR, BP, SNP information
      matched$A1<-matched$A1.x
      matched$A1.y<-NULL
      matched$A1.x<-NULL
      matched$A2<-matched$A2.x
      matched$A2.y<-NULL
      matched$A2.x<-NULL
      matched$IUPAC.y<-NULL
      matched$IUPAC.x<-NULL
      matched$SNP<-matched$REF.SNP
      matched$REF.SNP<-NULL
      matched$CHR<-matched$REF.CHR
      matched$REF.CHR<-NULL
      matched$BP<-matched$ORIGBP
      matched$ORIGBP<-NULL
      matched<-matched[,!grepl('REF.BP_GRCh',names(matched)), with=F]

      GWAS_matched_chr<-rbind(matched)
      GWAS_matched<-rbind(GWAS_matched, GWAS_matched_chr)
    }
  }
}

if(is.na(target_build) & rsid_avail){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Using SNP, A1 and A2 to merge with the reference.\n')
  sink()
  
  # Run per chromosome to reduce memory requirements
  chrs<-c(1:22)
  
  GWAS_matched<-NULL
  for(i in chrs){
    print(i)
    
    # Read reference data
    tmp<-readRDS(file = paste0(opt$ref_chr,i,'.rds'))
    
    # Subset relevent data
    tmp<-tmp[nchar(tmp$A1) == 1 & nchar(tmp$A2) == 1,]
    names(tmp)[names(tmp) == 'CHR']<-'REF.CHR'
    names(tmp)[names(tmp) == 'BP_GRCh36']<-'REF.BP_GRCh36'
    names(tmp)[names(tmp) == 'BP_GRCh37']<-'REF.BP_GRCh37'
    names(tmp)[names(tmp) == 'BP_GRCh38']<-'REF.BP_GRCh38'
    names(tmp)[names(tmp) == 'REF.FRQ']<-'REF.FREQ'
    
    # Merge target and reference by SNP ID
    ref_target<-merge(GWAS, tmp, by='SNP')
    
    # Identify SNPs that are opposite strands
    flipped<-ref_target[(ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'Y') | 
                          (ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'R') | 
                          (ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'M') |
                          (ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'K'),]

    # Change target alleles to compliment for flipped variants
    flipped$A1.x<-snp_allele_comp(flipped$A1.x)
    flipped$A2.x<-snp_allele_comp(flipped$A2.x)
    
    # Update IUPAC codes
    flipped$IUPAC.x<-snp_iupac(flipped$A1.x, flipped$A2.x)
    
    # Identify SNPs that have matched alleles
    matched<-ref_target[ref_target$IUPAC.x == ref_target$IUPAC.y,]
    matched<-rbind(matched, flipped)
    
    # Flip REF.FREQ if alleles are swapped
    matched$REF.FREQ[matched$A1.x != matched$A1.y]<-1-matched$REF.FREQ[matched$A1.x != matched$A1.y]
    
    # Retain reference CHR, BP, SNP information
    matched$A1<-matched$A1.x
    matched$A1.y<-NULL
    matched$A1.x<-NULL
    matched$A2<-matched$A2.x
    matched$A2.y<-NULL
    matched$A2.x<-NULL
    matched$IUPAC.y<-NULL
    matched$IUPAC.x<-NULL
    matched$BP<-matched$REF.BP_GRCh37
    matched$REF.BP_GRCh37<-NULL
    matched$REF.BP_GRCh38<-NULL
    matched$REF.BP_GRCh36<-NULL
    matched$ORIGBP<-NULL
    matched$CHR<-matched$REF.CHR
    matched$REF.CHR<-NULL
    
    GWAS_matched_chr<-rbind(matched)
    GWAS_matched<-rbind(GWAS_matched, GWAS_matched_chr)
  }
}

GWAS<-GWAS_matched

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After matching variants to the reference,',dim(GWAS)[1],'variants remain.\n')
sink()

#####
# Remove SNPs with INFO < opt$info
#####

if(sum(names(GWAS) == 'INFO') == 1){
  GWAS<-GWAS[GWAS$INFO >= opt$info,]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with INFO < ',opt$info,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('INFO column is not present\n', sep='')
  sink()
}

#####
# Remove SNPs with reported MAF < opt$maf
#####

if(sum(names(GWAS) == 'FREQ') == 1){
  GWAS<-GWAS[GWAS$FREQ >= opt$maf & GWAS$FREQ <= (1-opt$maf),]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with reported MAF < ',opt$maf,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Reported MAF column is not present\n', sep='')
  sink()
}

#####
# Remove SNPs with reference MAF < opt$maf
#####

if(sum(names(GWAS) == 'FREQ') == 1){
  GWAS<-GWAS[GWAS$REF.FREQ >= opt$maf & GWAS$REF.FREQ <= (1-opt$maf),]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with reported MAF < ',opt$maf,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Reported MAF column is not present\n', sep='')
  sink()
}

#####
# Remove SNPs with discordant MAF
#####

if(sum(names(GWAS) == 'FREQ') == 1){
  GWAS$diff<-abs(GWAS$FREQ-GWAS$REF.FREQ)
  
  bitmap(paste0(opt$output,'.MAF_plot.png'), unit='px', res=300, width=1200, height=1200)
    plot(GWAS$REF.FREQ[GWAS$diff > opt$maf_diff],GWAS$FREQ[GWAS$diff > opt$maf_diff], xlim=c(0,1), ylim=c(0,1), xlab='Reference Allele Frequency', ylab='Sumstat Allele Frequency')
    abline(coef = c(0,1))
  dev.off()
  
  GWAS<-GWAS[GWAS$diff < opt$maf_diff,]
  GWAS$diff<-NULL
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with absolute MAF difference of < ',opt$maf_diff,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Reported MAF column is not present, so discordance with reference cannot be determined.\n', sep='')
  sink()
}

#####
# Remove SNPs with out-of-bounds p-values
#####

GWAS<-GWAS[GWAS$P <= 1 & GWAS$P > 0,]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of SNPs with out-of-bound P values, ',dim(GWAS)[1],' variants remain.\n', sep='')
sink()

#####
# Remove SNPs with duplicated rs numbers
#####

dups<-GWAS$SNP[duplicated(GWAS$SNP)]
GWAS<-GWAS[!(GWAS$SNP %in% dups),]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of SNPs with duplicate IDs, ',dim(GWAS)[1],' variants remain.\n', sep='')
sink()

#####
# Remove SNPs with N < 3SD from median N
#####

if(length(unique(GWAS$N)) > 1){
  N_sd<-sd(GWAS$N)
  GWAS<-GWAS[GWAS$N < median(GWAS$N)+(3*N_sd) & GWAS$N > median(GWAS$N)-(3*N_sd),]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with N > ',median(GWAS$N)+(3*N_sd),' or < ',median(GWAS$N)-(3*N_sd),', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('N column is not present or invariant.\n', sep='')
  sink()
}

#####
# Check for genomic control
#####

if(sum(names(GWAS) == 'SE') == 1){
  if(sum(names(GWAS) == 'OR') == 1){
    GWAS$BETA<-log(GWAS$OR)
    GWAS$Z<-GWAS$BETA/GWAS$SE
    GWAS$P_check<-2*pnorm(-abs(GWAS$Z))
    GWAS$Z<-NULL
    GWAS$BETA<-NULL
    
    if(abs(mean(GWAS$P[!is.na(GWAS$P_check)]) - mean(GWAS$P_check[!is.na(GWAS$P_check)])) > 0.01){
       GWAS$P<-GWAS$P_check
       GWAS$P_check<-NULL
      
       sink(file = paste(opt$output,'.log',sep=''), append = T)
       cat('Genomic control detected. P-value recomputed using OR and SE.\n', sep='')
       sink()
    } else {
       sink(file = paste(opt$output,'.log',sep=''), append = T)
       cat('Genomic control was not detected.\n', sep='')
       sink()
       GWAS$P_check<-NULL
    }
  }
  
  if(sum(names(GWAS) == 'BETA') == 1){
    GWAS$Z<-GWAS$BETA/GWAS$SE
    GWAS$P_check<-2*pnorm(-abs(GWAS$Z))
    GWAS$Z<-NULL

    if(abs(mean(GWAS$P[!is.na(GWAS$P_check)]) - mean(GWAS$P_check[!is.na(GWAS$P_check)])) > 0.01){
      GWAS$P<-GWAS$P_check
      GWAS$P_check<-NULL
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Genomic control detected. P-value recomputed using BETA and SE.\n', sep='')
      sink()
    } else {
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Genomic control was not detected.\n', sep='')
      sink()
      GWAS$P_check<-NULL
    }
    
  }
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('SE column is not present, genomic control cannot be detected.\n', sep='')
  sink()
}

#####
# Insert SE column
#####

if(sum(names(GWAS) == 'SE') == 0){
  if(sum(names(GWAS) == 'BETA') == 1){
    GWAS$Z<-abs(qnorm(GWAS$P/2))
    GWAS$SE<-abs(GWAS$BETA/GWAS$Z)
  } else {
    GWAS$Z<-abs(qnorm(GWAS$P/2))
    GWAS$SE<-abs(log(GWAS$OR)/GWAS$Z)
  }
  
  GWAS$Z<-NULL
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('SE column inserted based on BETA/OR and P.\n', sep='')
  sink()
}

#####
# Remove SNPs with SE == 0
#####

GWAS<-GWAS[GWAS$SE != 0,]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of SNPs with SE == 0, ',dim(GWAS)[1],' variants remain.\n', sep='')
sink()

#####
# Write out results
#####

if(file.exists(paste0(opt$output,'.gz'))){
  system(paste0(paste0('rm ',opt$output,'.gz')))
}
if(file.exists(paste0(opt$output))){
  system(paste0(paste0('rm ',opt$output)))
}

if(opt$gz == T){
  fwrite(GWAS, paste0(opt$output,'.gz'), sep='\t')
} else {
  fwrite(GWAS, opt$output, sep='\t')
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
