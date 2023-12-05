#!/usr/bin/Rscript
library(data.table)

# Create snp_modifyBuild_offline
make_executable <- function(exe) {
  Sys.chmod(exe, mode = (file.info(exe)$mode | "111"))
}

snp_modifyBuild_offline<-function (info_snp, liftOver, chain, from = "hg18", to = "hg19"){
  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop2("Please use proper names for variables in 'info_snp'. Expected %s.",
          "'chr' and 'pos'")
  liftOver <- normalizePath(liftOver)
  make_executable(liftOver)
  BED <- tempfile(fileext = ".BED")
  info_BED <- with(info_snp, data.frame(paste0("chr", chr),
                                        pos0 = pos - 1L, pos, id = seq_len(nrow(info_snp))))
  fwrite(info_BED, BED, col.names = FALSE, sep = " ")
  lifted<-'tmp/tmp.lifted'
  unmapped<-'tmp/tmp.unmapped'
  system(paste(liftOver, BED, chain, lifted, unmapped))
  new_pos <- fread(lifted)
  bad <- grep("^#", readLines(unmapped), value = TRUE, invert = TRUE)
  print(paste0(length(bad)," variants have not been mapped."))
  info_snp$pos <- NA
  info_snp$pos[new_pos$V4] <- new_pos$V3
  info_snp
}

###################
# PLINK format 1000 genomes phase 3 data
###################
# Based on instructions from Hannah Meyer https://cran.r-project.org/web/packages/plinkQC/vignettes/Genomes1000.pdf
# With adaptations from Joni Coleman's reference resource

output_dir<-'resources/data/1kg'
dir.create(output_dir, recursive = T)

# Download the 1000Genomes data provided by PLINK (PLINK 2 format)
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.pgen.zst https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1'))
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.pvar.zst https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1'))
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.psam https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam?dl=1'))

# Decompress the pgen file
system(paste0('plink2 --zst-decompress ',output_dir,'/all_phase3.pgen.zst > ',output_dir,'/all_phase3.pgen'))

# Convert to plink 1 format, subset EUR, and restict to MAF 0.001
system(paste0('plink2 --pfile ',output_dir,'/all_phase3 vzs --max-alleles 2 --make-bed --out ',output_dir,'/1KG.Phase3'))

# Delete the plink2 format data
system(paste0('rm ',output_dir,'/all_phase3.pgen*'))
system(paste0('rm ',output_dir,'/all_phase3.pvar*'))

# Make FID == IID
fam<-read.table(paste0(output_dir,'/1KG.Phase3.fam'), header=F)
fam$V1<-fam$V2
write.table(fam, paste0(output_dir,'/1KG.Phase3.fam'), col.names=F, row.names=F, quote=F)

# Create keep file for each super population
dir.create(paste0(output_dir,'/keep_files'), recursive = T)
pop_data<-read.table(paste0(output_dir, '/all_phase3.psam'), header=F, stringsAsFactors=F)
for(pop in unique(pop_data$V5)){
  write.table(cbind(pop_data$V1[pop_data$V5 == pop],pop_data$V1[pop_data$V5 == pop]), paste0(output_dir,'/keep_files/',pop,'_samples.keep'), col.names=F, row.names=F, quote=F)
}

# Split the plink data by super population, by chromosome (autosome only), restict to variants with MAF > 0.001
for(pop in unique(pop_data$V5)){
  for(chr in 1:22){
    system(paste0('plink --bfile ',output_dir,'/1KG.Phase3 --keep ',output_dir,'/keep_files/',pop,'_samples.keep --maf 0.001 --chr ',chr,' --allow-extra-chr --make-bed --out ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr))
  }
}

# Delete genome wide data
system(paste0('rm ',output_dir,'/1KG.Phase3.bed'))
system(paste0('rm ',output_dir,'/1KG.Phase3.bim'))
system(paste0('rm ',output_dir,'/1KG.Phase3.fam'))
system(paste0('rm ',output_dir,'/1KG.Phase3.log'))
system(paste0('rm ',output_dir,'/all_phase3.psam'))
system(paste0('rm -r ',output_dir,'/keep_files'))

# For autosomes: Remove variants with 3+ alleles or multiple positions
for(pop in unique(pop_data$V5)){
  for(chr in 1:22){
    system(paste0('plink --bfile ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,' --bmerge ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,' --merge-mode 6 --make-bed --out ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,'_noDup'))
  }

  for(chr in 1:22){
    missnp<-fread(paste0(output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,'_noDup.missnp'), header=F)$V1
    log_chr<-readLines(paste0(output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,'_noDup.log'))
    log_chr<-log_chr[grepl("Warning: Multiple positions seen for variant", log_chr)]
    log_chr<-unique(gsub("'\\.","",gsub(".*variant '","",log_chr)))
    missnp<-c(missnp, log_chr)
    write.table(missnp, paste0(output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,'_noDup.missnp'), col.names=F, row.names=F, quote=F)
  }

  for(chr in 1:22){
    system(paste0('plink --bfile ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,' --exclude ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,'_noDup.missnp --make-bed --out ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr))
  }

  system(paste0('rm ',output_dir,'/*',pop,'*_noDup.missnp'))
  system(paste0('rm ',output_dir,'/*',pop,'*_noDup.log'))
  system(paste0('rm ',output_dir,'/*',pop,'*.bed~'))
  system(paste0('rm ',output_dir,'/*',pop,'*.bim~'))
  system(paste0('rm ',output_dir,'/*',pop,'*.fam~'))  
}

####################
# Compute allele frequencies across all individuals
####################

for(pop in unique(pop_data$V5)){
  dir.create(paste0(output_dir,'/frq_files/',pop), recursive = T)
  for(chr in 1:22){
    system(paste0('plink --bfile ',output_dir,'/frq_files/',pop,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,' --freq --out ',output_dir,'/frq_files/',pop,'/1KG.Phase3.',pop,'.MAF_001.chr',chr))
  }
  # Delete the log and nosex files
  system(paste0('rm ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr*.log'))
}

###################
# Prepare summary data for GWAS sumstat cleaning
###################

for(pop in unique(pop_data$V5)){

  dir.create('tmp')

  chrs<-c(1:22)
  for(chr in chrs){
    
    ######
    # Read in the reference data
    ######
    ref<-list()
    ref[['GRCh37']]<-fread(paste0(output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,'.bim'))
    ref[['GRCh37']]$V3<-NULL
    names(ref[['GRCh37']])<-c('chr','snp','pos','a1','a2')
    
    ######
    # Liftover from GRCh37 to GRCh38 and GRCh36
    ######
    
    # Liftover BP to GRCh38
    ref[['GRCh38']]<-snp_modifyBuild_offline(ref[['GRCh37']], liftOver='resources/software/liftover/liftover', chain='resources/data/liftover/hg19ToHg38.over.chain.gz', from = "hg19", to = "hg38")
    # Liftover BP to GRCh36
    ref[['GRCh36']]<-snp_modifyBuild_offline(ref[['GRCh37']], liftOver='resources/software/liftover/liftover', chain='resources/data/liftover/hg19ToHg18.over.chain.gz', from = "hg19", to = "hg18")
    
    # Combine the two builds 
    tmp<-ref[['GRCh37']]
    names(tmp)<-c('CHR','SNP','BP_GRCh37','A1','A2')
    tmp$BP_GRCh38<-ref[['GRCh38']]$pos
    tmp$BP_GRCh36<-ref[['GRCh36']]$pos
    rm(ref)
    
    # Insert IUPAC codes into ref
    tmp$IUPAC[tmp$A1 == 'A' & tmp$A2 =='T' | tmp$A1 == 'T' & tmp$A2 =='A']<-'W'
    tmp$IUPAC[tmp$A1 == 'C' & tmp$A2 =='G' | tmp$A1 == 'G' & tmp$A2 =='C']<-'S'
    tmp$IUPAC[tmp$A1 == 'A' & tmp$A2 =='G' | tmp$A1 == 'G' & tmp$A2 =='A']<-'R'
    tmp$IUPAC[tmp$A1 == 'C' & tmp$A2 =='T' | tmp$A1 == 'T' & tmp$A2 =='C']<-'Y'
    tmp$IUPAC[tmp$A1 == 'G' & tmp$A2 =='T' | tmp$A1 == 'T' & tmp$A2 =='G']<-'K'
    tmp$IUPAC[tmp$A1 == 'A' & tmp$A2 =='C' | tmp$A1 == 'C' & tmp$A2 =='A']<-'M'
    
    # Read in reference frequency data
    freq<-fread(paste0(output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,'.frq'))
    
    # The freq files have come from the reference files, so we can assume they are on the same strand
    freq_match<-merge(tmp, freq[,c('SNP','A1','A2','MAF'), with=F], by=c('SNP','A1','A2'))
    freq_swap<-merge(tmp, freq[,c('SNP','A1','A2','MAF'), with=F], by.x=c('SNP','A1','A2'), by.y=c('SNP','A2','A1'))
      freq_swap$MAF<-1-freq_swap$MAF
      tmp_freq<-rbind(freq_match, freq_swap)
      tmp_freq<-tmp_freq[match(tmp$SNP, tmp_freq$SNP),]
      
      tmp[['REF.FRQ']]<-tmp_freq$MAF

    tmp<-tmp[,c("CHR","SNP","BP_GRCh36","BP_GRCh37","BP_GRCh38","A1","A2","IUPAC","REF.FRQ"), with=F]
    saveRDS(tmp, file = paste0(output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr,'.rds'))
  }

  system('rm -r tmp')
}
