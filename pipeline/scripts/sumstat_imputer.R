#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics file [required]"),
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
              help="Path to per chromosome PLINK files [required]"),
  make_option("--fizi", action="store", default=NA, type='character',
              help="Path to FIZI binary [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
              help="Number of cores for parallel computing [optional]"),
  make_option("--min_prop", action="store", default=0.5, type='numeric',
              help="Minimum proportion of SNPs in region to impute [optional]"),
  make_option("--output", action="store", default='./Output', type='character',
              help="Path for output files [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.imputed.log',sep=''), append = F)
cat(
  '#################################################################
# sumstat_imputer.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Munge sumstats manually
#####

# Requires CHR, SNP, BP, A1, A2, Z, SE, and N.
# Read in sumstats
gwas<-fread(opt$sumstats)

# Insert BETA if not present
if(!('BETA' %in% names(gwas))){
  gwas$BETA<-log(gwas$OR)
}

# Calculate Z an sign based on BETA
gwas$Z<-abs(qnorm(gwas$P/2))
gwas$Z<-sign(gwas$BETA)*gwas$Z

# Subset relevent columns
gwas<-gwas[,c('CHR','SNP','BP','A1','A2','Z','SE','N')]

# Round numbers as is done in fizi munge function
gwas$Z<-round(gwas$Z,3)
gwas$SE<-round(gwas$SE,3)

# Remove any rows with missing data
gwas<-gwas[complete.cases(gwas),]

# Write out munged sumstats
fwrite(gwas, paste0(opt$output,'.munged.tmp.gz'), sep='\t')

sink(file = paste(opt$output,'.imputed.log',sep=''), append = T)
cat('Before imputation, sumstats contain',nrow(gwas),'variants.\n')
sink()

rm(gwas)

####
# Run imputation using ImpG algorithm
####

# This takes quite a while, so an option to run in parallel should be implemented.
if(opt$n_cores == 1){
  for(i in 1:22){
    system(paste0(opt$fizi,' impute ',opt$output,'.munged.tmp.gz ',opt$ref_plink_chr,i,' --chr ',i,' --min-prop ',opt$min_prop,' --out ', opt$output,'.imped.chr',i))
  }
} else {
  # Make a data.frame listing chromosome and phi combinations
  jobs<-NULL
  for(i in 1:22){
    jobs<-rbind(jobs,data.frame(CHR=i))
  }

  write.table(jobs, paste0(opt$output,'.job_list'), col.names=F, row.names=F, quote=F)

  # Write batch job
  writeLines(paste0("#!/bin/sh

#SBATCH -p neurohack_cpu,cpu
#SBATCH --mem 15G
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH -t 3:00:00
#SBATCH -J fizi

export MKL_NUM_THREADS=$SLURM_CPUS_ON_NODE
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_ON_NODE
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

echo $SLURM_CPUS_ON_NODE
chr=$(awk -v var=$SLURM_ARRAY_TASK_ID 'NR == var {print $1}' ", opt$output,".job_list)

echo ${chr}

",opt$fizi,' impute ',opt$output,'.munged.tmp.gz ',opt$ref_plink_chr,'${chr} --chr ${chr}  --min-prop ',opt$min_prop,' --out ', opt$output,".imped.chr${chr}

"), paste0(opt$output,'.batch.sh'))

  # Run batch job
  jobID<-system(paste0("sbatch --array ",1,"-",nrow(jobs),"%",opt$n_cores," ", opt$output,'.batch.sh'),intern=T)
  jobID<-gsub('.* ','', jobID)

  # Check whether finished
  Sys.sleep(30)
  while(i){
    system(paste0('sacct -j ',jobID,' > ',opt$output,'.sacct_log.txt'))
    sacct_log<-fread(paste0(opt$output,'.sacct_log.txt'), fill=T)
    sacct_log<-sacct_log[sacct_log$JobName == 'fizi',]

    print(sacct_log)

    if(sum(sacct_log$State == 'FAILED') > 0){
      sink(file = paste(opt$output,'.imputed.log',sep=''), append = T)
      cat('Job failed.\n')
      sink()
      q()
    }


    if(sum(sacct_log$State != 'COMPLETED') == 0){
      break
    } else {
      if(sum(sacct_log$State != 'RUNNING' & sacct_log$State != 'COMPLETED' & sacct_log$State != 'PENDING' & sacct_log$State != 'NODE_FAIL') > 0){
        stop()
      }
      Sys.sleep(60)
    }
  }
}

####
# Format imputed sumstats
####

imped<-NULL
for(i in 1:22){
  log<-readLines(paste0(opt$output,'.imped.chr',i,'.log'))

  if(!any(grepl('ERROR', log))){
    imped<-rbind(imped, fread(paste0(opt$output,'.imped.chr',i,'.sumstat')))
  } else {
    sink(file = paste(opt$output,'.imputed.log',sep=''), append = T)
    cat('At least one chromosome failed to complete. Check log files.\n')
    sink()

    stop()
  }
}

# Remove imputed variant with R2 < 0.8
imped<-imped[!(imped$R2.BLUP < 0.8),]

sink(file = paste(opt$output,'.imputed.log',sep=''), append = T)
cat('After applying R2 >= 0.8 threshold, sumstats contain',nrow(imped),'variants.\n')
cat(sum(imped$TYPE == 'imputed'),'variants are imputed.\n')
cat(sum(imped$TYPE == 'gwas'),'variants were present in the GWAS\n')
sink()

names(imped)<-c('CHR','SNP','BP','A1','A2','TYPE','Z','R2.BLUP','N','P')
imped$N<-round(imped$N)
imped<-imped[complete.cases(imped),]

fwrite(imped, paste0(opt$output,'.imputed.gz'), sep=' ')

####
# Clean up directory
####

system(paste0('rm ',opt$output,'.munged.tmp.gz'))
system(paste0('rm ',opt$output,'.job_list'))
system(paste0('rm ',opt$output,'.batch.sh'))
system(paste0('rm ',opt$output,'.sacct_log.txt'))
system(paste0('rm ',opt$output,'.imped.chr*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.imputed.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
