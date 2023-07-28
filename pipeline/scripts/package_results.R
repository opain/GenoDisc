#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
source('scripts/functions/package_results_functions.R')

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

##################
# GWAS QC
##################

gwas_qc<-list()

gwas_qc$cleaner_dat<-process_cleaner_log(config_file=opt$config_file, gwas=opt$gwas)
gwas_qc$focus_dat<-process_focus_log(config_file=opt$config_file, gwas=opt$gwas)
gwas_qc$ldsc_dat<-process_ldsc_log(config_file=opt$config_file, gwas=opt$gwas)

##################
# SNP associations
##################

snp_assoc<-list()

if(config[grepl('clump:',config)] == "clump: T"){
  snp_assoc$clump<-fread(paste0(outdir,'/results/',opt$gwas,'/clump/',opt$gwas,'.GW.clump.clean.csv'))
}

if(config[grepl('cojo:',config)] == "cojo: T"){
  snp_assoc$cojo<-fread(paste0(outdir,'/results/',opt$gwas,'/cojo/',opt$gwas,'.GW.cojo.clean.csv'))
}

if(config[grepl('finemap:',config)] == "finemap: T"){
  snp_assoc$susie<-list()
  snp_assoc$susie$L10<-process_susie(outdir=outdir, gwas=opt$gwas, L=10)
  snp_assoc$susie$L1<-process_susie(outdir=outdir, gwas=opt$gwas, L=1)
}

#################
# Molecular associations
#################

mol_assoc<-list()

######
# Expression
######

mol_assoc$exp<-list()

###
# FUSION
###

mol_assoc$exp$fusion<-read_fusion_exp(config_file=opt$config_file, gwas=opt$gwas)

###
# SMR
###

mol_assoc$exp$smr<-read_smr_exp(config_file=opt$config_file, gwas=opt$gwas)

######
# Protein
######

mol_assoc$protein<-list()

###
# FUSION
###

mol_assoc$protein$fusion<-read_fusion_protein(config_file=opt$config_file, gwas=opt$gwas)

###
# SMR
###

mol_assoc$exp$smr<-read_smr_protein(config_file=opt$config_file, gwas=opt$gwas)

######
# MAGMA
######

mol_assoc$magma<-read_magma(config_file=opt$config_file, gwas=opt$gwas)

######
# Nearest
######

mol_assoc$nearest<-list()

if(config[grepl('clump:',config)] == "clump: T"){
  mol_assoc$nearest$clump<-identify_nearest(snp_assoc$clump$NearestGene)
}

if(config[grepl('cojo:',config)] == "cojo: T"){
  mol_assoc$nearest$cojo<-identify_nearest(snp_assoc$cojo$NearestGene)

}

######
# Finemapping
######

mol_assoc$finemap<-list()

if(config[grepl('finemap:',config)] == "finemap: T"){
  mol_assoc$finemap$L1<-unlist(strsplit(snp_assoc$susie$L1$Gene, ', '))
  mol_assoc$finemap$L1<-mol_assoc$finemap$L1[mol_assoc$finemap$L1 != 'None']
  mol_assoc$finemap$L10<-unlist(strsplit(snp_assoc$susie$L10$Gene, ', '))
  mol_assoc$finemap$L10<-mol_assoc$finemap$L10[mol_assoc$finemap$L10 != 'None']
}

#################
# Drug repruposing
#################

tx<-list()

######
# Drug-specific
######

tx$drug<-list()

######
# ATC
######

tx$atc<-list()

#################
# Configuration
#################

config<-list()

config$repo<-system('git describe --tags', intern=T)
config$config<-config
config$gwas_list<-fread(gsub('gwas_list: ','', config[grepl('gwas_list: ',config)]))
