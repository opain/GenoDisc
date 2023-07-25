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

gwas_qc$cleaner_dat<-process_cleaner_log(outdir=outdir, gwas=opt$gwas)
gwas_qc$focus_dat<-process_focus_log(outdir=outdir, gwas=opt$gwas)
if(config[grepl('ldsc:',config)] == "ldsc: T"){
  gwas_qc$ldsc_dat<-process_ldsc_log(outdir=outdir, gwas=opt$gwas)
}

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

# Check whether TWAS was performed
twas_panel_psychencode_logical<-config[grepl('twas_panel_psychencode:',config)] == "twas_panel_psychencode: T"
twas_panel_fusion_logical<-config[grepl('twas_panel_fusion:',config)] == "twas_panel_fusion: T"
twas_logical<-any(twas_panel_psychencode_logical, twas_panel_fusion_logical)

if(twas_logical){
  mol_assoc$exp$fusion<-read_fusion_exp(config=config, gwas=opt$gwas, outdir=outdir)
}

###
# SMR
###

# Check whether SMR with expression data was performed
smr_expression_panel_psychencode_logical<-config[grepl('smr_expression_panel_psychencode:',config)] == "smr_expression_panel_psychencode: T"

smr_expression_panel_metabrain_basalganglia_logical<-config[grepl('smr_expression_panel_metabrain_basalganglia:',config)] == "smr_expression_panel_metabrain_basalganglia: T"
smr_expression_panel_metabrain_cerebellum_logical<-config[grepl('smr_expression_panel_metabrain_cerebellum:',config)] == "smr_expression_panel_metabrain_cerebellum: T"
smr_expression_panel_metabrain_cortex_logical<-config[grepl('smr_expression_panel_metabrain_cortex:',config)] == "smr_expression_panel_metabrain_cortex: T"
smr_expression_panel_metabrain_hippocampus_logical<-config[grepl('smr_expression_panel_metabrain_hippocampus:',config)] == "smr_expression_panel_metabrain_hippocampus: T"
smr_expression_panel_metabrain_spinalcord_logical<-config[grepl('smr_expression_panel_metabrain_spinalcord:',config)] == "smr_expression_panel_metabrain_spinalcord: T"

smr_expression_panel_eqtlgen_logical<-config[grepl('smr_expression_panel_eqtlgen:',config)] == "smr_expression_panel_eqtlgen: T"

metabrain_logical<-any(smr_expression_panel_metabrain_basalganglia_logical,
                        smr_expression_panel_metabrain_cerebellum_logical,
                        smr_expression_panel_metabrain_cortex_logical,
                        smr_expression_panel_metabrain_hippocampus_logical,
                        smr_expression_panel_metabrain_spinalcord_logical)

smr_expression_logical<-any(smr_expression_panel_psychencode_logical,
                         metabrain_logical,
                         smr_expression_panel_eqtlgen_logical)
if(smr_expression_logical){
  mol_assoc$exp$smr<-read_smr_exp(config=config, gwas=opt$gwas, outdir=outdir)
}

######
# Protein
######

mol_assoc$protein<-list()


