#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(  
  make_option("--config", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
source('scripts/functions/package_results_functions.R')

# Read in config file
config<-readLines(opt$config)

# Read in GWAS list
gwas_list<-fread(gsub('gwas_list: ','', config[grepl('gwas_list: ',config)]))

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

output<-list()
for(gwas_i in gwas_list$name){

  ##################
  # GWAS QC
  ##################

  gwas_qc<-list()

  gwas_qc$cleaner_dat<-process_cleaner_log(config=config, gwas=gwas_i)
  gwas_qc$focus_dat<-process_focus_log(config=config, gwas=gwas_i)
  gwas_qc$ldsc_dat<-process_ldsc_log(config=config, gwas=gwas_i)

  ##################
  # SNP associations
  ##################

  snp_assoc<-list()

  if(config[grepl('clump:',config)] == "clump: T"){
    snp_assoc$clump<-fread(paste0(outdir,'/results/',gwas_i,'/clump/',gwas_i,'.GW.clump.clean.csv'))
  }

  if(config[grepl('cojo:',config)] == "cojo: T"){
    snp_assoc$cojo<-fread(paste0(outdir,'/results/',gwas_i,'/cojo/',gwas_i,'.GW.cojo.clean.csv'))
  }

  if(config[grepl('finemap:',config)] == "finemap: T"){
    snp_assoc$susie<-list()
    snp_assoc$susie$L10<-process_susie(outdir=outdir, gwas=gwas_i, L=10)
    snp_assoc$susie$L1<-process_susie(outdir=outdir, gwas=gwas_i, L=1)
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

  mol_assoc$exp$fusion<-read_fusion_exp(config=config, gwas=gwas_i)

  ###
  # SMR
  ###

  mol_assoc$exp$smr<-read_smr_exp(config=config, gwas=gwas_i)

  ######
  # Protein
  ######

  mol_assoc$protein<-list()

  ###
  # FUSION
  ###

  mol_assoc$protein$fusion<-read_fusion_protein(config=config, gwas=gwas_i)

  ###
  # SMR
  ###

  mol_assoc$protein$smr<-read_smr_protein(config=config, gwas=gwas_i)

  ######
  # MAGMA
  ######

  mol_assoc$magma<-read_magma_gene(config=config, gwas=gwas_i)

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

  ###
  # MAGMA
  ###

  tx$drug$magma<-read_magma_drug(config=config, gwas=gwas_i)

  ###
  # GCSC
  ###

  tx$drug$gcsc<-read_gcsc(config=config, gwas=gwas_i)

  ###
  # TWAS-GSEA
  ###

  tx$drug$twas_gsea<-read_twas_gsea_drug(config=config, gwas=gwas_i)

  ######
  # ATC
  ######

  tx$atc<-list()

  ###
  # MAGMA
  ###

  tx$atc$magma<-read_magma_drug_atc(config=config, gwas=gwas_i)

  ###
  # GCSC
  ###

  tx$atc$gcsc<-read_gcsc_atc(config=config, gwas=gwas_i)

  ###
  # TWAS-GSEA
  ###

  tx$atc$twas_gsea<-read_twas_gsea_atc(config=config, gwas=gwas_i)

  ################
  # Package results
  ################

  output[[gwas_i]]<-list( gwas_qc=gwas_qc,
                          snp_assoc=snp_assoc,
                          mol_assoc=mol_assoc,
                          tx=tx)

}

#################
# Configuration
#################

output$configuration<-list()

output$configuration$repo<-list()
output$configuration$repo$remote<-gsub('.*@','',gsub(' .*','',system('git remote -v', intern=T)[1]))
output$configuration$repo$branch<-gsub('On branch ','', system('git status', intern=T)[1])
output$configuration$repo$commit<-system('git describe --tags --always', intern=T)
output$configuration$config<-config
output$configuration$gwas_list<-fread(gsub('gwas_list: ','', config[grepl('gwas_list: ',config)]))

################
# Save results as .RDS
################

saveRDS(output, file = paste0(outdir,'/results/results_package.rds'))

