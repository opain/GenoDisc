#!/usr/bin/Rscript
library("optparse")

option_list <- list(
  make_option("--config", action = "store", default = NA, type = "character",
              help = "Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))
source_all(file.path(opt$pipeline_dir, 'scripts', 'functions'))

# Read in config: merge default config with user config (user takes priority)
library(yaml)
default_config <- read_yaml(file.path(opt$pipeline_dir, 'config.yaml'))
user_config <- read_yaml(opt$config)
merged_config <- default_config
merged_config[names(user_config)] <- user_config

# Convert to readLines-style character vector for downstream compatibility
config <- vapply(names(merged_config), function(k) {
  v <- merged_config[[k]]
  if (is.null(v) || (length(v) == 1 && is.na(v))) {
    paste0(k, ": NA")
  } else if (length(v) > 1) {
    paste0(k, ": [", paste0("'", v, "'", collapse = ","), "]")
  } else {
    paste0(k, ": ", v)
  }
}, character(1), USE.NAMES = FALSE)

# Read in config parameters
gwas_list<-read_param(config = opt$config, param = 'gwas_list', return_obj = T)
outdir <- read_param(config = opt$config, param = 'outdir', return_obj = F)

output<-list()
for(gwas_i in gwas_list$name){

  ##################
  # GWAS QC
  ##################

  gwas_qc<-list()

  gwas_qc$cleaner_dat<-process_cleaner_log(config=opt$config, gwas=gwas_i)
  gwas_qc$focus_dat<-process_focus_log(config=opt$config, gwas=gwas_i)
  gwas_qc$ldsc_dat<-process_ldsc_log(config=opt$config, gwas=gwas_i)

  # Read MAF plot as base64 (if it exists)
  maf_plot_path <- paste0(outdir, '/results/', gwas_i, '/gwas_sumstat/', gwas_i, '.cleaned.MAF_plot.png')
  if (file.exists(maf_plot_path)) {
    gwas_qc$maf_plot_base64 <- base64enc::base64encode(maf_plot_path)
  } else {
    gwas_qc$maf_plot_base64 <- NULL
  }

  ##################
  # SNP associations
  ##################

  snp_assoc<-list()

  if(read_param(config = opt$config, param = 'clump', return_obj = F) == "T"){
    snp_assoc$clump<-fread(paste0(outdir,'/results/',gwas_i,'/clump/',gwas_i,'.GW.clump.clean.csv'))
  } else {
    snp_assoc$clump<-NULL
  }

  if(read_param(config = opt$config, param = 'cojo', return_obj = F) == "T"){
    snp_assoc$cojo<-fread(paste0(outdir,'/results/',gwas_i,'/cojo/',gwas_i,'.GW.cojo.clean.csv'))
    snp_assoc$cojo<-snp_assoc$cojo[order(snp_assoc$cojo$CHR, snp_assoc$cojo$BP),]
  } else {
    snp_assoc$cojo<-NULL
  }

  if(read_param(config = opt$config, param = 'finemap', return_obj = F) == "T"){
    snp_assoc$susie<-list()
    snp_assoc$susie$L10<-process_susie(outdir=outdir, gwas=gwas_i, L=10)
    snp_assoc$susie$L1<-process_susie(outdir=outdir, gwas=gwas_i, L=1)
  } else {
    snp_assoc$susie<-NULL
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

  mol_assoc$exp$fusion<-read_fusion_exp(config=opt$config, gwas=gwas_i)

  ###
  # SMR
  ###

  mol_assoc$exp$smr<-read_smr_exp(config=opt$config, gwas=gwas_i)

  ######
  # Protein
  ######

  mol_assoc$protein<-list()

  ###
  # FUSION
  ###

  mol_assoc$protein$fusion<-read_fusion_protein(config=opt$config, gwas=gwas_i)

  ###
  # SMR
  ###

  mol_assoc$protein$smr<-read_smr_protein(config=opt$config, gwas=gwas_i)

  ######
  # MAGMA
  ######

  mol_assoc$magma<-read_magma_gene(config=opt$config, gwas=gwas_i)

  ######
  # Nearest
  ######

  mol_assoc$nearest<-list()

  if(read_param(config = opt$config, param = 'clump', return_obj = F) == "T"){
    mol_assoc$nearest$clump<-identify_nearest(snp_assoc$clump$NearestGene)
  }

  if(read_param(config = opt$config, param = 'cojo', return_obj = F) == "T"){
    mol_assoc$nearest$cojo<-identify_nearest(snp_assoc$cojo$NearestGene)
  }

  ######
  # Finemapping
  ######

  mol_assoc$finemap<-list()

  if(read_param(config = opt$config, param = 'finemap', return_obj = F) == "T"){
    mol_assoc$finemap$L1<-unlist(strsplit(snp_assoc$susie$L1$Gene, ', '))
    mol_assoc$finemap$L1<-mol_assoc$finemap$L1[mol_assoc$finemap$L1 != 'None']
    mol_assoc$finemap$L10<-unlist(strsplit(snp_assoc$susie$L10$Gene, ', '))
    mol_assoc$finemap$L10<-mol_assoc$finemap$L10[mol_assoc$finemap$L10 != 'None']
  } else {
    mol_assoc$finemap<-NULL
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

  tx$drug$magma<-read_magma_drug(config=opt$config, gwas=gwas_i)

  ###
  # GCSC
  ###

  tx$drug$gcsc <- read_gcsc(config = opt$config, gwas = gwas_i)

  ###
  # TWAS-GSEA
  ###

  tx$drug$twas_gsea<-read_twas_gsea_drug(config=opt$config, gwas=gwas_i, mode='directional')
  tx$drug$twas_gsea_nondir<-read_twas_gsea_drug(config=opt$config, gwas=gwas_i, mode='nondirectional')

  ######
  # ATC
  ######

  tx$atc<-list()

  ###
  # MAGMA
  ###

  tx$atc$magma<-read_magma_drug_atc(config=opt$config, gwas=gwas_i)

  ###
  # GCSC
  ###

  tx$atc$gcsc<-read_gcsc_atc(config=opt$config, gwas=gwas_i)

  ###
  # TWAS-GSEA
  ###

  tx$atc$twas_gsea<-read_twas_gsea_atc(config=opt$config, gwas=gwas_i, mode='directional')
  tx$atc$twas_gsea_nondir<-read_twas_gsea_atc(config=opt$config, gwas=gwas_i, mode='nondirectional')

  ######
  # CMAP TWAS-GSEA (per-signature drug + per-MOA enrichment)
  ######

  tx$cmap<-list()
  tx$cmap$drug<-read_twas_gsea_cmap_drug(config=opt$config, gwas=gwas_i)
  tx$cmap$moa<-read_twas_gsea_cmap_moa(config=opt$config, gwas=gwas_i)

  #################
  # Tissue Enrichment
  #################

  tissue<-list()

  ######
  # Tissue-specific
  ######

  tissue$specific<-read_magma_tissue(config=opt$config, gwas=gwas_i, type='specific')

  ################
  # Package results
  ################

  output[[gwas_i]]<-list( gwas_qc=gwas_qc,
                          snp_assoc=snp_assoc,
                          mol_assoc=mol_assoc,
                          tx=tx,
                          tissue=tissue)

}

#################
# Configuration
#################

output$configuration<-list()

output$configuration$repo<-list()
output$configuration$repo$remote<-gsub('.*@','',gsub(' .*','',system('git remote -v', intern=T)[1])) # nolint # nolint: line_length_linter.
output$configuration$repo$branch<-gsub('On branch ','', system('git status', intern=T)[1])
output$configuration$repo$commit<-system('git describe --tags --always', intern=T)
output$configuration$config<-config
output$configuration$gwas_list<-gwas_list

################
# Save results as .RDS
################

saveRDS(output, file = paste0(outdir,'/results/results_package.rds'))

