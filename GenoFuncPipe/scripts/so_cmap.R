#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--lincs_level5_path", action="store", default=NA, type='character',
              help="Path to LINCS Level 5 data [required]"),
  make_option("--lincs_siginfo_path", action="store", default=NA, type='character',
              help="Path to LINCS signiture information data [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$lincs_level5_path<-'/mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/CMAP/2020_beta_01022022/level5_beta_all_n1201944x12328.gctx'
opt$lincs_siginfo_path<-'/mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/CMAP/2020_beta_01022022/siginfo_beta.txt'

library(dplyr)
library(Hmisc)
library(data.table)
library(cmapR) # Installation of these 
library(mygene)

# Installation of these bioconductor packages takes ages. Consider removing dependency.

# Read in siginfo file
signinfo<-fread(opt$lincs_siginfo_path)

# Extract that pass QC
signinfo<-signinfo[signinfo$qc_pass == 1,]

# Extract that are high techincal quality
signinfo<-signinfo[signinfo$is_hiq == 1,]

# Retain pert_type, pert_idose, pert_itime combinations with the lowest batch effect
signinfo<-signinfo[order(abs(signinfo$batch_effect_tstat)),]
signinfo<-signinfo[!duplicated(paste0(signinfo$pert_mfc_id,'-',signinfo$pert_itime,'-',signinfo$pert_idose,'-',signinfo$cell_iname))]

# Remove control pertubagens
signinfo<-signinfo[!grepl('ctl', signinfo$pert_type),]

# Remove trt_cp trt_lig  trt_oe  trt_sh  trt_xpr pertubagens
signinfo<-signinfo[(signinfo$pert_type %in% c('trt_cp','trt_lig','trt_oe','trt_sh','trt_xpr')),]

# Remove cell types with < 2000 signitures (chosen to remove unknown cell types whilst retaining NPC)
cell_count<-melt(table(signinfo$cell_iname))
signinfo<-signinfo[signinfo$cell_iname %in% as.character(cell_count$Var1[cell_count$value > 2000]),]

tmp<-signinfo[signinfo$cmap_name == 'citalopram',]

# Subset pert_type == trt_cp
signinfo_trt_cp<-signinfo[signinfo$pert_type == 'trt_cp',]

# Subset pert_iname=="citalopram"
signinfo_trt_cp_citalopram<-signinfo_trt_cp[signinfo_trt_cp$cmap_name == 'citalopram',]

table(signinfo_trt_cp_citalopram$pert_idose)
table(signinfo_trt_cp_citalopram$pert_itime)
table(signinfo_trt_cp_citalopram$cell_iname)

# Extract data for NPC/NEU cell lines
signinfo_trt_cp_citalopram_NEU_NPC<-signinfo_trt_cp_citalopram[  signinfo_trt_cp_citalopram$cell_iname == 'NEU' |
                                                                 signinfo_trt_cp_citalopram$cell_iname == 'NPC',]


table(signinfo_trt_cp_citalopram_NEU_NPC$pert_idose)
table(signinfo_trt_cp_citalopram_NEU_NPC$pert_itime)

# The most data is available for pert_itime == 24 hours
signinfo_trt_cp_citalopram_NEU_NPC[signinfo_trt_cp_citalopram_NEU_NPC$pert_idose == '0.74 uM',]