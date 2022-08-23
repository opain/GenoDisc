#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--lincs_level5_path", action="store", default=NA, type='character',
              help="Path to LINCS Level 5 data [required]"),
  make_option("--lincs_siginfo_path", action="store", default=NA, type='character',
              help="Path to LINCS signiture information data [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(cmapR)
library(biomaRt)

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
# signinfo<-signinfo[(signinfo$pert_type %in% c('trt_cp','trt_lig','trt_oe','trt_sh','trt_xpr')),]

# Retain only trt_cp
signinfo<-signinfo[signinfo$pert_type == 'trt_cp',]

# Remove cell types with < 2000 signitures (chosen to remove unknown cell types whilst retaining NPC)
cell_count<-melt(table(signinfo$cell_iname))
signinfo<-signinfo[signinfo$cell_iname %in% as.character(cell_count$Var1[cell_count$value > 2000]),]

# Retain only drug with 10uM and 24hr
signinfo<-signinfo[signinfo$pert_itime == '24 h' & signinfo$pert_idose == '10 uM',]

# Read in the lincs data
lincs = parse_gctx(opt$lincs_level5_path, cid=signinfo$sig_id)@mat
lincs<-data.table(entrez_id=rownames(lincs),
                 lincs)

lincs<-lincs[!duplicated(lincs$entrez_id),]

# Insert entrez IDs into twas results
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
biomartCacheClear()
Genes<-getBM(attributes=c('ensembl_gene_id','entrezgene_id'), mart = ensembl)

Genes<-Genes[!duplicated(Genes$entrezgene_id),]
Genes<-Genes[!duplicated(Genes$ensembl_gene_id),]

lincs<-merge(Genes,lincs, by.x='entrezgene_id', by.y='entrez_id')
lincs$entrezgene_id<-NULL

lincs<-data.table(lincs)

########
# Split data into ten batches
########

n_drug<-ncol(lincs)-1
n_batch<-10

batches<-ceiling(2:ncol(lincs)/ceiling(n_drug/n_batch))

lincs_IDs<-lincs[,1]
for(i in 1:n_batch){
  saveRDS(cbind(lincs_IDs, lincs[,which(batches == i) + 1, with=F]), file = paste0("resources/data/lincs/lincs_core_subset_batch_",i,".rds"))
}

file.create("resources/data/lincs/subset_lincs.done")

