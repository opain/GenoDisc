#!/usr/bin/Rscript

library(data.table)
library(optparse)

option_list = list(
  make_option("--resdir", type="character", default="resources")
)
option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

# Download the GTEx v8 TPM data
dir.create(paste0(opt$resdir, '/data/gtex/'))
system(paste0('wget -O ', opt$resdir, '/data/gtex/GTEx_v8_median_tpm.gct.gz https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz'))
system(paste0('wget -O ', opt$resdir, '/data/gtex/GTEx_v8_samp_att.txt https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'))

#####
# QC GTEx data
#####

GTeX<-fread(paste0(opt$resdir, '/data/gtex/GTEx_v8_median_tpm.gct.gz'))
samp_att<-fread(paste0(opt$resdir, '/data/gtex/GTEx_v8_samp_att.txt'))

# Identify sample size per tissue
samp_att$SBJID<-sapply(strsplit(samp_att$SAMPID,'-'), function(x) paste0(x[1],'-',x[2]))

# Assign each tissue to broader tissue group
tissue_groups<-samp_att[,c('SMTSD','SMTS'), with=F]
names(tissue_groups)<-c('Tissue','Group')
tissue_groups<-tissue_groups[!duplicated(tissue_groups$Tissue),]
tissue_groups<-tissue_groups[order(tissue_groups$Group),]

# Calculate sample size for each tissue
tissue_groups$N<-NA
for(i in tissue_groups$Tissue){
    tissue_groups$N[tissue_groups$Tissue == i]<-length(unique(samp_att$SBJID[samp_att$SMTSD == i & samp_att$SMAFRZE == 'RNASEQ']))
}

tissue_groups<-tissue_groups[tissue_groups$Tissue %in% names(GTeX),]

# Update the tissue names to be less error prone (make a record of what they are originally for plotting and tables)
GTeX_tissues<-data.frame(Original=names(GTeX)[-1:-2])

names(GTeX)<-gsub(' ','_',names(GTeX))
names(GTeX)<-gsub(' ','_',names(GTeX))
names(GTeX)<-gsub("\\(",'',names(GTeX))
names(GTeX)<-gsub("\\)",'',names(GTeX))
names(GTeX)<-gsub("-",'',names(GTeX))
names(GTeX)<-gsub("__",'_',names(GTeX))

GTeX_tissues$new<-names(GTeX)[-1:-2]

tissue_groups<-merge(tissue_groups, GTeX_tissues, by.x='Tissue', by.y='Original')
fwrite(tissue_groups, paste0(opt$resdir, '/data/gtex/Tissue_labels.tsv'), sep='\t')

# First extract genes that have an average TPM > 1 in at least one tissue.
GTeX$tpm1<-apply(GTeX[,-1:-2], 1, function(x) any(x >= 1))
GTeX_tpm1<-GTeX[GTeX$tpm1 == T,]
GTeX_tpm1$tpm1<-NULL

# Windsorise TPM to 50
GTeX_tpm1_ids<-GTeX_tpm1[,1:2]
GTeX_tpm1_tissue<-GTeX_tpm1[,-1:-2]
GTeX_tpm1_tissue[GTeX_tpm1_tissue > 50] <-50
GTeX_tpm1_winsorised<-cbind(GTeX_tpm1_ids, GTeX_tpm1_tissue)

# log transform TPM with pseudocount 1
GTeX_tpm1_winsorised_ids<-GTeX_tpm1_winsorised[,1:2]
GTeX_tpm1_winsorised_tissue<-GTeX_tpm1_winsorised[,-1:-2]
GTeX_tpm1_winsorised_tissue<-log2(GTeX_tpm1_winsorised_tissue+1)
GTeX_tpm1_winsorised_pseudo1_log2<-cbind(GTeX_tpm1_winsorised_ids, GTeX_tpm1_winsorised_tissue)

##
# Create gene property analysis file for analysis of 54 specific tissues
##

# Calculate average across tissues
GTeX_tpm1_winsorised_pseudo1_log2_ids<-GTeX_tpm1_winsorised_pseudo1_log2[,1:2]
GTeX_tpm1_winsorised_pseudo1_log2_tissue<-GTeX_tpm1_winsorised_pseudo1_log2[,-1:-2]
GTeX_tpm1_winsorised_pseudo1_log2_tissue$Average<-rowMeans(GTeX_tpm1_winsorised_pseudo1_log2_tissue)
GTeX_tpm1_winsorised_pseudo1_log2_withAv<-cbind(GTeX_tpm1_winsorised_pseudo1_log2_ids, GTeX_tpm1_winsorised_pseudo1_log2_tissue)

##
# Create gene property analysis file for tissue groups
##

# Calculate weighted mean of TPM data across groups
weighted_mean <- function(values, weights) {
  # Check if the length of values and weights are equal
  if(length(values) != length(weights)) {
    stop("Length of values and weights must be the same")
  }

  # Calculate the weighted sum
  weighted_sum <- sum(values * weights)

  # Calculate the sum of weights
  sum_weights <- sum(weights)

  # Return the weighted mean
  return(weighted_sum / sum_weights)
}

GTeX_tpm1_winsorised_pseudo1_log2_ids<-GTeX_tpm1_winsorised_pseudo1_log2[,1:2]
GTeX_tpm1_winsorised_pseudo1_log2_tissue<-GTeX_tpm1_winsorised_pseudo1_log2[,-1:-2]

GTeX_tpm1_winsorised_pseudo1_log2_group<-NULL
for(i in unique(tissue_groups$Group)){
    tissue<-tissue_groups$new[tissue_groups$Group == i]
    n<-tissue_groups$N[tissue_groups$Group == i]

    GTeX_tpm1_winsorised_pseudo1_log2_tissue_i<-GTeX_tpm1_winsorised_pseudo1_log2_tissue[, tissue, with=F]
    GTeX_tpm1_winsorised_pseudo1_log2_tissue_i_mean<-apply(GTeX_tpm1_winsorised_pseudo1_log2_tissue_i, 1, function(x){
        weighted_mean(
            values=as.numeric(x),
            weights=n)
    })

    GTeX_tpm1_winsorised_pseudo1_log2_group<-cbind(GTeX_tpm1_winsorised_pseudo1_log2_group, GTeX_tpm1_winsorised_pseudo1_log2_tissue_i_mean)
}
GTeX_tpm1_winsorised_pseudo1_log2_group<-data.frame(GTeX_tpm1_winsorised_pseudo1_log2_group)
names(GTeX_tpm1_winsorised_pseudo1_log2_group)<-gsub(' ','_',unique(tissue_groups$Group))

# Calculate average across tissues
GTeX_tpm1_winsorised_pseudo1_log2_group$Average<-rowMeans(GTeX_tpm1_winsorised_pseudo1_log2_group)
GTeX_tpm1_winsorised_pseudo1_log2_group_withAv<-cbind(GTeX_tpm1_winsorised_pseudo1_log2_ids, GTeX_tpm1_winsorised_pseudo1_log2_group)

########
# Convert gene IDs to entrez IDs for analysis in MAGMA
########

MAGMA_id<-fread(paste0(opt$resdir, '/data/magma/NCBI37.3.gene.loc'))
MAGMA_id<-MAGMA_id[,c(1,6)]

# Tissue specific
GTeX_tpm1_winsorised_pseudo1_log2_withAv_MAGMA<-merge(MAGMA_id,GTeX_tpm1_winsorised_pseudo1_log2_withAv,by.x='V6',by.y='Description')
names(GTeX_tpm1_winsorised_pseudo1_log2_withAv_MAGMA)[1:3]<-c('Name','entrez','ensembl')
GTeX_tpm1_winsorised_pseudo1_log2_withAv_MAGMA$Name<-NULL
GTeX_tpm1_winsorised_pseudo1_log2_withAv_MAGMA$ensembl<-NULL
GTeX_tpm1_winsorised_pseudo1_log2_withAv_MAGMA<-GTeX_tpm1_winsorised_pseudo1_log2_withAv_MAGMA[!duplicated(GTeX_tpm1_winsorised_pseudo1_log2_withAv_MAGMA$entrez),]

fwrite(GTeX_tpm1_winsorised_pseudo1_log2_withAv_MAGMA, paste0(opt$resdir, '/data/gtex/GTEx_v8_tissue.tsv'), sep='\t')

# Tissue groups
GTeX_tpm1_winsorised_pseudo1_log2_group_withAv_MAGMA<-merge(MAGMA_id,GTeX_tpm1_winsorised_pseudo1_log2_group_withAv,by.x='V6',by.y='Description')
names(GTeX_tpm1_winsorised_pseudo1_log2_group_withAv_MAGMA)[1:3]<-c('Name','entrez','ensembl')
GTeX_tpm1_winsorised_pseudo1_log2_group_withAv_MAGMA$Name<-NULL
GTeX_tpm1_winsorised_pseudo1_log2_group_withAv_MAGMA$ensembl<-NULL
GTeX_tpm1_winsorised_pseudo1_log2_group_withAv_MAGMA<-GTeX_tpm1_winsorised_pseudo1_log2_group_withAv_MAGMA[!duplicated(GTeX_tpm1_winsorised_pseudo1_log2_group_withAv_MAGMA$entrez),]

fwrite(GTeX_tpm1_winsorised_pseudo1_log2_group_withAv_MAGMA, paste0(opt$resdir, '/data/gtex/GTEx_v8_group.tsv'), sep='\t')

