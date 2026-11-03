#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),  
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

# Read in the sumstats
ss<-fread(paste0(outdir,'/results/',opt$gwas,'/gwas_sumstat/',opt$gwas,'.cleaned.gz'))

# Read in COJO results
clumped_res<-NULL
for(i in 1:22){
  if(file.exists(paste0(outdir,'/results/',opt$gwas,'/clump/',opt$gwas,'_chr',i,'.clumped'))){
    tmp<-fread(paste0(outdir,'/results/',opt$gwas,'/clump/',opt$gwas,'_chr',i,'.clumped'))
    clumped_res<-rbind(clumped_res, tmp)
  }
}

# Make table with original sumstats but containing only independent associations from clumping
# Insert the COJO p-value
ss_subset<-ss[(ss$SNP %in% clumped_res$SNP),]

# Tidy table
ss_subset<-ss_subset[,names(ss_subset) %in% c('CHR','BP','SNP','A1','A2','OR','BETA','SE','P','INFO','FREQ','REF.FREQ','N'), with=F]
col_order<-match(c('CHR','BP','SNP','A1','A2','OR','BETA','SE','P','INFO','FREQ','REF.FREQ','N'), names(ss_subset))
col_order<-col_order[!is.na(col_order)]
ss_subset<-ss_subset[,col_order,with=F]

# Insert nearest gene information
biomart<-read.delim(paste0(resdir, '/data/biomart/biomart_genes_grch37.tsv'), stringsAsFactors=FALSE)
Genes<-biomart[,c('external_gene_name','chromosome_name','start_position','end_position')]
Genes<-Genes[!duplicated(Genes),]

# Protein-coding-only gene set for the companion NearestGeneProteinCoding column.
# Prefer gene_biotype on the biomart table (added by download_biomart.R); until
# that re-download lands, borrow biotype from the sibling gene_locations.tsv.
Genes_pc<-NULL
if('gene_biotype' %in% names(biomart)){
  Genes_pc<-biomart[biomart$gene_biotype == 'protein_coding',
                    c('external_gene_name','chromosome_name','start_position','end_position')]
} else {
  loc_path<-paste0(resdir, '/data/biomart/gene_locations.tsv')
  if(file.exists(loc_path)){
    loc<-read.delim(loc_path, stringsAsFactors=FALSE)
    loc<-loc[!is.na(loc$gene_biotype) & loc$gene_biotype == 'protein_coding',]
    Genes_pc<-loc[,c('external_gene_name','chromosome_name','start_position','end_position')]
  }
}
if(!is.null(Genes_pc)) Genes_pc<-Genes_pc[!duplicated(Genes_pc),]

window<-50000

# Nearest-gene string for one variant against a gene table, matching the legacy
# format: genes within `window` of BP, each "SYMBOL (±dist kb)" (bare symbol if
# the variant is inside the gene), sorted nearest-first; 'None' if nothing is in
# range. Shared by the all-gene and protein-coding-only annotations.
nearest_gene_string<-function(bp, chr, genes, window){
  gi<-genes[genes$start_position < (bp + window) & genes$end_position > (bp - window) & genes$chromosome_name == chr,]
  if(nrow(gi) == 0) return('None')
  gene_string<-NULL
  for(j in 1:nrow(gi)){
    if(bp > gi$start_position[j] & bp < gi$end_position[j]){
      gene_string<-rbind(gene_string, data.frame(ID=gi$external_gene_name[j], Dist=0, Pos=NA))
    }
    if(bp < gi$start_position[j]){
      gene_string<-rbind(gene_string, data.frame(ID=gi$external_gene_name[j], Dist=abs(bp - gi$start_position[j]), Pos='-'))
    }
    if(bp > gi$end_position[j]){
      gene_string<-rbind(gene_string, data.frame(ID=gi$external_gene_name[j], Dist=abs(bp - gi$end_position[j]), Pos='+'))
    }
  }
  gene_string$Text<-paste0(gene_string$ID, " (", gene_string$Pos, round(gene_string$Dist/1000,2),"kb)")
  gene_string$Text[gene_string$Dist == 0]<-as.character(gene_string$ID[gene_string$Dist == 0])
  gene_string<-gene_string[order(gene_string$Dist),]
  paste(gene_string$Text, collapse=', ')
}

# Always emit the columns, even when there are no clumped SNPs, so an empty
# ss_subset (GWAS with no genome-wide-significant clumps) still writes a CSV with
# the columns downstream consumers expect (e.g. plot_manhattan.R
# fread(..., select=c(..., "NearestGene"))). Non-empty runs overwrite each row.
ss_subset$NearestGene<-NA_character_
ss_subset$NearestGeneProteinCoding<-NA_character_

for(i in seq_len(nrow(ss_subset))){
  ss_subset$NearestGene[i]<-nearest_gene_string(ss_subset$BP[i], ss_subset$CHR[i], Genes, window)
  if(!is.null(Genes_pc))
    ss_subset$NearestGeneProteinCoding[i]<-nearest_gene_string(ss_subset$BP[i], ss_subset$CHR[i], Genes_pc, window)
}

# Write out results
write.csv(ss_subset, paste0(outdir,'/results/',opt$gwas,'/clump/',opt$gwas,'.GW.clump.clean.csv'), row.names=F, quote=T)



