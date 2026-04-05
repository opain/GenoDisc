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
source('scripts/functions/utils_functions.R')

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

window<-50000

for(i in 1:nrow(ss_subset)){
  Genes_i<-Genes[Genes$start_position < (ss_subset$BP[i] + window) & Genes$end_position > (ss_subset$BP[i] - window) & Genes$chromosome_name == ss_subset$CHR[i],]
  if(nrow(Genes_i) != 0){
    gene_string<-NULL
    for(j in 1:nrow(Genes_i)){
      if(ss_subset$BP[i] > Genes_i$start_position[j] & ss_subset$BP[i] < Genes_i$end_position[j]){
        gene_string<-rbind(gene_string, data.frame(ID=Genes_i$external_gene_name[j],
                                                   Dist=0,
                                                   Pos=NA))
      }
      if(ss_subset$BP[i] < Genes_i$start_position[j]){
        gene_string<-rbind(gene_string, data.frame(ID=Genes_i$external_gene_name[j],
                                                   Dist=abs(ss_subset$BP[i] - Genes_i$start_position[j]),
                                                   Pos='-'))
      }
      if(ss_subset$BP[i] > Genes_i$end_position[j]){
        gene_string<-rbind(gene_string, data.frame(ID=Genes_i$external_gene_name[j],
                                                   Dist=abs(ss_subset$BP[i] - Genes_i$end_position[j]),
                                                   Pos='+'))
      }
    }
    gene_string$Text<-paste0(gene_string$ID, " (", gene_string$Pos, round(gene_string$Dist/1000,2),"kb)")
    gene_string$Text[gene_string$Dist == 0]<-as.character(gene_string$ID[gene_string$Dist == 0])
    gene_string<-gene_string[order(gene_string$Dist),]
    ss_subset$NearestGene[i]<-paste(gene_string$Text, collapse=', ')
  } else {
    ss_subset$NearestGene[i]<-'None'
  }
}

# Write out results
write.csv(ss_subset, paste0(outdir,'/results/',opt$gwas,'/clump/',opt$gwas,'.GW.clump.clean.csv'), row.names=F, quote=T)



