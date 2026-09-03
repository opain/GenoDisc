#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAs [required]"),  
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

# Read in COJO results, and classify each chromosome's COJO status.
# rule cojo writes a per-chromosome ".cojo.status" file: "ok" (GCTA succeeded, with or
# without signals) or "reference_too_small" (GCTA aborted with "too many SNPs" because the
# independent signals outnumber the LD reference sample size). A chromosome that produced a
# .jma.cojo has genome-wide independent signals; "ok" with no .jma.cojo means no signal.
cojo_res<-NULL
cojo_status<-NULL
for(i in 1:22){
  jma<-paste0(outdir,'/results/',opt$gwas,'/cojo/',opt$gwas,'_chr',i,'.jma.cojo')
  sf<-paste0(outdir,'/results/',opt$gwas,'/cojo/',opt$gwas,'_chr',i,'.cojo.status')
  status_word<-if(file.exists(sf)) trimws(readLines(sf, warn=FALSE)[1]) else NA

  if(!is.na(status_word) && status_word == 'reference_too_small'){
    chr_status<-'reference_too_small'
  } else if(file.exists(jma)){
    chr_status<-'ok_with_signals'
  } else if(!is.na(status_word) && status_word == 'ok'){
    chr_status<-'ok_no_signals'
  } else {
    chr_status<-'not_run'
  }
  cojo_status<-rbind(cojo_status, data.frame(chr=i, status=chr_status, stringsAsFactors=FALSE))

  if(file.exists(jma)){
    tmp<-fread(jma)
    cojo_res<-rbind(cojo_res, tmp)
  }
}

# Always record the per-chromosome status so packaging/Shiny can report which chromosomes
# (if any) could not be analysed.
write.csv(cojo_status, paste0(outdir,'/results/',opt$gwas,'/cojo/',opt$gwas,'.GW.cojo.status.csv'), row.names=F, quote=F)

# Make table with original sumstats but containing only independent associations from COJO
# (insert the COJO joint p-value pJ). When no chromosome produced any independent signal
# (all failed, or none significant), emit a valid 0-row table so downstream never breaks.
if(is.null(cojo_res) || nrow(cojo_res) == 0){
  ss_subset<-ss[0]
  ss_subset$pJ<-numeric(0)
} else {
  ss_subset<-merge(ss, cojo_res[,c('SNP','pJ'),], by='SNP')
}

# Tidy table
ss_subset<-ss_subset[,names(ss_subset) %in% c('CHR','BP','SNP','A1','A2','OR','BETA','SE','P','pJ','INFO','FREQ','REF.FREQ','N'), with=F]
col_order<-match(c('CHR','BP','SNP','A1','A2','OR','BETA','SE','P','pJ','INFO','FREQ','REF.FREQ','N'), names(ss_subset))
col_order<-col_order[!is.na(col_order)]
ss_subset<-ss_subset[,col_order,with=F]

# Insert nearest gene information
biomart<-read.delim(paste0(resdir, '/data/biomart/biomart_genes_grch37.tsv'), stringsAsFactors=FALSE)
Genes<-biomart[,c('external_gene_name','chromosome_name','start_position','end_position')]
Genes<-Genes[!duplicated(Genes),]

window<-50000

for(i in seq_len(nrow(ss_subset))){
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

# On an empty table the nearest-gene loop above was skipped, so ensure the column still exists.
if(!('NearestGene' %in% names(ss_subset))) ss_subset$NearestGene<-character(0)

# Write out results
write.csv(ss_subset, paste0(outdir,'/results/',opt$gwas,'/cojo/',opt$gwas,'.GW.cojo.clean.csv'), row.names=F, quote=T)



