#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$gwas<-'SCHI02'

library(data.table)
library(susieR)
library(stringr)

# Read in the sumstats
ss<-fread(paste0('resources/data/gwas_sumstat/',opt$gwas,'/',opt$gwas,'.cleaned.gz'))

# Read in gene locations
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('external_gene_name','chromosome_name','start_position','end_position'), mart = ensembl)

# Use 10kb window to define gene window
gene_window<-10000
Genes$start_position<-Genes$start_position-gene_window
Genes$end_position<-Genes$end_position+gene_window

# List all finemapping results
finemap_files<-list.files(path=paste0('results/',opt$gwas,'/finemap/'), pattern='.rds')
finemap_files<-paste0('results/',opt$gwas,'/finemap/',finemap_files)

finemap_summary<-NULL
for(i in 1:length(finemap_files)){
  lead<-gsub('.*\\.','',gsub('.rds','',finemap_files[i]))
  tmp<-readRDS(finemap_files[i])

  plot_file<-gsub('.rds','.png',finemap_files[i])
  png(plot_file, units='px', res=300, height=1500, width=1500)
  susie_plot(tmp, y="PIP", add_legend=T)
  dev.off()
  
  if(!is.null(summary(tmp)$cs)){
    
    tmp_sum<-data.frame(CHR = ss$CHR[ss$SNP == lead],
                        BP = ss$BP[ss$SNP == lead],
                        SNP = lead,
                        summary(tmp)$cs,
                        NSNP=NA,
                        TopPIP=NA,
                        Gene=NA)
    
    for(j in 1:nrow(summary(tmp)$cs)){
      snp_index<-as.numeric(unlist(str_split(tmp_sum$variable[j], ',')))
      
      tmp_sum$NSNP[j]<-length(snp_index)
      tmp_sum$TopPIP[j]<-max(tmp$pip[snp_index])
      tmp_sum$variable[j]<-paste(names(tmp$pip)[snp_index], collapse=', ')
      
      ss_subset<-ss[ss$SNP %in% names(tmp$pip)[snp_index],]
      min_bp<-min(ss_subset$BP)
      max_bp<-max(ss_subset$BP)
      chr<-ss_subset$CHR[1]
      
      Genes_subset<-Genes[Genes$start_position < min_bp & Genes$end_position > max_bp & Genes$chromosome_name == chr,]
      if(nrow(Genes_subset) > 0){
        tmp_sum$Gene[j]<-paste(Genes_subset$external_gene_name, collapse=', ')
      }
    }
    
    finemap_summary<-rbind(finemap_summary, tmp_sum)
  }
}

# Write out results
write.csv(finemap_summary, paste0('results/',opt$gwas,'/finemap/',opt$gwas,'.GW.finemap.csv'), row.names=F, quote=T)



