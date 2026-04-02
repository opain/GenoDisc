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
library(susieR)
library(stringr)

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

# Read in the sumstats
ss<-fread(paste0(outdir,'/data/gwas_sumstat/',opt$gwas,'/',opt$gwas,'.cleaned.gz'))

# Read in gene locations
biomart<-read.delim('resources/data/biomart/biomart_genes_grch37.tsv', stringsAsFactors=FALSE)
Genes<-biomart[,c('external_gene_name','chromosome_name','start_position','end_position')]
Genes<-Genes[!duplicated(Genes),]

# Use 10kb window to define gene window
gene_window<-0
Genes$start_position<-Genes$start_position-gene_window
Genes$end_position<-Genes$end_position+gene_window

# List all finemapping results
finemap_files<-list.files(path=paste0(outdir,'/results/',opt$gwas,'/finemap/'), pattern='.rds')
finemap_files<-paste0(outdir,'/results/',opt$gwas,'/finemap/',finemap_files)

finemap_files_L10<-finemap_files[!grepl('L1.rds',finemap_files)]
finemap_files_L1<-finemap_files[grepl('L1.rds',finemap_files)]

# Summarise unrestricted L analysis
finemap_summary<-NULL
if(length(finemap_files_L10) > 0){
  for(i in 1:length(finemap_files_L10)){
    lead<-gsub('.*\\.','',gsub('.rds','',finemap_files_L10[i]))
    tmp<-readRDS(finemap_files_L10[i])
  
    plot_file<-gsub('.rds','.png',finemap_files_L10[i])
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
  
  if(is.null(finemap_summary)){
    finemap_summary<-data.frame(CHR = NA,
                                BP = NA,
                                SNP = NA,
                                cs=NA,
                                NSNP=NA,
                                TopPIP=NA,
                                Gene=NA)
  }
  
  # Write out results
  write.csv(finemap_summary, paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.GW.finemap.csv'), row.names=F, quote=T)
}

# Summarise L=1 restricted analyses
finemap_summary<-NULL
for(i in 1:length(finemap_files_L1)){
  lead<-gsub('.*\\.','',gsub('.L1.rds','',finemap_files_L1[i]))
  tmp<-readRDS(finemap_files_L1[i])
  
  plot_file<-gsub('.rds','.png',finemap_files_L1[i])
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
  
  if(is.null(finemap_summary)){
    finemap_summary<-data.frame(CHR = NA,
                                BP = NA,
                                SNP = NA,
                                cs=NA,
                                NSNP=NA,
                                TopPIP=NA,
                                Gene=NA)
  }
  
}

# Write out results
write.csv(finemap_summary, paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.GW.finemap.L1.csv'), row.names=F, quote=T)



