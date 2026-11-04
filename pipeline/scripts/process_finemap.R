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
library(susieR)
library(stringr)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

# Build a per-SNP PIP table for a whole fine-mapped locus (for the locus-plot
# baseline). Every SNP in the fit gets its marginal PIP; `cs` is the credible set
# it most strongly belongs to (largest per-effect alpha), or 0 if in none.
build_locus_pip <- function(tmp, lead, Lval, ss){
  pip<-tmp$pip
  snps<-names(pip)
  n<-length(pip)
  cs_assign<-integer(n)
  best_alpha<-rep(-Inf, n)
  s<-summary(tmp)$cs
  if(!is.null(s) && nrow(s) > 0){
    for(j in 1:nrow(s)){
      idx<-as.numeric(unlist(str_split(s$variable[j], ',')))
      idx<-idx[!is.na(idx)]
      if(length(idx) == 0) next
      a<-as.numeric(tmp$alpha[s$cs[j], idx])
      better<-a > best_alpha[idx]
      cs_assign[idx[better]]<-as.integer(s$cs[j])
      best_alpha[idx[better]]<-a[better]
    }
  }
  d<-data.table(locus=lead, SNP=snps, pip=as.numeric(pip), cs=cs_assign, L=as.integer(Lval))
  join_cols<-intersect(c('SNP','CHR','BP','P'), names(ss))
  d<-merge(d, ss[, join_cols, with=FALSE], by='SNP', all.x=TRUE, sort=FALSE)
  out_cols<-intersect(c('locus','SNP','CHR','BP','pip','cs','L','P'), names(d))
  d[, out_cols, with=FALSE]
}

# Read in config file
config<-readLines(opt$config_file)

# Identify outdir
outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])
resdir <- read_param(config = opt$config_file, param = 'resdir', return_obj = F)

# Read in the sumstats
ss<-fread(paste0(outdir,'/results/',opt$gwas,'/gwas_sumstat/',opt$gwas,'.cleaned.gz'))

# Read in gene locations
biomart<-read.delim(paste0(resdir, '/data/biomart/biomart_genes_grch37.tsv'), stringsAsFactors=FALSE)
Genes<-biomart[,c('external_gene_name','chromosome_name','start_position','end_position')]
Genes<-Genes[!duplicated(Genes),]

# Use 10kb window to define gene window
gene_window<-0
Genes$start_position<-Genes$start_position-gene_window
Genes$end_position<-Genes$end_position+gene_window

# List all finemapping results
finemap_files<-list.files(path=paste0(outdir,'/results/',opt$gwas,'/finemap/'), pattern='.rds')
if(length(finemap_files) > 0){
  finemap_files<-paste0(outdir,'/results/',opt$gwas,'/finemap/',finemap_files)
}
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

      lead_idx <- which(ss$SNP == lead)
      if(length(lead_idx) != 1){
        message("Skipping stale finemap result for lead ", lead, ": ", length(lead_idx), " matches in cleaned sumstats.")
        next
      }

      tmp_sum<-data.frame(CHR = ss$CHR[lead_idx],
                          BP = ss$BP[lead_idx],
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

# Summarise L=1 restricted analyses
finemap_summary<-NULL
if(length(finemap_files_L1) > 0){
  for(i in 1:length(finemap_files_L1)){
    lead<-gsub('.*\\.','',gsub('.L1.rds','',finemap_files_L1[i]))
    tmp<-readRDS(finemap_files_L1[i])

    plot_file<-gsub('.rds','.png',finemap_files_L1[i])
    png(plot_file, units='px', res=300, height=1500, width=1500)
    susie_plot(tmp, y="PIP", add_legend=T)
    dev.off()

    if(!is.null(summary(tmp)$cs)){

      lead_idx <- which(ss$SNP == lead)
      if(length(lead_idx) != 1){
        message("Skipping stale finemap result for lead ", lead, ": ", length(lead_idx), " matches in cleaned sumstats.")
        next
      }

      tmp_sum<-data.frame(CHR = ss$CHR[lead_idx],
                          BP = ss$BP[lead_idx],
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
write.csv(finemap_summary, paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.GW.finemap.L1.csv'), row.names=F, quote=T)

# ---- Per-SNP locus PIP tables (locus-plot association/fine-mapping baseline) ----
# Separate pass over the same fit files (same stale-lead guard as the summaries).
build_all_locus_pip <- function(files, Lval){
  out<-NULL
  for(f in files){
    base<-if(Lval == 1) gsub('.L1.rds','',f) else gsub('.rds','',f)
    lead<-gsub('.*\\.','',base)
    tmp<-tryCatch(readRDS(f), error=function(e) NULL)
    if(is.null(tmp) || is.null(summary(tmp)$cs)) next
    if(length(which(ss$SNP == lead)) != 1) next
    out<-rbind(out, build_locus_pip(tmp, lead, Lval, ss))
  }
  if(is.null(out)) out<-data.table(locus=character(), SNP=character(), CHR=integer(),
                                   BP=numeric(), pip=numeric(), cs=integer(), L=integer(), P=numeric())
  out
}
fwrite(build_all_locus_pip(finemap_files_L10, 10L),
       paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.GW.finemap.locus_pip.csv'))
fwrite(build_all_locus_pip(finemap_files_L1, 1L),
       paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.GW.finemap.locus_pip.L1.csv'))

# ---- Aggregate per-lead LD-to-lead tables (written by finemap.R) ----
ld_files<-list.files(path=paste0(outdir,'/results/',opt$gwas,'/finemap/'),
                     pattern='\\.ld_to_lead\\.csv$', full.names=TRUE)
ld_to_lead<-NULL
for(f in ld_files){
  d<-tryCatch(fread(f), error=function(e) NULL)
  if(!is.null(d) && nrow(d) > 0) ld_to_lead<-rbind(ld_to_lead, d, fill=TRUE)
}
if(is.null(ld_to_lead)) ld_to_lead<-data.table(locus=character(), SNP=character(), r2=numeric())
fwrite(ld_to_lead, paste0(outdir,'/results/',opt$gwas,'/finemap/',opt$gwas,'.GW.finemap.ld_to_lead.csv'))

