#!/usr/bin/Rscript
# Save start time
start.time <- Sys.time()

# Load optparse package
library("optparse")

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="config file [required]")
)

# Parse parameters from command line
opt = parse_args(OptionParser(option_list=option_list))

# Load required packages
library(data.table)

# Read in config file
config_file<-readLines(opt$config_file)

# Idenitfy outdir
outdir<-gsub('outdir: ','',config_file[grepl('^outdir:',config_file)])

sink(file = paste(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.log',sep=''), append = F)
cat(
  '#################################################################
# magma_tissue_conditional.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

# Read in the MAGMA gene property enrichment results
property_enrich<-fread(cmd=paste0("grep -v '^#' ",outdir,"/results/",opt$gwas,'/magma/magma_tissue_spec.gsa.out'))

# Insert FULL_NAME column if not present
if(all(names(property_enrich) != 'FULL_NAME')){
    property_enrich$FULL_NAME<-property_enrich$VARIABLE
}

# Select FDR significant properties
property_enrich$P.FDR<-p.adjust(property_enrich$P, method = 'fdr')
property_enrich<-property_enrich[property_enrich$P.FDR <= 0.05,]

sink(file = paste(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.log',sep=''), append = T)
cat(paste0(nrow(property_enrich)," properties are FDR significant\n"))
sink()

# If more than 1 sig set, perform conditional analysis
if(nrow(property_enrich) > 1){
    
    sink(file = paste(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.log',sep=''), append = T)
    cat(paste0("Performing conditional analysis...\n"))
    sink()

    # Create object indicating tmpdir
    tmp_folder<-tempdir()
    
    # Read in property file
    property_annot<-fread("resources/data/gtex/GTEx_v8_tissue.tsv")

    # Subset property file to contain enriched properties
    property_annot<-property_annot[,names(property_annot) %in% c('entrez',property_enrich$FULL_NAME,'Average'), with=F]
    write.table(property_annot, paste0(tmp_folder,"/sig_property.txt"), col.names=T, row.names=F, quote=F)

    # Sort results by p-value
    property_enrich<-property_enrich[order(property_enrich$P),]
    
    # Now condition each set on the most significant properties until all are independently significant
    i<-1
    
    property_indep <- property_enrich
    chisq_mat<-
    while(1){
        if(nrow(property_indep) <= i){
            break
        }
        
        property_i<-property_indep$FULL_NAME[i]

        log<-system(paste0(
            "resources/software/magma/magma",
            " --gene-results ",outdir,"/results/",opt$gwas,"/magma/magma_gene_level.genes.raw",
            " --gene-covar ",tmp_folder,"/sig_property.txt",
            " --model direction-covar=greater condition-hide=",paste(c(property_i,'Average'),collapse=','),
            " --out ",tmp_folder,"/res"
        ), intern = T)
        
        if(any(grepl('ERROR - running gene-level regression: could not invert design matrix of conditioned-on variables; variables are collinear with each other', log))){
            print(log)
            print('ERROR: There was too much multicolinearity between sets.')
            q()
        }
        
        # Read in the results
        cond_res<-fread(cmd=paste0("grep -v '^#' ",tmp_folder,'/res.gsa.out'))
        
        # Insert FULL_NAME column if not present
        if(all(names(cond_res) != 'FULL_NAME')){
            cond_res$FULL_NAME<-cond_res$VARIABLE
        }
        
        # Remove properties already excluded
        cond_res<-cond_res[cond_res$FULL_NAME %in% property_indep$FULL_NAME,]
        
        # Calculate % explained
        for(j in cond_res$FULL_NAME){
            cond_chisq<-(cond_res$BETA[cond_res$FULL_NAME == j]/cond_res$SE[cond_res$FULL_NAME == j])^2
            marg_chisq<-(property_enrich$BETA[property_enrich$FULL_NAME == j]/property_enrich$SE[property_enrich$FULL_NAME == j])^2
            lead_chisq<-(property_enrich$BETA[property_enrich$FULL_NAME == property_i]/property_enrich$SE[property_enrich$FULL_NAME == property_i])^2
            cond_res$VAR.EXP[cond_res$FULL_NAME == j] <- 1 - (cond_chisq / marg_chisq)
            cond_res$R2_ratio[cond_res$FULL_NAME == j] <- marg_chisq / lead_chisq
        }
        
        sink(file = paste(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.log',sep=''), append = T)
        cat(paste0("Done!\n"))
        cat(paste0(nrow(property_indep), " independent properties remain.\n"))
        sink()

        # Retain tissues that:
        # - Are colinear with the lead tissue (returns NA)
        # - Are similar marginal significance (R2_ratio > 0.95)
        # - Are still significant (P<0.05)
        cond_res<-cond_res[
            is.na(cond_res$VAR.EXP) | 
            cond_res$R2_ratio > 0.95 | 
            cond_res$P < 0.05
        ,]

        property_indep<-property_indep[property_indep$FULL_NAME %in% c(property_i, cond_res$FULL_NAME),]

        i<-i+1
    }
    
    # Save file listing significant and independent property
    write.table(property_indep$FULL_NAME, paste0(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.indep.txt'), row.names=F, col.names=F, quote=F)
                
    sink(file = paste(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.log',sep=''), append = T)
    cat(paste0("Done!\n"))
    cat(paste0(nrow(property_indep), " independent properties remain.\n"))
    sink()

}

# If 1 sig set, no conditional analysis required
if(nrow(property_enrich) == 1){
    sink(file = paste(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.log',sep=''), append = T)
    cat(paste0("No conditional analysis required.\n"))
    sink()

    # Save file listing significant and independent properties
    write.table(property_enrich$FULL_NAME, paste0(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.indep.txt'), row.names=F, col.names=F, quote=F)
}

if(nrow(property_enrich) == 0){
    sink(file = paste(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.log',sep=''), append = T)
    cat(paste0("No conditional analysis required.\n"))
    sink()
}
  
end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(outdir,"/results/",opt$gwas,'/magma/magma_tissue_conditional.log',sep=''), append = T)
cat(paste0('Analysis finished at ',as.character(end.time),'\n'))
cat(paste0('Analysis duration was ',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n'))
sink()

