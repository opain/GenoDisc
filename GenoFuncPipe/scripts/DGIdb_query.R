#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="config_file location [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(rjson)

# Read in config file
config<-readLines(opt$config_file)

# Determine which analyses were requested
twas_panel_psychencode_logical<-config[grepl('twas_panel_psychencode:',config)] == "twas_panel_psychencode: T"
twas_panel_fusion_logical<-config[grepl('twas_panel_fusion:',config)] == "twas_panel_fusion: T"
twas_logical<-any(twas_panel_psychencode_logical, twas_panel_fusion_logical)

smr_expression_panel_psychencode_logical<-config[grepl('smr_expression_panel_psychencode:',config)] == "smr_expression_panel_psychencode: T"

smr_expression_panel_metabrain_basalganglia_logical<-config[grepl('smr_expression_panel_metabrain_basalganglia:',config)] == "smr_expression_panel_metabrain_basalganglia: T"
smr_expression_panel_metabrain_cerebellum_logical<-config[grepl('smr_expression_panel_metabrain_cerebellum:',config)] == "smr_expression_panel_metabrain_cerebellum: T"
smr_expression_panel_metabrain_cortex_logical<-config[grepl('smr_expression_panel_metabrain_cortex:',config)] == "smr_expression_panel_metabrain_cortex: T"
smr_expression_panel_metabrain_hippocampus_logical<-config[grepl('smr_expression_panel_metabrain_hippocampus:',config)] == "smr_expression_panel_metabrain_hippocampus: T"
smr_expression_panel_metabrain_spinalcord_logical<-config[grepl('smr_expression_panel_metabrain_spinalcord:',config)] == "smr_expression_panel_metabrain_spinalcord: T"

smr_expression_panel_eqtlgen_logical<-config[grepl('smr_expression_panel_eqtlgen:',config)] == "smr_expression_panel_eqtlgen: T"

metabrain_logical<-any(smr_expression_panel_metabrain_basalganglia_logical,
                       smr_expression_panel_metabrain_cerebellum_logical,
                       smr_expression_panel_metabrain_cortex_logical,
                       smr_expression_panel_metabrain_hippocampus_logical,
                       smr_expression_panel_metabrain_spinalcord_logical)

pwas_panel_rosmap_logical<-config[grepl('pwas_panel_rosmap:',config)] == "pwas_panel_rosmap: T"
pwas_panel_banner_logical<-config[grepl('pwas_panel_banner:',config)] == "pwas_panel_banner: T"

smr_protein_panel_rosmap_logical<-config[grepl('smr_protein_panel_rosmap:',config)] == "smr_protein_panel_rosmap: T"

###############
# Read in gene associations with directional information (i.e. PWAS/TWAS/SMR)
###############

gene_assoc<-NULL
if(twas_logical){
  twas_all<-fread(paste0('results/',opt$gwas,'/twas/',opt$gwas,'_twas_GW_clean.txt.gz'))
  twas_all$TWAS.P.FDR<-p.adjust(twas_all$TWAS.P, method='fdr')
  twas_all<-twas_all[twas_all$TWAS.P.FDR < 0.05 & (twas_all$COLOC.PP4-twas_all$COLOC.PP3)/twas_all$COLOC.PP4 > 0.8,]
  twas_all<-twas_all[!grepl('SPLIC',twas_all$PANEL),]
  
  gene_assoc<-rbind(gene_assoc, data.frame(ID=twas_all$ID,
                                           DIR=sign(twas_all$TWAS.Z)))
}

if(smr_expression_panel_psychencode_logical){
  smr_psychencode_files<-list.files(path=paste0('results/',opt$gwas,'/smr/psychencode/'), pattern=paste0(opt$gwas,'_smr_psychencode_chr'))
  smr_psychencode_files<-smr_psychencode_files[grepl('.smr$', smr_psychencode_files)]
  
  smr_psychencode<-NULL
  for(i in smr_psychencode_files){
    smr_psychencode<-rbind(smr_psychencode, fread(paste0('results/',opt$gwas,'/smr/psychencode/',i)))
  }
  
  smr_psychencode$p_SMR.FDR<-p.adjust(smr_psychencode$p_SMR, method = 'fdr')
  
  smr_psychencode<-smr_psychencode[smr_psychencode$p_SMR.FDR < 0.05 & smr_psychencode$p_HEIDI > 0.05,]
  
  gene_assoc<-rbind(gene_assoc, data.frame(ID=smr_psychencode$Gene,
                                           DIR=sign(smr_psychencode$b_SMR)))
}

if(metabrain_logical){
  smr_metabrain_all<-fread(paste0('results/',opt$gwas,'/smr/metabrain/',opt$gwas,'_smr_metabrain_GW.txt.gz'))
  
  smr_metabrain_all$p_SMR.FDR<-p.adjust(smr_metabrain_all$p_SMR, method = 'fdr')
  
  smr_metabrain_all<-smr_metabrain_all[smr_metabrain_all$p_SMR.FDR < 0.05 & smr_metabrain_all$p_HEIDI > 0.05,]
  
  gene_assoc<-rbind(gene_assoc, data.frame(ID=smr_metabrain_all$external_gene_name,
                                           DIR=sign(smr_metabrain_all$b_SMR)))
  
}

if(smr_expression_panel_eqtlgen_logical){
  smr_eqtlgen<-fread(paste0('results/',opt$gwas,'/smr/eqtlgen/',opt$gwas,'_smr_eqtlgen_GW.txt.gz'))
  
  smr_eqtlgen$p_SMR.FDR<-p.adjust(smr_eqtlgen$p_SMR, method = 'fdr')
  
  smr_eqtlgen<-smr_eqtlgen[smr_eqtlgen$p_SMR.FDR < 0.05 & smr_eqtlgen$p_HEIDI > 0.05,]
  
  gene_assoc<-rbind(gene_assoc, data.frame(ID=smr_eqtlgen$external_gene_name,
                                           DIR=sign(smr_eqtlgen$b_SMR)))
  
}

if(pwas_panel_rosmap_logical){
  
  pwas_rosmap_files<-list.files(path=paste0('results/',opt$gwas,'/pwas/rosmap/'), pattern=paste0(opt$gwas,'_pwas_rosmap_chr'))
  
  pwas_rosmap<-NULL
  for(i in pwas_rosmap_files){
    pwas_rosmap<-rbind(pwas_rosmap, fread(paste0('results/',opt$gwas,'/pwas/rosmap/',i)))
  }
  
  pwas_rosmap$TWAS.P.FDR<-p.adjust(pwas_rosmap$TWAS.P, method = 'fdr')
  
  pwas_rosmap$`Ensembl ID`<-gsub('\\..*','',pwas_rosmap$ID)
  pwas_rosmap$`Gene Symbol`<-gsub('.*\\.','',pwas_rosmap$ID)
  
  pwas_rosmap<-pwas_rosmap[,c('CHR','P0','P1','Ensembl ID','Gene Symbol','TWAS.Z','TWAS.P','TWAS.P.FDR','COLOC.PP3','COLOC.PP4'), with=T]
  
  pwas_rosmap<-pwas_rosmap[pwas_rosmap$TWAS.P.FDR < 0.05 & (pwas_rosmap$COLOC.PP4-pwas_rosmap$COLOC.PP3)/pwas_rosmap$COLOC.PP4 > 0.8,]

  gene_assoc<-rbind(gene_assoc, data.frame(ID=pwas_rosmap$`Gene Symbol`,
                                           DIR=sign(pwas_rosmap$TWAS.Z)))
  
  
}

if(pwas_panel_banner_logical){
  
  pwas_banner_files<-list.files(path=paste0('results/',opt$gwas,'/pwas/banner/'), pattern=paste0(opt$gwas,'_pwas_banner_chr'))
  
  pwas_banner<-NULL
  for(i in pwas_banner_files){
    pwas_banner<-rbind(pwas_banner, fread(paste0('results/',opt$gwas,'/pwas/banner/',i)))
  }
  
  pwas_banner$TWAS.P.FDR<-p.adjust(pwas_banner$TWAS.P, method = 'fdr')
  
  pwas_banner$`Ensembl ID`<-gsub('\\..*','',pwas_banner$ID)
  pwas_banner$`Gene Symbol`<-gsub('.*\\.','',pwas_banner$ID)
  
  pwas_banner<-pwas_banner[,c('CHR','P0','P1','Ensembl ID','Gene Symbol','TWAS.Z','TWAS.P','TWAS.P.FDR','COLOC.PP3','COLOC.PP4'), with=T]
  
  pwas_banner<-pwas_banner[pwas_banner$TWAS.P.FDR < 0.05 & (pwas_banner$COLOC.PP4-pwas_banner$COLOC.PP3)/pwas_banner$COLOC.PP4 > 0.8,]
  
  gene_assoc<-rbind(gene_assoc, data.frame(ID=pwas_banner$`Gene Symbol`,
                                           DIR=sign(pwas_banner$TWAS.Z)))
  
  
}

if(smr_protein_panel_rosmap_logical){
  
  smr_rosmap<-fread(paste0('results/',opt$gwas,'/smr/rosmap/',opt$gwas,'_smr_rosmap_GW.txt.gz'))
  
  smr_rosmap$p_SMR.FDR<-p.adjust(smr_rosmap$p_SMR, method = 'fdr')
  smr_rosmap<-smr_rosmap[,c('ProbeChr','Probe_bp','ensembl_gene_id','external_gene_name','b_SMR','p_SMR','p_SMR.FDR','p_HEIDI'), with=T]
  
  smr_rosmap<-smr_rosmap[smr_rosmap$p_SMR.FDR < 0.05 & smr_rosmap$p_HEIDI > 0.05,]
  
  gene_assoc<-rbind(gene_assoc, data.frame(ID=smr_rosmap$external_gene_name,
                                           DIR=sign(smr_rosmap$b_SMR)))
  
  
}

gene_assoc<-gene_assoc[!duplicated(gene_assoc),]
gene_assoc<-gene_assoc[!is.na(gene_assoc$ID),]

# The predicted direction varies across some panels
# Restricing to genes with consistent direction may be a little harsh, as proxy tissues may be included
# Provide drugs matching either direction of effect
# Perhaps highlight which panel/method implicated each association

# Send query to DGIdb using API
dir.create(paste0('results/',opt$gwas,'/DGIdb'))
chunks<-split(unique(gene_assoc$ID), ceiling(seq_along(unique(gene_assoc$ID))/100))

gene_assoc_drug_all<-NULL
for(chunk_i in 1:length(chunks)){
  
  system(paste0('curl -k https://dgidb.org/api/v2/interactions.json?genes=',paste(chunks[[chunk_i]], collapse=','),' | python -mjson.tool > results/',opt$gwas,'/DGIdb/query_',chunk_i,'.json'))

  query<-fromJSON(file=paste0('results/',opt$gwas,'/DGIdb/query_',chunk_i,'.json'))
  
  # List agonists/antagonists of genes
  interactions<-NULL
  for(i in 1:length(query$matchedTerms)){
    if(length(query$matchedTerms[[i]]$interactions) > 0){
      for(k in 1:length(query$matchedTerms[[i]]$interactions)){
        if(length(query$matchedTerms[[i]]$interactions[[k]]$interactionTypes) > 0){
          interactions<-rbind(interactions, data.frame(ID=query$matchedTerms[[i]]$searchTerm,
                                                       Drug=query$matchedTerms[[i]]$interactions[[k]]$drugName,
                                                       Type=query$matchedTerms[[i]]$interactions[[k]]$interactionTypes[1],
                                                       Sources=query$matchedTerms[[i]]$interactions[[k]]$sources,
                                                       Score=query$matchedTerms[[i]]$interactions[[k]]$score,
                                                       IntID=query$matchedTerms[[i]]$interactions[[k]]$interactionId))
        }
      }
    }
  }
  
  if(!is.null(interactions)){
    # Filter results to identify drugs that reverse the gene's association with the GWAS phenotype
    interactions$Direction[(interactions$Type %in% c('activator','agonist','chaperone','cofactor','inducer','partial agonist','positive modulator','stimulator','vaccine'))]<-'Activating'
    interactions$Direction[(interactions$Type %in% c('antagonist','antibody','antisense oligonucleotide','blocker','cleavage','inhibitor','inhibitory allosteric modulator','inverse agonist','negative modulator','partial antagonist','suppressor'))]<-'Inhibitory'
    
    gene_assoc_drug<-merge(gene_assoc, interactions, by='ID')
    gene_assoc_drug<-gene_assoc_drug[(gene_assoc_drug$DIR < 0 & gene_assoc_drug$Direction == 'Activating') | (gene_assoc_drug$DIR > 0 & gene_assoc_drug$Direction == 'Inhibitory'),]
    
    gene_assoc_drug$Sources<-NULL
    gene_assoc_drug$Type<-NULL
    gene_assoc_drug<-gene_assoc_drug[!duplicated(gene_assoc_drug),]
    gene_assoc_drug$GenePhenoAssoc[gene_assoc_drug$DIR > 0]<-'Upregulated'
    gene_assoc_drug$GenePhenoAssoc[gene_assoc_drug$DIR < 0]<-'Downregulated'
    gene_assoc_drug$DrugGeneInt<-gene_assoc_drug$Direction
    
    gene_assoc_drug<-gene_assoc_drug[,c('ID','GenePhenoAssoc','Drug','DrugGeneInt','Score','IntID')]
  
    gene_assoc_drug_all<-rbind(gene_assoc_drug_all, gene_assoc_drug)
  }
}

write.csv(gene_assoc_drug_all, paste0('results/',opt$gwas,'/DGIdb/DGIdb_opposing_clean.csv'), row.names=F, quote=F)

