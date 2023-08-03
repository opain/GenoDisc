#!/usr/bin/Rscript

process_cleaner_log<-function(config, gwas){

  # Identify outdir
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  # Create empty list to store output
  dat<-list()

  # Read in log file
  dat$log<-readLines(paste0(outdir,'/data/gwas_sumstat/',gwas,'/',gwas,'.cleaned.log'))

  # Create val list to store specific values of interest from the log file
  dat$val<-list()

  # Identify number of variants in original GWAS file
  dat$val$n_var_orig<-as.numeric(gsub(' variants.','',gsub('GWAS contains ','',dat$log[grepl('GWAS contains ', dat$log)])))

  # Identify the genome build identified in the GWAS
  dat$val$build<-extract_build(dat$log)

  # Identify the number of SNPs in the GWAS after QC
  dat$val$n_snp_final<-as.numeric(gsub(' variants remain.','',gsub('After removal of SNPs with SE == 0, ','',dat$log[grepl('After removal of SNPs with SE == 0, ', dat$log)])))

  return(dat)
}

extract_build<-function(x){
  build_match<-data.table(t(matrix(unlist(strsplit(x[grepl('match: ', x)], ' match: ')), nrow=2)))
  build_match$V2<-as.numeric(gsub('%','',build_match$V2))/100

  best_match<-list()
  best_match$build<-build_match$V1[build_match$V2 == max(build_match$V2)]
  best_match$overlap<-build_match$V2[build_match$V2 == max(build_match$V2)]

  return(best_match)
}

process_focus_log<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-list()

  dat$log<-readLines(paste0(outdir,'/data/gwas_sumstat/',gwas,'/',gwas,'.cleaned.munged.log'))

  dat$val<-list()
  dat$val$lambda_gc<-as.numeric(gsub('.*Lambda GC = ','',dat$log[grepl('Lambda GC', dat$log)]))
  dat$val$max_chi2<-as.numeric(gsub('.*Max chi\\^2 = ','',dat$log[grepl('Max chi\\^2', dat$log)]))
  dat$val$n_sig_snp<-as.numeric(gsub('.* ','',gsub(' Genome-wide significant SNPs.*','',dat$log[grepl('Genome-wide significant SNPs', dat$log)])))

  return(dat)
}

process_ldsc_log<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(config[grepl('ldsc:',config)] == "ldsc: T"){

    dat<-list()

    dat$log<-readLines(paste0(outdir,'/results/',gwas,'/ldsc/',gwas,'_ldsc_res.log'))

    dat$val<-list()
    tmp<-gsub('Total Observed scale h2: ','',dat$log[grepl('Total Observed scale h2:', dat$log)])
    dat$val$obs_h2_est<-as.numeric(gsub(' .*','',tmp))
    dat$val$obs_h2_se<-as.numeric(gsub('\\)','',gsub('.*\\(','',tmp)))

    tmp<-gsub('Intercept: ','',dat$log[grepl('Intercept: ', dat$log)])
    dat$val$int_est<-as.numeric(gsub(' .*','',tmp))
    dat$val$int_se<-as.numeric(gsub('\\)','',gsub('.*\\(','',tmp)))

  }

  return(dat)
}

process_susie<-function(outdir, gwas, L){
  if(L==1){
  dat<-fread(paste0(outdir,'/results/',gwas,'/finemap/',gwas,'.GW.finemap.L1.csv'))
  }
  if(L==10){
  dat<-fread(paste0(outdir,'/results/',gwas,'/finemap/',gwas,'.GW.finemap.csv'))
  }

  dat$variable<-NULL
  dat$Gene[is.na(dat$Gene)]<-'None'
  dat<-dat[order(dat$CHR, dat$BP),]

  return(dat)
}

tidy_panel_names<-function(x){
  
  panel_names<-data.table(t(matrix(c("CMC.BRAIN.RNASEQ_SPLICING","Brain - DLPFC (CMC)","Brain_Anterior_cingulate_cortex_BA24","Brain - Anterior cingulate cortex BA24 (GTEx)","psychencode","Brain - DLPFC (PsychENCODE)","Brain_Cerebellar_Hemisphere","Brain - Cerebellar hemisphere (GTEx)","NTR.BLOOD.RNAARR","Blood (NTR)","Brain_Cortex","Brain - Cortex (GTEx)","Brain_Caudate_basal_ganglia","Brain - Caudate basal ganglia (GTEx)","YFS.BLOOD.RNAARR","Blood (YFS)","Brain_Hippocampus","Brain - Hippocampus (GTEx)","Brain_Substantia_nigra","Brain - Substantia nigra (GTEx)","Brain_Nucleus_accumbens_basal_ganglia","Brain - Nucleus accumbens basal ganglia (GTEx)","Brain_Putamen_basal_ganglia","Brain - Putamen basal ganglia (GTEx)","Brain_Frontal_Cortex_BA9","Brain - Frontal cortex BA9 (GTEx)","Whole_Blood","Blood (GTEx)","CMC.BRAIN.RNASEQ","Brain - DLPFC (CMC)","Brain_Cerebellum","Brain - Cerebellum (GTEx)","Brain_Spinal_cord_cervical_c-1","Brain - Spinal cord vervical c-1 (GTEx)","Brain_Hypothalamus","Brain - Hypothalamus (GTEx)","kcl_brainbank_motor_cortex","Brain - Motor cortex (KCL Brain Bank)","Brain_Amygdala","Brain - Amygdala (GTEx)"), nrow=2)))

  names(panel_names)<-c('original','clean')
  x_tab<-data.table(original=x)

  x_tab<-merge(x_tab, panel_names, by='original')
  x_tab<-x_tab[match(x, x_tab$original),]

  return(x_tab$clean)
  
}

read_fusion_exp<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  # Check whether TWAS was performed
  twas_panel_psychencode_logical<-config[grepl('twas_panel_psychencode:',config)] == "twas_panel_psychencode: T"
  twas_panel_fusion_logical<-config[grepl('twas_panel_fusion:',config)] == "twas_panel_fusion: T"
  twas_logical<-any(twas_panel_psychencode_logical, twas_panel_fusion_logical)

  dat<-NULL

  if(twas_logical){

    dat<-list()

    # Identify TWAS panels included
    gtex_weights<-config[grepl('^gtex_weights', config)]
    gtex_weights<-unlist(strsplit(gsub('"','',gsub('\\]','',gsub('.*\\[','',gtex_weights))),','))

    non_gtex_weights<-config[grepl('^non_gtex_weights', config)]
    non_gtex_weights<-unlist(strsplit(gsub('"','',gsub('\\]','',gsub('.*\\[','',non_gtex_weights))),','))

    twas_weights<-c(gtex_weights, non_gtex_weights)

    psychencode_weights_log<-config[grepl('^twas_panel_psychencode:', config)]
    if(psychencode_weights_log == "twas_panel_psychencode: T"){
        twas_weights<-c(twas_weights, 'psychencode')
    }

    external_weights_log<-config[grepl('^external_weights:', config)]
    if(external_weights_log == "external_weights: T"){
        external_weights<-config[grepl('^external_weights_pos_path', config)]
        external_weights<-gsub('.pos','',basename(unlist(strsplit(gsub('"','',gsub('\\]','',gsub('.*\\[','',external_weights))),','))))
        twas_weights<-c(twas_weights, external_weights)
    }

    # Create table of TWAS panels included
    dat$panels<-NULL
    for(i in twas_weights){
        pos<-fread(paste0('resources/data/fusion_snp_weights/',i,'/',i,'.pos'))
        dat$panels<-rbind(dat$panels, data.frame( Type='Expression',
                                                  Software='FUSION',
                                                  Panel=i,
                                                  N_indiv=median(pos$N),
                                                  N_gene=nrow(pos)))
    }

    # Update panel names to be easier to read
    dat$panels$Panel<-tidy_panel_names(dat$panels$Panel)

    # Read in TWAS results
    dat$res<-fread(paste0(outdir,'/results/',gwas,'/twas/',gwas,'_twas_GW_clean.txt.gz'))
    dat$res$FILE<-NULL
    dat$res$MODEL<-NULL
    dat$res<-dat$res[!is.na(dat$res$TWAS.Z),]

    # Calculate P accounting for multiple testing
    dat$res$TWAS.P.FDR<-p.adjust(dat$res$TWAS.P, method='fdr')
    dat$res<-dat$res[,c('CHR','P0','P1','ensembl_gene_id','external_gene_name','PANEL','TWAS.Z','TWAS.P','TWAS.P.FDR','COLOC.PP3','COLOC.PP4'), with=T]
    names(dat$res)[names(dat$res) == 'ensembl_gene_id']<-'Ensembl ID'
    names(dat$res)[names(dat$res) == 'external_gene_name']<-'Gene Symbol'
    dat$res$COLOC_logical<-(dat$res$COLOC.PP4-dat$res$COLOC.PP3)/dat$res$COLOC.PP4 > 0.8

    dat$res<-dat$res[order(dat$res$CHR, dat$res$P0),]

    dat$res$PANEL<-tidy_panel_names(dat$res$PANEL)

  }

  return(dat)
}

# Read in PWAS results
read_pwas<-function(outdir, gwas, panel){
  pwas_all<-NULL
  for(panel_i in panel){

    pwas_files<-list.files(path=paste0(outdir,'/results/',gwas,'/pwas/',panel,'/'), pattern=paste0(gwas,'_pwas_',panel,'_chr'))

    pwas<-NULL
    for(i in pwas_files){
      pwas<-rbind(pwas, fread(paste0(outdir,'/results/',gwas,'/pwas/rosmap/',i)))
    }

    pwas$PANEL<-NULL

    if(panel == 'rosmap'){
      pwas$PANEL<-"Brain - DLPFC (ROSMAP)"
    }
    if(panel == 'banner'){
      pwas$PANEL<-"Brain - DLPFC (Banner)"
    }

    pwas_all<-rbind(pwas, pwas)
  }

  pwas_all<-pwas_all[!is.na(pwas_all$TWAS.P),]
  pwas_all$TWAS.P.FDR<-p.adjust(pwas_all$TWAS.P, method = 'fdr')
  pwas_all$`Ensembl ID`<-gsub('\\..*','',pwas_all$ID)
  pwas_all$`Gene Symbol`<-gsub('.*\\.','',pwas_all$ID)
  pwas_all<-pwas_all[,c('PANEL','CHR','P0','P1','Ensembl ID','Gene Symbol','TWAS.Z','TWAS.P','TWAS.P.FDR','COLOC.PP3','COLOC.PP4'), with=T]
  pwas_all$COLOC_logical<-(pwas_all$COLOC.PP4-pwas_all$COLOC.PP3)/pwas_all$COLOC.PP4 > 0.8
  pwas_all<-pwas_all[order(pwas_all$CHR, pwas_all$P0),]
  names(pwas_all)<-gsub('TWAS','pwas_all',names(pwas_all))

  return(pwas_all)
}

read_fusion_protein<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  pwas_panel_rosmap_logical<-config[grepl('pwas_panel_rosmap:',config)] == "pwas_panel_rosmap: T"
  pwas_panel_banner_logical<-config[grepl('pwas_panel_banner:',config)] == "pwas_panel_banner: T"

  dat<-NULL

  if(any(pwas_panel_rosmap_logical, pwas_panel_banner_logical)){

    dat<-list()

    # List panel information
    if(pwas_panel_rosmap_logical){
      dat$panels<-rbind(dat$panels, data.frame(Type='Protein',
                                                                    Software='FUSION',
                                                                    Panel="Brain - DLPFC (ROSMAP)",
                                                                    N_indiv=376,
                                                                    N_gene=1477))
    }

    if(pwas_panel_banner_logical){
      dat$panels<-rbind(dat$panels, data.frame(Type='Protein',
                                                                    Software='FUSION',
                                                                    Panel="Brain - DLPFC (Banner)",
                                                                    N_indiv=152,
                                                                    N_gene=1148))
    }

    # Read in the results from each PANEL
    panels<-NULL
    if(pwas_panel_rosmap_logical){
      panels<-c(panels,'rosmap')
    }
    if(pwas_panel_banner_logical){
      panels<-c(panels,'banner')
    }

    pwas<-read_pwas(outdir=outdir, gwas=gwas, panel=panels)

    dat$results<-pwas
  }

  return(dat)
}

read_smr_exp<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  # Check whether SMR with expression data was performed
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

  smr_expression_logical<-any(smr_expression_panel_psychencode_logical,
                          metabrain_logical,
                          smr_expression_panel_eqtlgen_logical)

  panels<-NULL
  res<-NULL

  if(smr_expression_logical){
    # Create a table listing SMR panels
    panels<-NULL

    if(smr_expression_panel_psychencode_logical){
      panels<-rbind(panels, data.frame(Type='Expression',
                                                                    Software='SMR',
                                                                    Panel="Brain - DLPFC (PsychENCODE)",
                                                                    N_indiv=1321,
                                                                    N_gene=24560))
    }

    if(smr_expression_panel_metabrain_basalganglia_logical){
      panels<-rbind(panels, data.frame(Type='Expression',
                                                                    Software='SMR',
                                                                    Panel="Brain - Basalganglia (MetaBrain)",
                                                                    N_indiv=574,
                                                                    N_gene=18406))
    }

    if(smr_expression_panel_metabrain_cerebellum_logical){
      panels<-rbind(panels, data.frame(Type='Expression',
                                                                    Software='SMR',
                                                                    Panel="Brain - Cerebellum (MetaBrain)",
                                                                    N_indiv=723,
                                                                    N_gene=18417))
    }

    if(smr_expression_panel_metabrain_cortex_logical){
      panels<-rbind(panels, data.frame(Type='Expression',
                                                                    Software='SMR',
                                                                    Panel="Brain - Cortex (MetaBrain)",
                                                                    N_indiv=6601,
                                                                    N_gene=18414))
    }

    if(smr_expression_panel_metabrain_hippocampus_logical){
      panels<-rbind(panels, data.frame(Type='Expression',
                                                                    Software='SMR',
                                                                    Panel="Brain - Hippocampus (MetaBrain)",
                                                                    N_indiv=206,
                                                                    N_gene=18406))
    }
                          
    if(smr_expression_panel_metabrain_spinalcord_logical){
      panels<-rbind(panels, data.frame(Type='Expression',
                                                                    Software='SMR',
                                                                    Panel="Brain - Spinalcord (MetaBrain)",
                                                                    N_indiv=285,
                                                                    N_gene=18417))
    }

    if(smr_expression_panel_eqtlgen_logical){
      panels<-rbind(panels, data.frame(Type='Expression',
                                                                    Software='SMR',
                                                                    Panel="Blood (eQTLGen)",
                                                                    N_indiv=31684,
                                                                    N_gene=19250))
    }

    # Read in the results
    res<-NULL

    if(smr_expression_panel_psychencode_logical){
      
      smr_psychencode_files<-list.files(path=paste0(outdir,'/results/',gwas,'/smr/psychencode/'), pattern=paste0(gwas,'_smr_psychencode_chr'))
      smr_psychencode_files<-smr_psychencode_files[grepl('.smr$', smr_psychencode_files)]
      
      smr_psychencode<-NULL
      for(i in smr_psychencode_files){
        smr_psychencode<-rbind(smr_psychencode, fread(paste0(outdir,'/results/',gwas,'/smr/psychencode/',i)))
      }
      
      smr_psychencode_tmp<-smr_psychencode
      
      smr_psychencode_tmp$PANEL<-"Brain - DLPFC (PsychENCODE)"
      
      smr_psychencode_tmp<-smr_psychencode_tmp[,c('PANEL','ProbeChr','Probe_bp','probeID','Gene','b_SMR','se_SMR','p_SMR','p_HEIDI'), with=T]
      
      names(smr_psychencode_tmp)<-c('PANEL','CHR','BP','Ensembl ID','Gene Symbol','b_SMR','se_SMR','p_SMR','p_HEIDI')
      
      res<-rbind(res, smr_psychencode_tmp)
    }

    if(smr_expression_panel_eqtlgen_logical){
      
      smr_eqtlgen<-fread(paste0(outdir,'/results/',gwas,'/smr/eqtlgen/',gwas,'_smr_eqtlgen_GW.txt.gz'))
      
      smr_eqtlgen_tmp<-smr_eqtlgen
      
      smr_eqtlgen_tmp$PANEL<-'Blood (eQTLGen)'
      
      smr_eqtlgen_tmp<-smr_eqtlgen_tmp[,c('PANEL','ProbeChr','Probe_bp','probeID','external_gene_name','b_SMR','se_SMR','p_SMR','p_HEIDI'), with=T]
      
      names(smr_eqtlgen_tmp)<-c('PANEL','CHR','BP','Ensembl ID','Gene Symbol','b_SMR','se_SMR','p_SMR','p_HEIDI')

      res<-rbind(res, smr_eqtlgen_tmp)

    }

    if(metabrain_logical){

      smr_metabrain_all<-fread(paste0(outdir,'/results/',gwas,'/smr/metabrain/',gwas,'_smr_metabrain_GW.txt.gz'))
      
      smr_metabrain_tmp<-smr_metabrain_all
      
      smr_metabrain_tmp$PANEL<-paste0('Brain - ', smr_metabrain_tmp$PANEL, " (MetaBrain)")
      
      smr_metabrain_tmp$probeID<-gsub('\\..*','',smr_metabrain_tmp$probeID)

      smr_metabrain_tmp<-smr_metabrain_tmp[,c('PANEL','ProbeChr','Probe_bp','probeID','external_gene_name','b_SMR','se_SMR','p_SMR','p_HEIDI'), with=T]
      
      names(smr_metabrain_tmp)<-c('PANEL','CHR','BP','Ensembl ID','Gene Symbol','b_SMR','se_SMR','p_SMR','p_HEIDI')

      res<-rbind(res, smr_metabrain_tmp)

    }

    # Adjust for multiple testing
    res$p_SMR.FDR<-p.adjust(res$p_SMR, method='fdr')

    res<-res[,c('PANEL','CHR','BP','Ensembl ID','Gene Symbol','b_SMR','se_SMR','p_SMR','p_SMR.FDR','p_HEIDI'), with=T]

    res<-res[order(res$CHR, res$BP),]
  
  }

  return(list(panels=panels,
              results=res))

}


read_smr_protein<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(config[grepl('smr_protein_panel_rosmap:',config)] == "smr_protein_panel_rosmap: T"){
    # Insert panel information
    dat$panels<-data.frame( Type='Protein',
                            Software='SMR',
                            Panel='Brain - DLPFC (ROSMAP)',
                            N_indiv=376,
                            N_gene=7809)
  
    # Insert results
    dat$results<-fread(paste0(outdir,'/results/',gwas,'/smr/rosmap/',gwas,'_smr_rosmap_GW.txt.gz'))

    dat$results$p_SMR.FDR<-p.adjust(dat$results$p_SMR, method = 'fdr')
    dat$results<-dat$results[,c('ProbeChr','Probe_bp','ensembl_gene_id','external_gene_name','b_SMR','se_SMR','p_SMR','p_SMR.FDR','p_HEIDI'), with=T]

    names(dat$results)<-c('CHR','BP','Ensembl ID','Gene Symbol','b_SMR','se_SMR','p_SMR','p_SMR.FDR','p_HEIDI')

    dat$results$PANEL<-"Brain - DLPFC (ROSMAP)"

    dat$results<-dat$results[order(dat$results$CHR, dat$results$BP),]

  }

  return(dat)
}

read_magma_gene<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(config[grepl('magma_gene:',config)] == "magma_gene: T"){
    dat<-fread(paste0(outdir,'/results/',gwas,'/magma/magma_gene_level.clean.csv'))
    dat$P.FDR<-p.adjust(dat$P, method = 'fdr')
    dat<-dat[order(dat$CHR, dat$START),]
  }
  
  return(dat)
}

identify_nearest<-function(x){
  if(!is.null(x)){
    # Identify genes containing lead SNP (Some are within multiple genes)
    nearest_list<-unlist(strsplit(x,','))
    nearest_list<-nearest_list[!grepl("kb", nearest_list)]
    nearest_list<-gsub(' ','',nearest_list)
    
    # Identify nearest genes
    nearest_list_2<-gsub(' .*','',gsub(',.*','',x))
  
    nearest<-unique(c(nearest_list,nearest_list_2))
    nearest<-nearest[nearest != 'None']
  } else{
    nearest<-NULL
  }

  return(nearest)
}


read_twas_gsea_drug<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(twas_gsea_drugtargetor_logical<-config[grepl('twas_gsea_drugtargetor:',config)] == "twas_gsea_drugtargetor: T"){

    atc<-fread('resources/data/atc/atc_20220201.txt', sep='!')
    names(atc)<-c('Code','Name')
    atc$Name<-tolower(atc$Name)

    atc_labels<-atc[nchar(atc$Code) == 4,]

    weights<-read.table(paste0(outdir,'/results/',gwas,'/twas/list_of_weights.txt'))$V1
    weights<-weights[!grepl('SPLIC',weights)]

    dat<-NULL
    for(i in weights){
      res<-fread(paste0(outdir,'/results/',gwas,'/twas/drugtargetor/twas_gsea_drugtargetor_',i,'.competitive.clean.csv'))
      res$Panel<-i
      dat<-rbind(dat, res)
    }

    dat$P.FDR<-p.adjust(dat$P, method='fdr')

    dat$NAME<-paste(toupper(substr(dat$NAME, 1, 1)), substr(dat$NAME, 2, nchar(dat$NAME)), sep="")
    dat$NAME<-gsub('\\.',' ',dat$NAME)

    # Format ATC codes and insert ATC descriptions
    tmp<-lapply(strsplit(dat$ATC, ','), function(x) substr(x, 1, 4))
    tmp2<-lapply(tmp, insert_atc_desc, atc_labels)
    dat$ATC_code<-unlist(lapply(tmp, function(x) paste0(x, collapse=';')))
    dat$ATC_desc<-unlist(lapply(tmp2, function(x) paste0(x, collapse=';')))

    dat<-dat[order(dat$P),]

    dat$Panel<-tidy_panel_names(dat$Panel)

    dat<-dat[,c("NAME","Panel","Estimate","SE","P","P.FDR","ATC_code","ATC_desc") , with=F]
    names(dat)<-c('Drug','Panel','Estimate','SE','P','P.FDR','ATC Code','ATC Description')

    dat$ChEMBL<-paste0('<a href="https://www.ebi.ac.uk/chembl/g/#search_results/all/query=',dat$Name,'">','Link','</a>')

  }
  return(dat)
}

read_twas_gsea_atc<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(twas_gsea_drugtargetor_logical<-config[grepl('twas_gsea_drugtargetor:',config)] == "twas_gsea_drugtargetor: T"){

    weights<-read.table(paste0(outdir,'/results/',gwas,'/twas/list_of_weights.txt'))$V1
    weights<-weights[!grepl('SPLIC',weights)]

    dat<-NULL
    for(i in weights){
      res<-fread(paste0(outdir,'/results/',gwas,'/twas/drugtargetor/twas_gsea_',i,'_res_atc_res.csv'))
      res$Panel<-i
      dat<-rbind(dat, res)
    }

    dat$P.FDR_all<-p.adjust(dat$P, method = 'fdr')

    dat<-dat[,c('Panel','ATC','Name','N','Estimate','Class_Median','Non_Class_Median','P','P.FDR_all'),with=T]
    names(dat)<-c('Panel','ATC Code','ATC Description','N Drugs','Estimate','Class Median T','Non-class Median T','P','P.FDR')

    dat<-dat[order(dat$P.onesided),]

    # Update labels
    dat$Panel<-tidy_panel_names(dat$Panel)
    
  }

  return(dat)
}

read_magma_drug<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(config[grepl('magma_drugtargetor:',config)] == "magma_drugtargetor: T"){

    atc<-fread('resources/data/atc/atc_20220201.txt', sep='!')
    names(atc)<-c('Code','Name')
    atc$Name<-tolower(atc$Name)

    atc_labels<-atc[nchar(atc$Code) == 4,]

    dat<-fread(paste0(outdir,'/results/',gwas,'/magma/magma_drug_targetor.clean.csv'))
    dat<-dat[order(dat$P),]
    dat$P.FDR<-p.adjust(dat$P, method='fdr')

    dat$NAME<-paste(toupper(substr(dat$NAME, 1, 1)), substr(dat$NAME, 2, nchar(dat$NAME)), sep="")
    
    # Format ATC codes and insert ATC descriptions
    tmp<-lapply(strsplit(dat$ATC, ','), function(x) substr(x, 1, 4))
    tmp2<-lapply(tmp, insert_atc_desc, atc_labels)
    dat$ATC_code<-unlist(lapply(tmp, function(x) paste0(x, collapse=';')))
    dat$ATC_desc<-unlist(lapply(tmp2, function(x) paste0(x, collapse=';')))

    dat<-dat[,c("NAME","NGENES","BETA","SE","P","P.FDR",'ATC_code','ATC_desc') , with=F]
    names(dat)<-c('Name','N Genes',"BETA","SE",'P','P.FDR','ATC Code','ATC Description')

    dat$ChEMBL<-paste0('<a href="https://www.ebi.ac.uk/chembl/g/#search_results/all/query=',dat$Name,'">','Link','</a>')

  }

  return(dat)
}

insert_atc_desc <- function(x, replacement_df) {
  idx <- match(x, replacement_df$Code)
  x <- replacement_df$Name[idx]
  return(x)
}

read_gcsc<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(config[grepl('gcsc:',config)] == "gcsc: T"){

    atc<-fread('resources/data/atc/atc_20220201.txt', sep='!')
    names(atc)<-c('Code','Name')
    atc$Name<-tolower(atc$Name)

    atc_labels<-atc[nchar(atc$Code) == 4,]

    dat<-fread(paste0('/users/k1806347/oliverpainfel/Analyses/GCSC/GCSC_Brain_Blood_ALS_only_results.csv'))
    dat$Drug<-paste(toupper(substr(dat$Drug, 1, 1)), substr(dat$Drug, 2, nchar(dat$Drug)), sep="")

    dat<-dat[order(dat$P),]
    dat$P.FDR<-p.adjust(dat$P, method='fdr')

    # Format ATC codes and insert ATC descriptions
    tmp<-lapply(strsplit(dat$ATC, '_'), function(x) substr(x, 1, 4))
    tmp2<-lapply(tmp, insert_atc_desc, atc_labels)
    dat$ATC_code<-unlist(lapply(tmp, function(x) paste0(x, collapse=';')))
    dat$ATC_desc<-unlist(lapply(tmp2, function(x) paste0(x, collapse=';')))

    dat$Z<--qnorm(dat$P/2)
    dat$Z<-dat$Z*sign(dat$Enrichment)

    dat<-dat[,c("Drug","Enrichment","SE","Z","P","P.FDR",'ATC_code','ATC_desc') , with=F]
    names(dat)<-c('Name',"Enrichment","SE","Z","P","P.FDR",'ATC Code','ATC Description')

    dat$ChEMBL<-paste0('<a href="https://www.ebi.ac.uk/chembl/g/#search_results/all/query=',dat$Name,'">','Link','</a>')
  }
  return(dat)
}

read_magma_drug_atc<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(config[grepl('magma_drugtargetor:',config)] == "magma_drugtargetor: T"){

    dat<-fread(paste0(outdir,'/results/',gwas,'/magma/magma_drug_targetor_atc_res.csv'))
    dat$P.FDR<-p.adjust(dat$P, method = 'fdr')

    dat<-dat[,c("ATC","N","P","P.FDR","Name"), with=F]
    names(dat)<-c("ATC Code","N Drugs","P","P.FDR","ATC Description")

    dat<-dat[order(dat$P),]

  }
  return(dat)
}

read_gcsc_atc<-function(config, gwas){
  
  outdir<-gsub('outdir: ','', config[grepl('outdir: ',config)])

  dat<-NULL

  if(config[grepl('gcsc:',config)] == "gcsc: T"){

    dat<-fread(paste0(outdir,'/results/',gwas,'/gcsc/',gwas,'_drugtargetor_gcsc_res_atc.csv'))
    
    dat$P.FDR<-p.adjust(dat$P, method = 'fdr')

    dat<-dat[,c("ATC","N","P","P.FDR","Name"), with=F]
    names(dat)<-c("ATC Code","N Drugs","P","P.FDR","ATC Description")

    dat<-dat[order(dat$P),]

  }
  return(dat)
}
