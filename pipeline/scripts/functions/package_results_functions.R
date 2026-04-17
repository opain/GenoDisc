#!/usr/bin/Rscript

process_cleaner_log<-function(config, gwas){

  # Identify outdir
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  # Create empty list to store output
  dat<-list()

  # Read in log file
  dat$log<-readLines(paste0(outdir,'/results/',gwas,'/gwas_sumstat/',gwas,'.cleaned.log'))

  # Create val list to store specific values of interest from the log file
  dat$val<-list()

  # Identify number of variants in original GWAS file
  dat$val$n_var_orig<-as.numeric(gsub(' variants.','',gsub('GWAS contains ','',dat$log[grepl('GWAS contains \\d+ variants\\.', dat$log)])))

  # Identify the genome build identified in the GWAS
  dat$val$build<-extract_build(dat$log)

  # Identify the number of SNPs in the GWAS after QC
  dat$val$n_snp_final<-as.numeric(gsub(' variants remain.','',gsub('After removal of SNPs with SE == 0, ','',dat$log[grepl('After removal of SNPs with SE == 0, ', dat$log)])))

  return(dat)
}

extract_build<-function(x){
  match_lines<-x[grepl('match: ', x)]
  build_names<-gsub(' match:.*','',match_lines)
  target_pct<-as.numeric(gsub('%.*','',gsub('.* match: ','',match_lines)))/100

  best_match<-list()
  best_match$build<-build_names[which.max(target_pct)]
  best_match$overlap<-max(target_pct)

  return(best_match)
}

process_focus_log<-function(config, gwas){

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-list()

  dat$log<-readLines(paste0(outdir,'/results/',gwas,'/gwas_sumstat/',gwas,'.cleaned.munged.log'))

  dat$val<-list()
  dat$val$lambda_gc<-as.numeric(gsub('.*Lambda GC = ','',dat$log[grepl('Lambda GC', dat$log)]))
  dat$val$max_chi2<-as.numeric(gsub('.*Max chi\\^2 = ','',dat$log[grepl('Max chi\\^2', dat$log)]))
  dat$val$n_sig_snp<-as.numeric(gsub('.* ','',gsub(' Genome-wide significant SNPs.*','',dat$log[grepl('Genome-wide significant SNPs', dat$log)])))

  return(dat)
}

process_ldsc_log<-function(config, gwas){

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'ldsc', return_obj = F) == "T"){

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

  panel_names<-data.table(t(matrix(c("CMC.BRAIN.RNASEQ_SPLICING","Brain - DLPFC - Splice (CMC)","Brain_Anterior_cingulate_cortex_BA24","Brain - Anterior cingulate cortex BA24 (GTEx)","psychencode","Brain - DLPFC (PsychENCODE)","Brain_Cerebellar_Hemisphere","Brain - Cerebellar hemisphere (GTEx)","NTR.BLOOD.RNAARR","Blood (NTR)","Brain_Cortex","Brain - Cortex (GTEx)","Brain_Caudate_basal_ganglia","Brain - Caudate basal ganglia (GTEx)","YFS.BLOOD.RNAARR","Blood (YFS)","Brain_Hippocampus","Brain - Hippocampus (GTEx)","Brain_Substantia_nigra","Brain - Substantia nigra (GTEx)","Brain_Nucleus_accumbens_basal_ganglia","Brain - Nucleus accumbens basal ganglia (GTEx)","Brain_Putamen_basal_ganglia","Brain - Putamen basal ganglia (GTEx)","Brain_Frontal_Cortex_BA9","Brain - Frontal cortex BA9 (GTEx)","Whole_Blood","Blood (GTEx)","CMC.BRAIN.RNASEQ","Brain - DLPFC (CMC)","Brain_Cerebellum","Brain - Cerebellum (GTEx)","Brain_Spinal_cord_cervical_c-1","Brain - Spinal cord vervical c-1 (GTEx)","Brain_Hypothalamus","Brain - Hypothalamus (GTEx)","kcl_brainbank_motor_cortex","Brain - Motor cortex (KCL Brain Bank)","Brain_Amygdala","Brain - Amygdala (GTEx)"), nrow=2)))

  names(panel_names)<-c('original','clean')
  x_tab<-data.table(original=x)

  x_tab<-merge(x_tab, panel_names, by='original')
  x_tab<-x_tab[match(x, x_tab$original),]

  return(x_tab$clean)

}

read_fusion_exp<-function(config, gwas){

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)
  resdir <- read_param(config = config, param = 'resdir', return_obj = F)

  # Check whether TWAS was performed
  twas_panel_psychencode_logical<-read_param(config = config, param = 'twas_panel_psychencode', return_obj = F) == "T"
  twas_panel_fusion_logical<-read_param(config = config, param = 'twas_panel_fusion', return_obj = F) == "T"
  twas_logical<-any(twas_panel_psychencode_logical, twas_panel_fusion_logical)

  dat<-NULL

  if(twas_logical){

    dat<-list()

    # Identify TWAS panels included
    gtex_weights <- read_param(config = config, param = 'gtex_weights', return_obj = F)
    non_gtex_weights <- read_param(config = config, param = 'non_gtex_weights', return_obj = F)

    twas_weights<-c(gtex_weights, non_gtex_weights)

    if(twas_panel_psychencode_logical){
        twas_weights<-c(twas_weights, 'psychencode')
    }

    external_weights_flag <- read_param(config = config, param = 'external_weights', return_obj = F)
    if(external_weights_flag == "T"){
        external_weights_pos_path <- read_param(config = config, param = 'external_weights_pos_path', return_obj = F)
        external_weights <- gsub('.pos','',basename(external_weights_pos_path))
        twas_weights<-c(twas_weights, external_weights)
    }

    # Create table of TWAS panels included
    dat$panels<-NULL
    for(i in twas_weights){
        pos<-fread(paste0(resdir, '/data/fusion_snp_weights/',i,'/',i,'.pos'))
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

    pwas_files<-list.files(path=paste0(outdir,'/results/',gwas,'/pwas/',panel_i,'/'), pattern=paste0(gwas,'_pwas_',panel_i,'_chr'))

    pwas<-NULL
    for(i in pwas_files){
      pwas<-rbind(pwas, fread(paste0(outdir,'/results/',gwas,'/pwas/',panel_i,'/',i)))
    }

    pwas$PANEL<-NULL

    if(panel_i == 'rosmap'){
      pwas$PANEL<-"Brain - DLPFC (ROSMAP)"
    }
    if(panel_i == 'banner'){
      pwas$PANEL<-"Brain - DLPFC (Banner)"
    }

    pwas_all<-rbind(pwas_all, pwas)
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

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  pwas_panel_rosmap_logical<-read_param(config = config, param = 'pwas_panel_rosmap', return_obj = F) == "T"
  pwas_panel_banner_logical<-read_param(config = config, param = 'pwas_panel_banner', return_obj = F) == "T"

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

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  # Check whether SMR with expression data was performed
  smr_expression_panel_psychencode_logical<-read_param(config = config, param = 'smr_expression_panel_psychencode', return_obj = F) == "T"

  smr_expression_panel_metabrain_basalganglia_logical<-read_param(config = config, param = 'smr_expression_panel_metabrain_basalganglia', return_obj = F) == "T"
  smr_expression_panel_metabrain_cerebellum_logical<-read_param(config = config, param = 'smr_expression_panel_metabrain_cerebellum', return_obj = F) == "T"
  smr_expression_panel_metabrain_cortex_logical<-read_param(config = config, param = 'smr_expression_panel_metabrain_cortex', return_obj = F) == "T"
  smr_expression_panel_metabrain_hippocampus_logical<-read_param(config = config, param = 'smr_expression_panel_metabrain_hippocampus', return_obj = F) == "T"
  smr_expression_panel_metabrain_spinalcord_logical<-read_param(config = config, param = 'smr_expression_panel_metabrain_spinalcord', return_obj = F) == "T"

  smr_expression_panel_eqtlgen_logical<-read_param(config = config, param = 'smr_expression_panel_eqtlgen', return_obj = F) == "T"

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

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'smr_protein_panel_rosmap', return_obj = F) == "T"){
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

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'magma_gene', return_obj = F) == "T"){
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


read_twas_gsea_drug<-function(config, gwas, mode = 'directional'){

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)
  resdir <- read_param(config = config, param = 'resdir', return_obj = F)

  dat<-NULL

  flag_param <- if(mode == 'nondirectional') 'twas_gsea_drugtargetor_nondirectional' else 'twas_gsea_drugtargetor'
  suffix <- if(mode == 'nondirectional') '_nondir' else ''

  if(read_param(config = config, param = flag_param, return_obj = F) == "T"){

    atc<-fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep='!')
    names(atc)<-c('Code','Name')
    atc$Name<-tolower(atc$Name)

    atc_labels<-atc[nchar(atc$Code) == 4,]

    weights<-read.table(paste0(outdir,'/results/',gwas,'/twas/list_of_weights.txt'))$V1
    weights<-weights[!grepl('SPLIC',weights)]

    dat<-NULL
    for(i in weights){
      res<-fread(paste0(outdir,'/results/',gwas,'/twas/drugtargetor/twas_gsea_drugtargetor',suffix,'_',i,'.competitive.clean.csv'))
      res$Panel<-i
      dat<-rbind(dat, res)
    }

    dat$P.FDR<-p.adjust(dat$P, method='fdr')

    dat$NAME<-paste(toupper(substr(dat$NAME, 1, 1)), substr(dat$NAME, 2, nchar(dat$NAME)), sep="")
    dat$NAME<-gsub('\\.',' ',dat$NAME)

    # Format ATC codes and insert ATC descriptions
    tmp<-lapply(strsplit(dat$ATC, '\\.'), function(x) substr(x, 1, 4))
    tmp2<-lapply(tmp, insert_atc_desc, atc_labels)
    dat$ATC_code<-unlist(lapply(tmp, function(x) paste0(x, collapse=';')))
    dat$ATC_desc<-unlist(lapply(tmp2, function(x) paste0(x, collapse=';')))

    dat<-dat[order(dat$P),]

    dat$Panel<-tidy_panel_names(dat$Panel)

    dat<-dat[,c("NAME","Panel","N_Mem_Avail","Estimate","SE","P","P.FDR","ATC_code","ATC_desc") , with=F]
    names(dat)<-c('Name','Panel','N Genes','Estimate','SE','P','P.FDR','ATC Code','ATC Description')

    dat$ChEMBL<-paste0('<a href="https://www.ebi.ac.uk/chembl/g/#search_results/all/query=',dat$Name,'">','Link','</a>')

  }
  return(dat)
}

read_twas_gsea_atc<-function(config, gwas, mode = 'directional'){

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-NULL

  flag_param <- if(mode == 'nondirectional') 'twas_gsea_drugtargetor_nondirectional' else 'twas_gsea_drugtargetor'
  suffix <- if(mode == 'nondirectional') '_nondir' else ''

  if(read_param(config = config, param = flag_param, return_obj = F) == "T"){

    weights<-read.table(paste0(outdir,'/results/',gwas,'/twas/list_of_weights.txt'))$V1
    weights<-weights[!grepl('SPLIC',weights)]

    dat<-NULL
    for(i in weights){
      res<-fread(paste0(outdir,'/results/',gwas,'/twas/drugtargetor/twas_gsea',suffix,'_',i,'_res_atc_res.csv'))
      res$Panel<-i
      dat<-rbind(dat, res)
    }

    dat$P.FDR_all<-p.adjust(dat$P, method = 'fdr')

    dat<-dat[,c('Panel','ATC','Name','N','Estimate','Class_Median','Non_Class_Median','P','P.FDR_all'),with=T]
    names(dat)<-c('Panel','ATC Code','ATC Description','N Drugs','Estimate','Class Median T','Non-class Median T','P','P.FDR')

    dat<-dat[order(dat$P),]

    # Update labels
    dat$Panel<-tidy_panel_names(dat$Panel)

  }

  return(dat)
}

read_twas_gsea_cmap_drug<-function(config, gwas){
  # Per-signature CMAP TWAS-GSEA results, one row per
  # (cmap_name x cell_iname x pert_itime x pert_idose x weight panel).
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'twas_gsea_cmap', return_obj = F) == "T"){
    weights<-read.table(paste0(outdir,'/results/',gwas,'/twas/list_of_weights.txt'))$V1
    weights<-weights[!grepl('SPLIC',weights)]

    for(i in weights){
      f <- paste0(outdir,'/results/',gwas,'/twas/cmap/twas_gsea_cmap_',i,'_drug_res.csv')
      if(!file.exists(f)) next
      res<-fread(f)
      dat<-rbind(dat, res, fill = TRUE)
    }
    if(is.null(dat) || nrow(dat) == 0) return(NULL)

    # FDR across all (signature x panel) rows for cross-panel ranking.
    dat[, P.FDR_all := p.adjust(P, method = 'fdr')]
    dat[, Panel := tidy_panel_names(Panel)]
    dat[, Name := paste0(toupper(substr(cmap_name, 1, 1)), substr(cmap_name, 2, nchar(cmap_name)))]
    setorder(dat, P)

    dat <- dat[, .(Name, cmap_name, cell_iname, pert_itime, pert_idose, moa,
                   Panel, N_Mem_Avail, Estimate, SE, Z, P,
                   `P.FDR (per panel)` = P.FDR,
                   P.FDR = P.FDR_all)]
  }
  return(dat)
}

read_twas_gsea_cmap_moa<-function(config, gwas){
  # Per-MOA enrichment of CMAP TWAS-GSEA results.
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'twas_gsea_cmap', return_obj = F) == "T"){
    weights<-read.table(paste0(outdir,'/results/',gwas,'/twas/list_of_weights.txt'))$V1
    weights<-weights[!grepl('SPLIC',weights)]

    for(i in weights){
      f <- paste0(outdir,'/results/',gwas,'/twas/cmap/twas_gsea_cmap_',i,'_moa_res.csv')
      if(!file.exists(f)) next
      res<-fread(f)
      dat<-rbind(dat, res, fill = TRUE)
    }
    if(is.null(dat) || nrow(dat) == 0) return(NULL)

    dat[, P.FDR_all := p.adjust(P, method = 'fdr')]
    dat[, Panel := tidy_panel_names(Panel)]
    setorder(dat, P)

    dat <- dat[, .(Panel, MOA, Cell_Line, `N Drugs` = N, Estimate,
                   `Class Median T` = Class_Median,
                   `Non-class Median T` = Non_Class_Median,
                   P, P.FDR = P.FDR_all)]
  }
  return(dat)
}

read_magma_drug<-function(config, gwas){

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)
  resdir <- read_param(config = config, param = 'resdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'magma_drugtargetor', return_obj = F) == "T"){

    atc<-fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep='!')
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

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)
  resdir <- read_param(config = config, param = 'resdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'gcsc', return_obj = F) == "T"){

    atc<-fread(paste0(resdir, '/data/atc/atc_20220201.txt'), sep='!')
    names(atc)<-c('Code','Name')
    atc$Name<-tolower(atc$Name)

    atc_labels<-atc[nchar(atc$Code) == 4,]

    dat<-fread(paste0(outdir,'/results/',gwas,'/gcsc/',gwas,'_drugtargetor_gcsc_res.txt'))
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

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'magma_drugtargetor', return_obj = F) == "T"){

    dat<-fread(paste0(outdir,'/results/',gwas,'/magma/magma_drug_targetor_atc_res.csv'))
    dat$P.FDR<-p.adjust(dat$P, method = 'fdr')

    dat<-dat[,c("ATC","N","P","P.FDR","Name"), with=F]
    names(dat)<-c("ATC Code","N Drugs","P","P.FDR","ATC Description")

    dat<-dat[order(dat$P),]

  }
  return(dat)
}

read_gcsc_atc<-function(config, gwas){

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'gcsc', return_obj = F) == "T"){

    dat<-fread(paste0(outdir,'/results/',gwas,'/gcsc/',gwas,'_drugtargetor_gcsc_res_atc.csv'))

    dat$P.FDR<-p.adjust(dat$P, method = 'fdr')

    dat<-dat[,c("ATC","N","P","P.FDR","Name"), with=F]
    names(dat)<-c("ATC Code","N Drugs","P","P.FDR","ATC Description")

    dat<-dat[order(dat$P),]

  }
  return(dat)
}

read_magma_tissue<-function(config, gwas, type){

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)
  resdir <- read_param(config = config, param = 'resdir', return_obj = F)

  dat<-NULL

  if(read_param(config = config, param = 'tissue_magma', return_obj = F) == "T"){
    dat<-list()

    if(type == 'specific'){

      # Read in the MAGMA gene property enrichment results
      property_enrich<-fread(cmd=paste0("grep -v '^#' ",outdir,"/results/",gwas,'/magma/magma_tissue_spec.gsa.out'))

      # Insert FULL_NAME column if not present
      if(all(names(property_enrich) != 'FULL_NAME')){
          property_enrich$FULL_NAME<-property_enrich$VARIABLE
      }

      # Calculate FDR-corrected p-value
      property_enrich$P.FDR<-p.adjust(property_enrich$P, method = 'fdr')

      # Read in the tissue names
      tissue_groups<-fread(paste0(resdir, '/data/gtex/Tissue_labels.tsv'))

      # Insert original tissue names
      property_enrich<-merge(property_enrich, tissue_groups, by.x='FULL_NAME', by.y='new')

      # Remove unwanted columns
      property_enrich<-property_enrich[, c('Tissue','NGENES','BETA','SE','P','P.FDR'), with=F]
      names(property_enrich)<-c('Tissue','N Gene','BETA','SE','P','P.FDR')

      dat$res<-property_enrich

      # Read in list of retained tissues (file is only written when >=1
      # tissue is FDR-significant; treat absence as "no retained tissues")
      indep_file<-paste0(outdir,"/results/",gwas,'/magma/magma_tissue_conditional.indep.txt')
      if(file.exists(indep_file)){
        property_keep<-fread(indep_file, header=F)$V1
        property_keep<-tissue_groups$Tissue[tissue_groups$new %in% property_keep]
      } else {
        property_keep<-character(0)
      }

      dat$keep<-property_keep

    }
  }
  return(dat)
}
