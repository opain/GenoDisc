#!/usr/bin/Rscript

process_cleaner_log<-function(outdir, gwas){
    dat<-list()

    dat$log<-readLines(paste0(outdir,'/data/gwas_sumstat/',gwas,'/',gwas,'.cleaned.log'))

    dat$val<-list()
    dat$val$n_snp_orig<-as.numeric(gsub(' variants.','',gsub('GWAS contains ','',dat$log[grepl('GWAS contains ', dat$log)])))
    dat$val$build<-extract_build(dat$log)
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

process_focus_log<-function(outdir, gwas){
    dat<-list()

    dat$log<-readLines(paste0(outdir,'/data/gwas_sumstat/',gwas,'/',gwas,'.cleaned.munged.log'))

    dat$val<-list()
    dat$val$lambda_gc<-as.numeric(gsub('.*Lambda GC = ','',dat$log[grepl('Lambda GC', dat$log)]))
    dat$val$max_chi2<-as.numeric(gsub('.*Max chi\\^2 = ','',dat$log[grepl('Max chi\\^2', dat$log)]))
    dat$val$n_sig_snp<-as.numeric(gsub('.* ','',gsub(' Genome-wide significant SNPs.*','',dat$log[grepl('Genome-wide significant SNPs', dat$log)])))

    return(dat)
}

process_ldsc_log<-function(outdir, gwas){
    dat<-list()

    dat$log<-readLines(paste0(outdir,'/results/',gwas,'/ldsc/',gwas,'_ldsc_res.log'))

    dat$val<-list()
    tmp<-gsub('Total Observed scale h2: ','',dat$log[grepl('Total Observed scale h2:', dat$log)])
    dat$val$obs_h2_est<-as.numeric(gsub(' .*','',tmp))
    dat$val$obs_h2_se<-as.numeric(gsub('\\)','',gsub('.*\\(','',tmp)))

    tmp<-gsub('Intercept: ','',dat$log[grepl('Intercept: ', dat$log)])
    dat$val$int_est<-as.numeric(gsub(' .*','',tmp))
    dat$val$int_se<-as.numeric(gsub('\\)','',gsub('.*\\(','',tmp)))

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

read_fusion_exp<-function(config, gwas, outdir){

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
    panels<-NULL
    for(i in twas_weights){
        pos<-fread(paste0('resources/data/fusion_snp_weights/',i,'/',i,'.pos'))
        panels<-rbind(panels, data.frame( Type='Expression',
                                        Software='FUSION',
                                        Panel=i,
                                        N_indiv=median(pos$N),
                                        N_gene=nrow(pos)))
    }

    # Update panel names to be easier to read
    panels$Panel<-tidy_panel_names(panels$Panel)

    # Read in TWAS results
    res<-fread(paste0(outdir,'/results/',gwas,'/twas/',gwas,'_twas_GW_clean.txt.gz'))
    res$FILE<-NULL
    res$MODEL<-NULL
    res<-res[!is.na(res$TWAS.Z),]

    # Calculate P accounting for multiple testing
    res$TWAS.P.FDR<-p.adjust(res$TWAS.P, method='fdr')
    res<-res[,c('CHR','P0','P1','ensembl_gene_id','external_gene_name','PANEL','TWAS.Z','TWAS.P','TWAS.P.FDR','COLOC.PP3','COLOC.PP4'), with=T]
    names(res)[names(res) == 'ensembl_gene_id']<-'Ensembl ID'
    names(res)[names(res) == 'external_gene_name']<-'Gene Symbol'
    res$COLOC_logical<-(res$COLOC.PP4-res$COLOC.PP3)/res$COLOC.PP4 > 0.8

    res<-res[order(res$CHR, res$P0),]

    res$PANEL<-tidy_panel_names(res$PANEL)

    return(list(panels=panels,
                results=res))
}

read_smr_exp<-function(config, gwas, outdir){
  # Check what was run
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

  return(list(panels=panels,
              results=res))

}
