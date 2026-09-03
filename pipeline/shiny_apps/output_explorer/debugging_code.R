################
# Code for debugging GenoDisc output explorer
################

gwas_data <- function(){
  rds<-readRDS("C:\\Users\\ollie\\OneDrive - King's College London\\KCL_Fellowship\\Programme\\GenoDisc\\shiny\\results_package.rds")
  rds
}

gwas_data <- function(){
  rds<-readRDS("C:\\Users\\ollie\\Downloads\\results_package.rds")
  rds
}

gwas_data <- function(){
  rds<-readRDS("~/oliverpainfel/Analyses/GenoDisc_ALS/output/results/results_package.rds")
  rds
}

input<-list()
input$gwas_selector<-'ALS_only'

selected_gwas<-reactive({
  req(gwas_data(), input$gwas_selector)
  gwas_selected <- ifelse(input$gwas_selector %in% names(gwas_data()),
                          input$gwas_selector,
                          names(gwas_data())[names(gwas_data()) != 'configuration'][1])
  return(gwas_selected)
})

mol_assoc_summary_data <- function(){

  all_func_res<-NULL

  if(finemap_logical){

    finemap_L1_tmp<-data.frame(Panel = "SuSie (L=1)",
                               ID=gwas_data()[[gwas_selected]]$mol_assoc$finemap$L1,
                               Z=1,
                               Sig=F,
                               Coloc=F,
                               Method="SNP\nFine-mapping",
                               Type='')

    all_func_res<-rbind(all_func_res, finemap_L1_tmp)

  }

  if(twas_logical){

    twas_tmp<-data.table(Panel=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$PANEL,
                         ID=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$`Gene Symbol`,
                         Z=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$TWAS.Z,
                         Sig=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$TWAS.P.FDR < 0.05,
                         Coloc=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$COLOC_logical)

    twas_tmp$Method<-'FUSION'
    twas_tmp$Type<-'Expr.'
    twas_tmp$Type[grepl('SPLIC',twas_tmp$Panel)]<-'Splice'

    # Retain only the most significant assoc for each gene within PANEL (only relevent for splice panel)
    twas_tmp<-twas_tmp[order(-abs(twas_tmp$Z)),]
    twas_tmp<-twas_tmp[!duplicated(paste0(twas_tmp$Panel, twas_tmp$ID)),]

    all_func_res<-rbind(all_func_res, twas_tmp)
  }

  if(smr_expression_logical){
    # SMR expression
    smr_expression_res_tmp<-data.table(Panel=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$PANEL,
                                       ID=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$`Gene Symbol`,
                                       Z=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$b_SMR/gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$se_SMR,
                                       Sig=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$p_SMR.FDR < 0.05,
                                       Coloc=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$p_HEIDI > 0.05)

    smr_expression_res_tmp$Method<-'SMR'
    smr_expression_res_tmp$Type<-'Expr.'

    all_func_res<-rbind(all_func_res, smr_expression_res_tmp)
  }

  # PWAS
  if(pwas_panel_rosmap_logical & !(pwas_panel_rosmap_logical & pwas_panel_banner_logical)){

    pwas_tmp<-data.table( Panel=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$PANEL,
                          ID=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$`Gene Symbol`,
                          Z=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$pwas_all.Z,
                          Sig=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$pwas_all.P.FDR < 0.05,
                          Coloc=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$COLOC_logical)

    pwas_tmp<-pwas_tmp[order(-abs(pwas_tmp$Z)),]
    pwas_tmp<-pwas_tmp[!duplicated(paste0(pwas_tmp$Panel, pwas_tmp$ID)),]
    pwas_tmp$Method<-'FUSION'
    pwas_tmp$Type<-'Protein'

    all_func_res<-rbind(all_func_res, pwas_tmp)
  }

  if(smr_protein_panel_rosmap_logical){

    # SMR protein
    smr_protein_res_tmp<-data.table(Panel=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$PANEL,
                                    ID=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$`Gene Symbol`,
                                    Z=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$b_SMR/gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$se_SMR,
                                    Sig=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$p_SMR.FDR < 0.05,
                                    Coloc=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$p_HEIDI > 0.05)

    smr_protein_res_tmp<-smr_protein_res_tmp[order(-abs(smr_protein_res_tmp$Z)),]
    smr_protein_res_tmp<-smr_protein_res_tmp[!duplicated(paste0(smr_protein_res_tmp$Panel, smr_protein_res_tmp$ID)),]
    smr_protein_res_tmp$Method<-'SMR'
    smr_protein_res_tmp$Type<-'Protein'

    all_func_res<-rbind(all_func_res, smr_protein_res_tmp)

  }

  if(magma_gene_logical){

    magma_tmp<-data.frame(Panel = 'MAGMA',
                          ID=gwas_data()[[gwas_selected]]$mol_assoc$magma$ID,
                          Z=abs(qnorm(as.numeric(gwas_data()[[gwas_selected]]$mol_assoc$magma$P))),
                          Sig=as.numeric(gwas_data()[[gwas_selected]]$mol_assoc$magma$P.FDR) < 0.05,
                          Coloc=F,
                          Method='MAGMA',
                          Type='')

    all_func_res<-rbind(all_func_res, magma_tmp)

  }

  if(clump_logical){

    nearest_tmp<-data.frame(Panel = 'NearestGene',
                            ID=gwas_data()[[gwas_selected]]$mol_assoc$nearest$clump,
                            Z=1,
                            Sig=F,
                            Coloc=F,
                            Method='Nearest\nGene',
                            Type='')

    all_func_res<-rbind(all_func_res, nearest_tmp)

  }

  return(all_func_res)
}

mol_assoc_summary_data_filtered<-function(){
  all_func_res<-mol_assoc_summary_data()

  # Filter results table by user specified methods
  all_func_res<-all_func_res[all_func_res$Method %in% input$selected_methods_mol,]

  # Filter results table by user specified expression and protein panels
  if(any(all_func_res$Type == 'Expr.')){
    all_func_res<-all_func_res[!(all_func_res$Type == 'Expr.' & !(all_func_res$Panel %in% input$selected_expr_panels_mol)),]
  }
  if(any(all_func_res$Type == 'Protein')){
    all_func_res<-all_func_res[!(all_func_res$Type == 'Protein' & !(all_func_res$Panel %in% input$selected_protein_panels_mol)),]
  }

  # Filter results table if user specifies high confidence genes only
  if(input$conf_only_mol){
    all_func_res<-all_func_res[all_func_res$ID %in% all_func_res$ID[which((all_func_res$Sig == T & all_func_res$Coloc == T) | all_func_res$Panel == "SuSie (L=1)")],]
  }

  input_genes <- unlist(strsplit(input$geneInput_mol, "[, ]"))
  selected_genes <- input_genes[input_genes != ""]

  # Insert NA rows for all panels and methods so when filtering by gene, all selected panels and methods remain
  na_rows<-all_func_res[!(duplicated(paste0(all_func_res$Panel, all_func_res$Method))),]
  na_rows$ID<-NA
  na_rows$Z<-NA
  na_rows$Sig<-NA
  na_rows$Coloc<-NA

  all_func_res<-rbind(na_rows, all_func_res)

  if(length(selected_genes) > 0){
    if(sum(selected_genes %in% all_func_res$ID) > 0){
      all_func_res<-all_func_res[all_func_res$ID %in% selected_genes | is.na(all_func_res$ID),]
    } else {
      all_func_res<-data.frame(matrix(nrow=0, ncol=5))
    }
  }

  return(all_func_res)
}

tx_drug_summary_data <- function(){

  ###
  # MAGMA
  ###
  magma_gs<-gwas_data()[[gwas_selected]]$tx$drug$magma

  # Convert one-sided p to a Z score
  magma_gs$Z<--qnorm(magma_gs$P)
  magma_gs<-magma_gs[,c('Name','Z','P','P.FDR','ATC Code')]
  magma_gs$Method<-'MAGMA'
  magma_gs$Panel<-'MAGMA'

  ###
  # GCSC
  ###
  gcsc_gs<-gwas_data()[[gwas_selected]]$tx$drug$gcsc

  gcsc_gs<-gcsc_gs[,c('Name','Z','P','P.FDR','ATC Code')]
  gcsc_gs$Method<-'GCSC'
  gcsc_gs$Panel<-'Brain and Blood'

  ###
  # TWAS-GSEA
  ###

  gsea_gs<-gwas_data()[[gwas_selected]]$tx$drug$twas_gsea
  gsea_gs$Method<-'TWAS-GSEA'
  gsea_gs$Z<-gsea_gs$Estimate/gsea_gs$SE
  gsea_gs<-gsea_gs[,c('Name','Z','P','P.FDR','Method','Panel','ATC Code')]
  # Flip Z so >0 indicates reversal of GWAS outcome.
  gsea_gs$Z<--gsea_gs$Z

  # Insert missing values
  gsea_gs_all<-gsea_gs
  for(i in unique(gsea_gs_all$Panel)){
    gsea_gs_i<-gsea_gs[gsea_gs$Panel == i,]
    gsea_gs_other<-gsea_gs[gsea_gs$Panel != i,]
    gsea_gs_rest<-data.frame(Name=unique(gsea_gs_other$Name[!(gsea_gs_other$Name %in% gsea_gs_i$Name)]),
                             Z=NA,
                             P=NA,
                             P.FDR=NA,
                             Method='TWAS-GSEA',
                             Panel=i,
                             ATC_Code=NA)
    names(gsea_gs_rest)<-gsub('ATC_Code','ATC Code',names(gsea_gs_rest))

    gsea_gs_all<-rbind(gsea_gs_all, gsea_gs_rest)

  }
  gsea_gs<-gsea_gs_all

  ###
  # Combine results
  ###

  all_gs<-do.call(rbind, list(magma_gs, gcsc_gs, gsea_gs))

  return(all_gs)
}

tx_drug_summary_data_filtered<-function(){
  all_gs<-tx_drug_summary_data()

  # Filter results table by user specified methods
  all_gs<-all_gs[all_gs$Method %in% input$selected_methods_drug,]

  # Filter results table by user specified expression
  if(any(all_gs$Method == 'TWAS-GSEA')){
    all_gs<-all_gs[!(all_gs$Method == 'TWAS-GSEA' & !(all_gs$Panel %in% input$selected_expr_panels_drug)),]
  }

  # Filter results table if user specifies high confidence genes only
  if(input$conf_only_drug){
    all_gs<-all_gs[all_gs$Name %in% all_gs$Name[which(all_gs$P.FDR < 0.05)],]
  }

  input_drugs <- unlist(strsplit(input$drugInput_drug, "[, ]"))
  selected_drugs <- input_drugs[input_drugs != ""]

  input_atc <- unlist(strsplit(input$atcInput_drug, "[, ]"))
  selected_atc <- input_atc[input_atc != ""]

  # Insert NA rows for all panels and methods so when filtering by drug, all selected panels and methods remain
  na_rows<-all_gs[!(duplicated(paste0(all_gs$Panel, all_gs$Method))),]
  na_rows$Name<-NA
  na_rows$Z<-NA
  na_rows$P<-NA
  na_rows$P.FDR<-NA

  all_gs<-rbind(na_rows, all_gs)

  if(length(selected_drugs) > 0){
    if(sum(grepl(paste(selected_drugs, collapse='|'), all_gs$Name, ignore.case = T)) > 0){
      all_gs<-all_gs[grepl(paste(selected_drugs, collapse='|'), all_gs$Name, ignore.case = T) | is.na(all_gs$Name),]
    } else {
      all_gs<-data.frame(matrix(nrow=0, ncol=5))
    }
  }

  if(length(selected_atc) > 0){
    if(sum(grepl(paste(selected_atc, collapse='|'), all_gs$`ATC Code`, ignore.case = T)) > 0){
      all_gs<-all_gs[grepl(paste(selected_atc, collapse='|'), all_gs$`ATC Code`, ignore.case = T) | is.na(all_gs$`ATC Code`),]
    } else {
      all_gs<-data.frame(matrix(nrow=0, ncol=5))
    }
  }

  return(all_gs)
}

plot_dim_drug<-function(){
  all_gs<-tx_drug_summary_data_filtered()

  if(nrow(all_gs) > 0){
    num_row <- length(unique(all_gs$Name))
    plot_height<-(max(nchar(all_gs$Panel))*3)+(num_row * 20)+100
    num_col <- length(unique(paste0(all_gs$Panel,'_',all_gs$Method,'_')))
    plot_width<-(max(nchar(all_gs$Name), na.rm=T)*1.2)+(num_col * 60)
    plot_width<-max(plot_width,(length(unique(all_gs$Method))*100))
  } else {
    plot_height<-100
    plot_width<-100
  }

  return(list(height=plot_height,
              width=plot_width))
}

tx_atc_summary_data <- function(){

  ###
  # MAGMA
  ###

  magma_gs_atc<-gwas_data()[[gwas_selected]]$tx$atc$magma

  magma_gs_atc$Z<--qnorm(magma_gs_atc$P)
  magma_gs_atc$FDR_Sig<-magma_gs_atc$P.FDR < 0.05
  magma_gs_atc$Nom_Sig<-magma_gs_atc$P < 0.05
  magma_gs_atc$Name<-paste0(magma_gs_atc$`ATC Code`,': ',magma_gs_atc$`ATC Description`)
  magma_gs_atc$Method<-'MAGMA'
  magma_gs_atc$Panel<-'MAGMA'

  magma_gs_atc<-magma_gs_atc[,c("Name","Z","FDR_Sig","Nom_Sig","Method","Panel"), with=F]

  ###
  # GCSC
  ###

  gcsc_gs_atc<-gwas_data()[[gwas_selected]]$tx$atc$gcsc

  gcsc_gs_atc$Z<--qnorm(gcsc_gs_atc$P)
  gcsc_gs_atc$FDR_Sig<-gcsc_gs_atc$P.FDR < 0.05
  gcsc_gs_atc$Nom_Sig<-gcsc_gs_atc$P < 0.05
  gcsc_gs_atc$Name<-paste0(gcsc_gs_atc$`ATC Code`,': ',gcsc_gs_atc$`ATC Description`)
  gcsc_gs_atc$Method<-'GCSC'
  gcsc_gs_atc$Panel<-'GCSC'

  gcsc_gs_atc<-gcsc_gs_atc[,c("Name","Z","FDR_Sig","Nom_Sig","Method","Panel"), with=F]

  ###
  # TWAS-GSEA
  ###

  gsea_gs_atc<-gwas_data()[[gwas_selected]]$tx$atc$twas_gsea

  gsea_gs_atc$P.FDR_all<-p.adjust(gsea_gs_atc$P, method = 'fdr')
  gsea_gs_atc$P.FDR.onside_all<-p.adjust(gsea_gs_atc$P.oneside, method = 'fdr')

  gsea_gs_atc$Z<--qnorm(gsea_gs_atc$P)
  gsea_gs_atc$Z<-gsea_gs_atc$Z*sign(gsea_gs_atc$Estimate)
  gsea_gs_atc$FDR_Sig<-gsea_gs_atc$P.FDR_all < 0.05
  gsea_gs_atc$Nom_Sig<-gsea_gs_atc$P < 0.05
  gsea_gs_atc$Name<-paste0(gsea_gs_atc$`ATC Code`,': ',gsea_gs_atc$`ATC Description`)

  gsea_gs_atc$Method<-'TWAS-GSEA'

  gsea_gs_atc<-gsea_gs_atc[,c("Name","Z","FDR_Sig","Nom_Sig","Method","Panel"), with=F]

  # Insert missing values
  gsea_gs_atc_all<-gsea_gs_atc
  for(i in unique(gsea_gs_atc_all$Panel)){
    gsea_gs_atc_i<-gsea_gs_atc[gsea_gs_atc$Panel == i,]
    gsea_gs_atc_other<-gsea_gs_atc[gsea_gs_atc$Panel != i,]
    gsea_gs_atc_rest<-data.frame(Name=unique(gsea_gs_atc_other$Name[!(gsea_gs_atc_other$Name %in% gsea_gs_atc_i$Name)]),
                                 Z =NA,
                                 FDR_Sig=NA,
                                 Nom_Sig=NA,
                                 Method='TWAS-GSEA',
                                 Panel=i)

    gsea_gs_atc_all<-rbind(gsea_gs_atc_all, gsea_gs_atc_rest)

  }
  gsea_gs_atc<-gsea_gs_atc_all

  all_gs_atc<-do.call(rbind, list(magma_gs_atc, gcsc_gs_atc, gsea_gs_atc))

  return(all_gs_atc)
}

tx_atc_summary_data_filtered<-function(){
  all_gs_atc<-tx_atc_summary_data()

  # Filter results table by user specified methods
  all_gs_atc<-all_gs_atc[all_gs_atc$Method %in% input$selected_methods_atc,]

  # Filter results table by user specified expression
  if(any(all_gs_atc$Method == 'TWAS-GSEA')){
    all_gs_atc<-all_gs_atc[!(all_gs_atc$Method == 'TWAS-GSEA' & !(all_gs_atc$Panel %in% input$selected_expr_panels_atc)),]
  }

  # Filter results table if user specifies high confidence genes only
  if(input$conf_only_atc){
    all_gs_atc<-all_gs_atc[all_gs_atc$Name %in% all_gs_atc$Name[which(all_gs_atc$FDR_Sig == T)],]
  }

  input_atcs <- unlist(strsplit(input$atcInput_atc, "[, ]"))
  selected_atcs <- input_atcs[input_atcs != ""]

  # Insert NA rows for all panels and methods so when filtering by atc, all selected panels and methods remain
  na_rows<-all_gs_atc[!(duplicated(paste0(all_gs_atc$Panel, all_gs_atc$Method))),]
  na_rows$Name<-NA
  na_rows$Z<-NA
  na_rows$FDR_Sig<-NA
  na_rows$Nom_Sig<-NA

  all_gs_atc<-rbind(na_rows, all_gs_atc)

  if(length(selected_atcs) > 0){
    if(sum(grepl(paste(selected_atcs, collapse='|'), all_gs_atc$Name, ignore.case = T)) > 0){
      all_gs_atc<-all_gs_atc[grepl(paste(selected_atcs, collapse='|'), all_gs_atc$Name, ignore.case = T) & !is.na(all_gs_atc$Name),]
    } else {
      all_gs_atc<-data.frame(matrix(nrow=0, ncol=5))
    }
  }

  return(all_gs_atc)
}

read_ss<-function(){
  req(input$sumstats)

  print(input$sumstats$datapath)

  # Read in the header and interpret column names
  sub_ss<-fread(input$sumstats$datapath, nrows = 1000)

  return(sub_ss)
}

head_interp<-function(){
  sub_ss<-read_ss()

  print(head(sub_ss))

  # Read in the header and interpret column names
  sub_header<-toupper(names(sub_ss))

  # Remove columns that are all NA
  sub_ss_comp<-sub_ss[,apply(sub_ss, 2, function(x) !(all(is.na(x)))), with=F]
  sub_header_comp<-toupper(names(sub_ss_comp))

  int_header <- sub_header_comp
  for(i in names(ss_head_dict)){
    int_header[int_header %in% ss_head_dict[[i]]] <- i
  }
  int_header[!(int_header %in% unlist(ss_head_dict))]<-NA

  # Show original and interpreted header
  header_interp <- data.frame(Original = sub_header_comp,
                              Interpreted = int_header)

  # Show columns that are ignored due to be irrelevant, duplicated, or all NA
  header_interp$Keep<-!(
    !(int_header %in% names(ss_head_dict)) |
      duplicated(int_header)
  )

  # Insert reason it was ignored
  header_interp$Reason<-NA
  header_interp$Reason[!(int_header %in% names(ss_head_dict))]<-'Not recognised'
  header_interp$Reason[duplicated(int_header)]<-'Duplicated'

  # Show columns ignored due to missingness
  for(i in sub_header[!(sub_header %in% header_interp$Original)]){
    header_interp<-rbind(header_interp, data.frame(Original = i,
                                                   Interpreted = NA,
                                                   Keep=F,
                                                   Reason = 'First 1000 rows NA'))
  }

  header_interp[is.na(header_interp)]<-'NA'
  header_interp$Keep<-as.character(header_interp$Keep)

  # Remove ignored columns
  header_interp_keep<-header_interp[header_interp$Keep == T,]
  header_interp_keep$Keep<-NULL
  header_interp_keep$Reason<-NULL

  # Insert description of each column after interpretation
  header_labels<-data.frame(Interpreted=c('SNP','CHR','BP','A1','A2','P','OR','BETA','Z','SE','N','N_CAS','N_CON','NEF','FRQ','FRQ_A','FREQ_U','INFO'),
                            Description=c("RSID for variant",
                                          "Chromosome number",
                                          "Base pair position",
                                          "Allele 1 (effect allele)",
                                          "Allele 2",
                                          "P-value of association",
                                          "Odds ratio effect size",
                                          "BETA effect size",
                                          "Z-score",
                                          "Standard error of log(OR) or BETA",
                                          "Total sample size",
                                          "Number of cases",
                                          "Number of controls",
                                          "Effective sample size",
                                          "Allele frequency",
                                          "Allele frequency in cases",
                                          "Allele frequency in controls",
                                          "Imputation quality"))

  header_interp_keep<-merge(header_interp_keep, header_labels, by='Interpreted')
  header_interp_keep<-header_interp_keep[match(names(ss_head_dict), header_interp_keep$Interpreted),]
  header_interp_keep<-header_interp_keep[complete.cases(header_interp_keep),]
  header_interp_keep<-header_interp_keep[,c('Original','Interpreted','Description')]

  return(list(
    sub_ss_head=header_interp_keep,
    sub_ss_head_full=header_interp
  ))
}

input<-list()
input$sumstats$datapath<-"C:\\Users\\ollie\\Downloads\\23358156-GCST001837-EFO_0004337.h.tsv.gz"
input$selected_methods_drug<-c('MAGMA','GCSC')
input$selected_expr_panels_drug<-c("Brain - DLPFC (PsychENCODE)", "Brain - Frontal cortex BA9 (GTEx)")
input$conf_only_drug<-F
input$drugInput_drug<-''
input$selected_sort_drug <-'All-GSEA - Z'
input$atcInput_drug<-''
input$atcInput_drug<-'C10A'

input$selected_methods_atc<-c('MAGMA','GCSC')
input$selected_expr_panels_atc<-c("Brain - DLPFC (PsychENCODE)", "Brain - Frontal cortex BA9 (GTEx)")
input$conf_only_atc<-F
input$selected_sort_atc <-'All - Z'
input$atcInput_atc<-''

input$selected_methods_mol<-unique(mol_assoc_summary_data()$Method)
input$conf_only_mol<-T
input$geneInput_mol<-''
input$selected_expr_panels_mol<-unique(mol_assoc_summary_data()$Panel[mol_assoc_summary_data()$Type == 'Expr.' | mol_assoc_summary_data()$Type == 'Splice'])
input$selected_protein_panels_mol<-unique(mol_assoc_summary_data()$Panel[mol_assoc_summary_data()$Type == 'Protein'])

# Read in config file
config<-gwas_data()$configuration$config

# Read in GWAS list
gwas_list<-gwas_data()$configuration$gwas_list

# Determine which analyses were requested
clump_logical<-config[grepl('clump:',config)] == "clump: T"

cojo_logical<-config[grepl('cojo:',config)] == "cojo: T"

finemap_logical<-config[grepl('finemap:',config)] == "finemap: T"

ldsc_logical<-config[grepl('ldsc:',config)] == "ldsc: T"

magma_gene_logical<-config[grepl('magma_gene:',config)] == "magma_gene: T"

twas_panel_psychencode_logical<-config[grepl('twas_panel_psychencode:',config)] == "twas_panel_psychencode: T"
twas_panel_fusion_logical<-config[grepl('twas_panel_fusion:',config)] == "twas_panel_fusion: T"

twas_logical<-any(twas_panel_psychencode_logical, twas_panel_fusion_logical)

twas_conditional_logical<-config[grepl('twas_conditional:',config)] == "twas_conditional: T" & twas_logical == T

smr_expression_panel_psychencode_logical<-config[grepl('smr_expression_panel_psychencode:',config)] == "smr_expression_panel_psychencode: T"

smr_expression_panel_metabrain_basalganglia_logical<-config[grepl('smr_expression_panel_metabrain_basalganglia:',config)] == "smr_expression_panel_metabrain_basalganglia: T"
smr_expression_panel_metabrain_cerebellum_logical<-config[grepl('smr_expression_panel_metabrain_cerebellum:',config)] == "smr_expression_panel_metabrain_cerebellum: T"
smr_expression_panel_metabrain_cortex_logical<-config[grepl('smr_expression_panel_metabrain_cortex:',config)] == "smr_expression_panel_metabrain_cortex: T"
smr_expression_panel_metabrain_hippocampus_logical<-config[grepl('smr_expression_panel_metabrain_hippocampus:',config)] == "smr_expression_panel_metabrain_hippocampus: T"
smr_expression_panel_metabrain_spinalcord_logical<-config[grepl('smr_expression_panel_metabrain_spinalcord:',config)] == "smr_expression_panel_metabrain_spinalcord: T"

smr_expression_panel_eqtlgen_logical<-config[grepl('smr_expression_panel_eqtlgen:',config)] == "smr_expression_panel_eqtlgen: T"

pwas_panel_rosmap_logical<-config[grepl('pwas_panel_rosmap:',config)] == "pwas_panel_rosmap: T"
pwas_panel_banner_logical<-config[grepl('pwas_panel_banner:',config)] == "pwas_panel_banner: T"

smr_protein_panel_rosmap_logical<-config[grepl('smr_protein_panel_rosmap:',config)] == "smr_protein_panel_rosmap: T"

magma_drugtargetor_logical<-config[grepl('magma_drugtargetor:',config)] == "magma_drugtargetor: T"
twas_gsea_lincs_logical<-config[grepl('twas_gsea_lincs:',config)] == "twas_gsea_lincs: T"
twas_so_lincs_logical<-config[grepl('twas_so_lincs:',config)] == "twas_so_lincs: T"

dgi_db_comp_logical<-config[grepl('dgi_db_comp:',config)] == "dgi_db_comp: T"

twas_gsea_drugtargetor_logical<-config[grepl('twas_gsea_drugtargetor:',config)] == "twas_gsea_drugtargetor: T"

# Summarise this information
mol_assoc_logical<-any(magma_gene_logical,
                       twas_logical,
                       smr_expression_panel_psychencode_logical,
                       smr_expression_panel_metabrain_basalganglia_logical,
                       smr_expression_panel_metabrain_cerebellum_logical,
                       smr_expression_panel_metabrain_cortex_logical,
                       smr_expression_panel_metabrain_hippocampus_logical,
                       smr_expression_panel_metabrain_spinalcord_logical,
                       smr_expression_panel_eqtlgen_logical,
                       pwas_panel_rosmap_logical,
                       pwas_panel_banner_logical,
                       smr_protein_panel_rosmap_logical)

metabrain_logical<-any(smr_expression_panel_metabrain_basalganglia_logical,
                       smr_expression_panel_metabrain_cerebellum_logical,
                       smr_expression_panel_metabrain_cortex_logical,
                       smr_expression_panel_metabrain_hippocampus_logical,
                       smr_expression_panel_metabrain_spinalcord_logical)

smr_expression_logical<-any(smr_expression_panel_psychencode_logical,
                            metabrain_logical,
                            smr_expression_panel_eqtlgen_logical)

expression_logical<-any( twas_logical,
                         smr_expression_logical)

psychencode_logical<-any(twas_panel_psychencode_logical,
                         smr_expression_panel_psychencode_logical)

protein_logical<-any( pwas_panel_rosmap_logical,
                      pwas_panel_banner_logical,
                      smr_protein_panel_rosmap_logical)

drug_logical<-any(dgi_db_comp_logical,
                  magma_drugtargetor_logical,
                  twas_gsea_lincs_logical,
                  twas_so_lincs_logical,
                  twas_gsea_drugtargetor_logical)

# Assign select GWAS name to variable
gwas_selected <- input$gwas_selector

