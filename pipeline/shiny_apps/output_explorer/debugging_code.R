################
# Code for debugging GenoDiscover output explorer
################

gwas_data <- function(){
  rds<-readRDS("C:\\Users\\ollie\\OneDrive - King's College London\\KCL_Fellowship\\Programme\\GenoDiscover\\shiny\\results_package.rds")
  rds
}

gwas_selected<-'COAD01'

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

input<-list()
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
