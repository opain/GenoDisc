gwasQcUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title="GWAS QC",
    br(),
    p("This tab shows key quality control statistics for your selected GWAS."),
    fluidRow(
      column(width=4,
             dataTableOutput(ns("qc_table")),
      )
    )
  )
}

gwasQcServer <- function(id, gwas_data, selected_gwas, gwas_list) {
  moduleServer(id, function(input, output, session) {

    # Create a table showing key statistics
    qc_val <- reactive({
      req(gwas_data(), selected_gwas(), gwas_list())

      qc_val<-data.table(name=selected_gwas(),
                         label=gwas_list()$label[gwas_list()$name == selected_gwas()],
                         n_var_orig=gwas_data()[[selected_gwas()]]$gwas_qc$cleaner_dat$val$n_var_orig,
                         build=gwas_data()[[selected_gwas()]]$gwas_qc$cleaner_dat$val$build$build,
                         n_snp_final=gwas_data()[[selected_gwas()]]$gwas_qc$cleaner_dat$val$n_snp_final,
                         lambda_gc=gwas_data()[[selected_gwas()]]$gwas_qc$focus_dat$val$lambda_gc,
                         max_chi2=gwas_data()[[selected_gwas()]]$gwas_qc$focus_dat$val$max_chi2,
                         n_sig_snp=gwas_data()[[selected_gwas()]]$gwas_qc$focus_dat$val$n_sig_snp,
                         obs_h2=paste0(round(gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_dat$val$obs_h2_est,3), " (",round(gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_dat$val$obs_h2_se,3),")"),
                         int=paste0(round(gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_dat$val$int_est,3), " (",round(gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_dat$val$int_se,3),")"))

      names(qc_val)<-c('GWAS Name',
                       'GWAS Label',
                       'N variants pre-QC',
                       'Identified genome build',
                       'N variants post-QC',
                       'Lambda GC',
                       'Max. chi^2',
                       'N genome-wide significant variants',
                       "LDSC SNP-heritability (SE; observed scale)",
                       "LDSC intercept (SE)")

      qc_val<-t(qc_val)
      qc_val<-data.table(Parameter=dimnames(qc_val)[[1]],
                         Value=qc_val[,1])

      qc_val
    })

    output$qc_table <- renderDataTable({
      datatable(qc_val(), options = list(dom = 't',
                                     ordering=F),
                selection = 'none',
                rownames = F,
                colnames = '') %>%
        formatStyle(columns = c("Parameter"), fontWeight = 'bold', textAlign = "center") %>%
        formatStyle(columns = c("Value"), textAlign = "center")
    })
  })
}
