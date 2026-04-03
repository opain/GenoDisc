gwasQcUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title="GWAS QC",
    br(),
    p("This tab shows key quality control statistics for your selected GWAS."),

    # --- Report Card: Top Row ---
    fluidRow(
      # Left column: Summary table card
      column(
        width = 4,
        div(
          class = "panel panel-default",
          div(class = "panel-heading", tags$strong("QC Summary")),
          div(class = "panel-body", dataTableOutput(ns("qc_table")))
        )
      ),
      # Right column: MAF plot card
      column(
        width = 8,
        div(
          class = "panel panel-default",
          div(class = "panel-heading", tags$strong("Allele Frequency Plot")),
          div(class = "panel-body", uiOutput(ns("maf_plot_ui")))
        )
      )
    ),

    # --- Technical Appendix ---
    br(),
    div(
      class = "panel panel-default",
      div(class = "panel-heading", tags$strong("Technical Appendix: Sumstat Cleaner Log")),
      div(class = "panel-body", uiOutput(ns("cleaner_log_ui")))
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

    # MAF plot rendering
    output$maf_plot_ui <- renderUI({
      req(gwas_data(), selected_gwas())
      b64 <- gwas_data()[[selected_gwas()]]$gwas_qc$maf_plot_base64

      if (!is.null(b64)) {
        tags$img(
          src = paste0("data:image/png;base64,", b64),
          style = "max-width: 50%; height: auto;"
        )
      } else {
        div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          icon("chart-bar", style = "font-size: 2em;"),
          br(), br(),
          tags$strong("MAF Plot Unavailable"),
          br(),
          "Allele frequency data was not provided in the input GWAS."
        )
      }
    })

    # Cleaner log rendering
    output$cleaner_log_ui <- renderUI({
      req(gwas_data(), selected_gwas())
      log_lines <- gwas_data()[[selected_gwas()]]$gwas_qc$cleaner_dat$log

      if (is.null(log_lines) || length(log_lines) == 0) {
        return(p("Log not available."))
      }

      pre(
        style = "background-color: #1e1e1e; color: #00ff00; font-family: 'Courier New', monospace; font-size: 12px; padding: 15px; border-radius: 4px; max-height: 400px; overflow-y: auto; white-space: pre-wrap;",
        paste(log_lines, collapse = "\n")
      )
    })
  })
}
