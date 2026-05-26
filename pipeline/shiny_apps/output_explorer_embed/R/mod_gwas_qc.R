gwasQcUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "GWAS QC",
    br(),
    p("This tab shows key quality control statistics for your selected GWAS."),
    hr(),
    tabsetPanel(
      tabPanel(
        title = "QC Summary",
        br(),
        div(style = "max-width: 700px;", dataTableOutput(ns("qc_table")))
      ),
      tabPanel(
        title = "Allele Frequency Plot",
        br(),
        uiOutput(ns("maf_plot_ui"))
      ),
      tabPanel(
        title = "QQ Plot",
        br(),
        uiOutput(ns("qq_plot_ui"))
      ),
      tabPanel(
        title = "Sumstat QC Log",
        br(),
        uiOutput(ns("cleaner_log_ui"))
      )
    )
  )
}

gwasQcServer <- function(id, gwas_data, selected_gwas, gwas_list, config_flags) {
  moduleServer(id, function(input, output, session) {

    # Create a table showing key statistics
    qc_val <- reactive({
      req(gwas_data(), selected_gwas(), gwas_list(), config_flags())

      g <- gwas_data()[[selected_gwas()]]
      cf <- config_flags()

      qc_val<-data.table(name=selected_gwas(),
                         label=gwas_list()$label[gwas_list()$name == selected_gwas()],
                         n_var_orig=g$gwas_qc$cleaner_dat$val$n_var_orig,
                         build=g$gwas_qc$cleaner_dat$val$build$build,
                         n_snp_final=g$gwas_qc$cleaner_dat$val$n_snp_final,
                         lambda_gc=g$gwas_qc$focus_dat$val$lambda_gc,
                         max_chi2=g$gwas_qc$focus_dat$val$max_chi2,
                         n_sig_snp=g$gwas_qc$focus_dat$val$n_sig_snp)

      col_labels <- c('GWAS Name',
                      'GWAS Label',
                      'N variants pre-QC',
                      'Identified genome build',
                      'N variants post-QC',
                      'Lambda GC',
                      'Max. chi^2',
                      'N genome-wide significant variants')

      if (isTRUE(cf$ldsc)) {
        qc_val[, obs_h2 := paste0(round(g$gwas_qc$ldsc_dat$val$obs_h2_est,3), " (",round(g$gwas_qc$ldsc_dat$val$obs_h2_se,3),")")]
        qc_val[, int := paste0(round(g$gwas_qc$ldsc_dat$val$int_est,3), " (",round(g$gwas_qc$ldsc_dat$val$int_se,3),")")]
        col_labels <- c(col_labels,
                        "LDSC SNP-heritability (SE; observed scale)",
                        "LDSC intercept (SE)")
      }

      names(qc_val) <- col_labels

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
        tagList(
          tags$img(
            src = paste0("data:image/png;base64,", b64),
            style = "max-height: 430px; width: auto;"
          ),
          tags$p(
            style = "font-size: 0.85em; color: #6c757d; margin-top: 8px;",
            "Only variants with an allele frequency difference greater than 0.2 from the reference are shown; these are the variants removed by the QC script."
          )
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

    # QQ plot rendering
    output$qq_plot_ui <- renderUI({
      req(gwas_data(), selected_gwas())
      b64 <- gwas_data()[[selected_gwas()]]$gwas_qc$qq_plot_base64

      if (!is.null(b64)) {
        tags$img(
          src = paste0("data:image/png;base64,", b64),
          style = "max-height: 500px; width: auto;"
        )
      } else {
        div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          icon("chart-line", style = "font-size: 2em;"),
          br(), br(),
          tags$strong("QQ Plot Unavailable"),
          br(),
          "This results package was produced before the QQ plot rule was added."
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
