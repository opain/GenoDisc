
gwasQcUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "GWAS QC",
    br(),
    p("This tab shows key quality control statistics for your selected GWAS."),
    hr(),
    tabsetPanel(
      id = ns("gwas_qc_tabs"),
      tabPanel(
        title = "QC Summary",
        value = "qc_summary",
        br(),
        div(style = "max-width: 700px;", tableOutput(ns("qc_table"))),
        uiOutput(ns("qc_legend"))
      ),
      tabPanel(
        title = "Allele Frequency Plot",
        value = "maf_plot",
        br(),
        uiOutput(ns("maf_plot_ui"))
      ),
      tabPanel(
        title = "QQ Plot",
        value = "qq_plot",
        br(),
        uiOutput(ns("qq_plot_ui"))
      ),
      # SNP-h² and rG have moved to their own top-level tab
      # (SNP-h² & rG). Removing them keeps this tab focused on QC.
      tabPanel(
        title = "Sumstat QC Log",
        value = "cleaner_log",
        br(),
        uiOutput(ns("cleaner_log_ui"))
      )
    )
  )
}

gwasQcServer <- function(id, gwas_data, selected_gwas, gwas_list, config_flags,
                          selected_gwas_multi = NULL,
                          comparison_mode = NULL) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # Every sub-tab in this module is single-GWAS-only. The outer GWAS QC
    # tab itself is hidden in compare mode from app.R (users drop back to
    # single-GWAS mode to inspect these); no per-sub-tab hide is needed.

    # Plain-text explanations for the QC metrics, used by the legend under the table.
    qc_help <- list(
      "N variants pre-QC" = "Number of variants in the input GWAS before quality control.",
      "Identified genome build" = paste0(
        "Genome build (e.g. GRCh37/GRCh38) inferred from the variant positions; ",
        "shown as 'Unknown' if it could not be determined."),
      "N variants post-QC" = "Number of variants remaining after quality-control filtering.",
      "Lambda GC" = paste0(
        "Genomic inflation factor: median chi-square divided by 0.455. About 1.0 means the ",
        "test statistics are well calibrated; much above 1.1 can indicate confounding (or, ",
        "for a well-powered polygenic trait, genuine widespread signal)."),
      "Max. chi^2" = paste0(
        "The strongest single-variant association in the study, as a chi-square statistic ",
        "derived from the p-value. Larger values mean a stronger top signal."),
      "N genome-wide significant variants" = paste0(
        "Number of variants reaching genome-wide significance (p < 5e-8) after QC. ",
        "Zero is common for underpowered GWAS.")
    )

    # Create a table showing key statistics
    qc_val <- reactive({
      req(gwas_data(), selected_gwas(), gwas_list(), config_flags())

      gwas_qc <- gd_read(gwas_data(), selected_gwas(), "gwas_qc")
      cf <- config_flags()

      # Build isn't always identified from CHR/BP (the cleaner can fall back to
      # matching by SNP ID instead, which never determines a build - see
      # extract_build() in package_results_functions.R). Coalesce to a scalar
      # so a single missing QC field can't collapse the whole table via
      # data.table()'s zero-length recycling.
      build_val <- gwas_qc$cleaner_dat$val$build$build
      if (length(build_val) == 0 || is.na(build_val)) build_val <- "Unknown"

      qc_val<-data.table(name=selected_gwas(),
                         label=gwas_list()$label[gwas_list()$name == selected_gwas()],
                         n_var_orig=gwas_qc$cleaner_dat$val$n_var_orig,
                         build=build_val,
                         n_snp_final=gwas_qc$cleaner_dat$val$n_snp_final,
                         lambda_gc=gd_qc_stat(gwas_qc, "lambda_gc"),
                         max_chi2=gd_qc_stat(gwas_qc, "max_chi2"),
                         n_sig_snp=gd_qc_stat(gwas_qc, "n_sig_snp"))

      col_labels <- c('GWAS Name',
                      'GWAS Label',
                      'N variants pre-QC',
                      'Identified genome build',
                      'N variants post-QC',
                      'Lambda GC',
                      'Max. chi^2',
                      'N genome-wide significant variants')

      # LDSC SNP-heritability and intercept have moved to the SNP-h² & rG
      # tab; they used to live here. Lambda GC / Max chi² / N sig SNPs
      # stay as QC metrics.

      names(qc_val) <- col_labels

      qc_val<-t(qc_val)
      qc_val<-data.table(Parameter=dimnames(qc_val)[[1]],
                         Value=qc_val[,1])

      qc_val
    })

    # Small, static 2-column summary - plain shiny::renderTable instead of a
    # DT widget, since this needs no sorting/searching/pagination (avoids the
    # DT/DataTables client-side rendering issues seen with this table).
    output$qc_table <- renderTable({
      dat <- qc_val()
      dat$Parameter <- paste0("<strong>", dat$Parameter, "</strong>")
      dat
    }, sanitize.text.function = function(x) x, colnames = FALSE, align = 'cc', rownames = FALSE)

    output$qc_legend <- renderUI({
      req(config_flags())
      gd_legend(list(
        "Lambda GC" = qc_help[["Lambda GC"]],
        "Max. chi^2" = qc_help[["Max. chi^2"]],
        "N genome-wide significant variants" = qc_help[["N genome-wide significant variants"]]
      ))
    })

    # MAF plot rendering
    output$maf_plot_ui <- renderUI({
      req(gwas_data(), selected_gwas())
      b64 <- gd_read(gwas_data(), selected_gwas(), "gwas_qc")$maf_plot_base64

      if (!is.null(b64)) {
        tagList(
          tags$img(
            src = paste0("data:image/png;base64,", b64),
            style = "max-height: 430px; width: auto;"
          ),
          tags$p(
            style = "font-size: 0.85em; color: #6c757d; margin-top: 8px; max-width: 800px;",
            "This is a quality-control check. It compares the allele frequency reported in your ",
            "GWAS summary statistics against the allele frequency of the same variant in the ",
            "reference panel (1000 Genomes). Large disagreements usually indicate allele ",
            "mis-coding, strand issues, an ancestry mismatch between your GWAS and the reference, ",
            "or data errors, so these variants are removed before downstream analysis."
          ),
          gd_legend(list(
            "X-axis" = "Reference-panel allele frequency.",
            "Y-axis" = "Allele frequency reported in your GWAS summary statistics.",
            "Diagonal line" = "y = x: points on this line agree between the GWAS and the reference.",
            "Points shown" = paste0(
              "Only variants whose two frequencies differ by more than 0.2 are plotted — these ",
              "are the discordant variants removed during QC. An empty or nearly-empty plot is a ",
              "good sign (few or no frequency mismatches).")
          ), heading = "How to read this plot")
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
      b64 <- gd_read(gwas_data(), selected_gwas(), "gwas_qc")$qq_plot_base64

      if (!is.null(b64)) {
        tagList(
          tags$img(
            src = paste0("data:image/png;base64,", b64),
            style = "max-height: 500px; width: auto;"
          ),
          gd_legend(list(
            "Diagonal line" = "Expected p-value distribution under the null hypothesis of no association.",
            "Points on the line" = "Test statistics are well calibrated (no inflation).",
            "Early upward departure" = paste0(
              "Points lifting above the line across most of the range suggests inflation or ",
              "confounding (compare with Lambda GC on the QC Summary tab)."),
            "Upward tail at the far right only" = paste0(
              "The strongest associations rising above the line, while the rest follows it, ",
              "is the expected signature of true genetic signal.")
          ), heading = "How to read this plot")
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
      log_lines <- gd_read(gwas_data(), selected_gwas(), "gwas_qc")$cleaner_dat$log

      if (is.null(log_lines) || length(log_lines) == 0) {
        return(p("Log not available."))
      }

      pre(
        style = "background-color: #1e1e1e; color: #00ff00; font-family: 'Courier New', monospace; font-size: 12px; padding: 15px; border-radius: 4px; white-space: pre-wrap;",
        paste(log_lines, collapse = "\n")
      )
    })

  })
}
