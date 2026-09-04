
gwasQcUI <- function(id) {
  ns <- NS(id)
  # Small helper for the per-tab GWAS picker. Populated server-side from
  # selected_gwas_multi(); default = the first-selected GWAS. Rendered
  # inside each per-GWAS tab so the user can flick between GWAS without
  # leaving the tab.
  gwas_picker <- function(input_id, label = "GWAS:") {
    div(style = "max-width: 260px; margin-bottom: 10px;",
      selectInput(ns(input_id), label, choices = NULL, multiple = FALSE)
    )
  }
  tabPanel(
    title = "GWAS QC",
    br(),
    p("This tab shows key quality-control statistics for the loaded GWAS. ",
      "The QC Summary table covers every selected GWAS; the other sub-tabs ",
      "are per-GWAS and expose a GWAS selector at the top so you can flick ",
      "between them without changing your bundle-wide selection."),
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
        gwas_picker("maf_gwas"),
        uiOutput(ns("maf_plot_ui"))
      ),
      tabPanel(
        title = "QQ Plot",
        value = "qq_plot",
        br(),
        gwas_picker("qq_gwas"),
        uiOutput(ns("qq_plot_ui"))
      ),
      # SNP-h² and rG have moved to their own top-level tab
      # (SNP-h² & rG). Removing them keeps this tab focused on QC.
      tabPanel(
        title = "Sumstat QC Log",
        value = "cleaner_log",
        br(),
        gwas_picker("log_gwas"),
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

    # Populate the per-tab GWAS selectors from selected_gwas_multi().
    # Default = the first-selected GWAS. Preserves the user's per-tab
    # pick across bundle-selection changes when the pick is still valid.
    .populate_per_gwas_picker <- function(input_id) {
      observe({
        choices <- selected_gwas_multi()
        req(length(choices) >= 1)
        cur <- isolate(input[[input_id]])
        keep <- if (!is.null(cur) && cur %in% choices) cur else choices[1L]
        updateSelectInput(session, input_id, choices = choices, selected = keep)
      })
    }
    if (!is.null(selected_gwas_multi)) {
      .populate_per_gwas_picker("maf_gwas")
      .populate_per_gwas_picker("qq_gwas")
      .populate_per_gwas_picker("log_gwas")
    }

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
        "Zero is common for underpowered GWAS."),
      "LDSC intercept (SE)" = paste0(
        "LD-score regression intercept. About 1.0 indicates little confounding or ",
        "sample overlap; values above 1 suggest residual confounding or population ",
        "stratification rather than polygenic signal (contrast with Lambda GC). ",
        "Standard error in brackets.")
    )

    # Which GWAS to summarise: the scalar selected GWAS in single mode,
    # the full multi-select vector in compare mode.
    qc_gwas_vec <- reactive({
      req(gwas_data())
      if (!is.null(comparison_mode) && isTRUE(comparison_mode())) {
        selected_gwas_multi()
      } else {
        selected_gwas()
      }
    })

    # Parameter rows shown in the QC Summary table, in display order.
    # (LDSC SNP-h² has moved to the SNP-h² & rG tab; the LDSC intercept
    # stays here as a QC-relevant metric — values close to 1 imply
    # polygenic signal, > 1 suggests confounding.)
    .qc_row_labels <- c(
      "GWAS Label"                         = "label",
      "N variants pre-QC"                  = "n_var_orig",
      "Identified genome build"            = "build",
      "N variants post-QC"                 = "n_snp_final",
      "Lambda GC"                          = "lambda_gc",
      "Max. chi^2"                         = "max_chi2",
      "N genome-wide significant variants" = "n_sig_snp",
      "LDSC intercept (SE)"                = "int_paren"
    )

    .fmt_qc_cell <- function(x, kind) {
      if (all(is.na(x))) return("—")
      switch(kind,
        n_var_orig  = ,
        n_snp_final = ,
        n_sig_snp   = ifelse(is.na(x), "—", formatC(as.integer(x), format = "d", big.mark = ",")),
        lambda_gc   = ,
        max_chi2    = ifelse(is.na(x), "—", sprintf("%.3f", as.numeric(x))),
        # label / build / int_paren / anything else: as-is with NA -> em-dash
        ifelse(is.na(x) | !nzchar(as.character(x)), "—", as.character(x))
      )
    }

    # Wide table (single row per GWAS from the data-layer helper) reshaped
    # into the display shape: rows are parameters, columns are GWAS.
    qc_val <- reactive({
      req(gwas_data(), config_flags())
      req(length(qc_gwas_vec()) >= 1)
      cf <- config_flags()

      long <- build_qc_summary_long(gd = gwas_data(),
                                      gwas_vec = qc_gwas_vec(),
                                      gwas_list = gwas_list())
      if (nrow(long) == 0) return(NULL)

      long <- long[match(qc_gwas_vec(), gwas)]

      # LDSC intercept lives on the QC tab as a QC-relevant stat
      # (values > 1 suggest confounding vs. polygenic signal). Compose
      # the "estimate (SE)" string once so the display renders directly.
      long[, int_paren := ifelse(is.na(int_est), NA_character_,
                                   sprintf("%.3f (%.3f)", int_est, int_se))]

      # Drop the intercept row if LDSC wasn't run — no need to show a
      # column of em-dashes.
      row_labels <- .qc_row_labels
      if (!isTRUE(cf$ldsc)) {
        row_labels <- row_labels[row_labels != "int_paren"]
      }

      cols <- lapply(seq_len(nrow(long)), function(i) {
        vapply(names(row_labels), function(lab) {
          field <- row_labels[[lab]]
          .fmt_qc_cell(long[[field]][i], field)
        }, character(1))
      })
      out <- data.frame(
        Parameter = names(row_labels),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      for (i in seq_along(cols)) out[[long$gwas[i]]] <- cols[[i]]
      out
    })

    # Small static summary — plain shiny::renderTable (avoids DT client-side
    # rendering quirks seen with this table). Parameter column is bold.
    output$qc_table <- renderTable({
      dat <- qc_val()
      req(dat)
      dat$Parameter <- paste0("<strong>", dat$Parameter, "</strong>")
      # colnames only meaningful in compare mode (single mode = "GWAS_NAME"
      # column header is redundant with the label row). Show colnames in
      # both modes for consistency.
      dat
    }, sanitize.text.function = function(x) x,
       colnames = TRUE,
       align = "l",
       rownames = FALSE)

    output$qc_legend <- renderUI({
      req(config_flags())
      cf <- config_flags()
      items <- list(
        "Lambda GC" = qc_help[["Lambda GC"]],
        "Max. chi^2" = qc_help[["Max. chi^2"]],
        "N genome-wide significant variants" = qc_help[["N genome-wide significant variants"]]
      )
      if (isTRUE(cf$ldsc)) {
        items[["LDSC intercept (SE)"]] <- qc_help[["LDSC intercept (SE)"]]
      }
      gd_legend(items)
    })

    # Falls back to selected_gwas() when the per-tab picker hasn't been
    # populated yet (e.g. first render before the observe fires).
    .pick_gwas <- function(input_id) {
      v <- input[[input_id]]
      if (is.null(v) || !nzchar(v)) selected_gwas() else v
    }

    # MAF plot rendering
    output$maf_plot_ui <- renderUI({
      req(gwas_data())
      g <- .pick_gwas("maf_gwas")
      req(g)
      b64 <- gd_read(gwas_data(), g, "gwas_qc")$maf_plot_base64

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
      req(gwas_data())
      g <- .pick_gwas("qq_gwas")
      req(g)
      b64 <- gd_read(gwas_data(), g, "gwas_qc")$qq_plot_base64

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
      req(gwas_data())
      g <- .pick_gwas("log_gwas")
      req(g)
      log_lines <- gd_read(gwas_data(), g, "gwas_qc")$cleaner_dat$log

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
