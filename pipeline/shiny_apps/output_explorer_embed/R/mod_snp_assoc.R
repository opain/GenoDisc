#' SNP Associations module UI
#'
#' @param id Module namespace id
#' @return UI elements for the SNP Associations tab
snpAssocUI <- function(id) {
  ns <- NS(id)

  tabPanel(
    title = "SNP Associations",
    br(),
    p("This tab shows SNP association results. Select the Lead variant tab below to view information for independent lead variants identified by either LD-based clumping or COJO. Select the Fine-mapping tab below to view SuSiE Finemaping results."),
    hr(),
    tabsetPanel(
      id = ns("snp_assoc_tabs"),
      tabPanel(
        title = "Manhattan plot",
        value = "manhattan",
        br(),
        fluidPage(
          sidebarPanel(
            radioButtons(ns("manhattan_label_choice"), "Variant labelling:",
                         choices = c("Without labels" = "unlabelled",
                                     "With nearest-gene labels" = "labelled"),
                         selected = "unlabelled"),
            width = 3
          ),

          mainPanel(
            uiOutput(ns("manhattan_plot_ui")),
            width = 9
          )
        )
      ),
      tabPanel(
        title = "Lead variants",
        value = "lead_variants",
        br(),
        fluidPage(
          sidebarPanel(
            radioButtons(ns("clumping_type"), "Select method:",
                         choices = c("COJO" = "cojo_analysis",
                                     "LD-based clumping" = "ld_clumping"),
                         selected = "cojo_analysis"),
            radioButtons(ns("pvalue_threshold"), "Select P-value Threshold:",
                         choices = c("Genome-wide significance (p < 5e-8)" = 5e-8),
                         selected = 5e-8),
            width = 3
          ),

          mainPanel(
            uiOutput(ns("cojo_status_message")),
            dataTableOutput(ns("snp_assoc_lead_table")),
            br(),
            uiOutput(ns("locus_plot_ui")),
            width = 9
          )
        )
      ),
      tabPanel(
        title = "Fine-mapping",
        value = "finemapping",
        br(),
        fluidPage(
          sidebarPanel(
            radioButtons(ns("l_param"), "Select L parameter:",
                         choices = c("L1" = "L1",
                                     "L10" = "L10"),
                         selected = "L1"),
            width = 3
          ),

          mainPanel(
            dataTableOutput(ns("snp_assoc_finemap_table")),
            width = 6
          )
        )
      )
    )
  )
}

#' SNP Associations module server
#'
#' @param id Module namespace id
#' @param gwas_data Reactive returning the loaded GWAS data list
#' @param selected_gwas Reactive returning the currently selected GWAS name
snpAssocServer <- function(id, gwas_data, selected_gwas, config_flags) {
  moduleServer(id, function(input, output, session) {

    observeEvent(config_flags(), {
      cf <- config_flags()

      if (any(cf$clump, cf$cojo)) {
        showTab("snp_assoc_tabs", "lead_variants", session = session)
      } else {
        hideTab("snp_assoc_tabs", "lead_variants", session = session)
      }

      if (isTRUE(cf$finemap)) {
        showTab("snp_assoc_tabs", "finemapping", session = session)
      } else {
        hideTab("snp_assoc_tabs", "finemapping", session = session)
      }

      lead_choices <- c()
      if (isTRUE(cf$cojo)) lead_choices <- c(lead_choices, "COJO" = "cojo_analysis")
      if (isTRUE(cf$clump)) lead_choices <- c(lead_choices, "LD-based clumping" = "ld_clumping")
      if (length(lead_choices) > 0) {
        updateRadioButtons(session, "clumping_type",
                           choices = lead_choices,
                           selected = unname(lead_choices[1]))
      }
    })

    # COJO output is already thresholded at 5e-8, so the suggestive (1e-5) filter is only
    # meaningful for LD-based clumping (which runs at --clump-p1 1e-5) - offer it there only.
    observeEvent(input$clumping_type, {
      if (input$clumping_type == "ld_clumping") {
        updateRadioButtons(session, "pvalue_threshold",
                           choices = c("Genome-wide significance (p < 5e-8)" = 5e-8,
                                       "Suggestive significance (p < 1e-5)" = 1e-5),
                           selected = 5e-8)
      } else {
        updateRadioButtons(session, "pvalue_threshold",
                           choices = c("Genome-wide significance (p < 5e-8)" = 5e-8),
                           selected = 5e-8)
      }
    })

    output$manhattan_plot_ui <- renderUI({
      req(gwas_data(), selected_gwas(), input$manhattan_label_choice)
      b64_field <- if (input$manhattan_label_choice == "labelled")
        "manhattan_plot_labelled_base64"
      else
        "manhattan_plot_unlabelled_base64"
      b64 <- gd_read(gwas_data(), selected_gwas(), "gwas_qc")[[b64_field]]

      if (!is.null(b64)) {
        tags$img(
          src = paste0("data:image/png;base64,", b64),
          style = "max-width: 100%; height: auto;"
        )
      } else {
        div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          icon("chart-bar", style = "font-size: 2em;"),
          br(), br(),
          tags$strong("Manhattan plot Unavailable"),
          br(),
          "This results package was produced before the Manhattan plot rule was added."
        )
      }
    })

    snp_assoc_lead_data <- reactive({
      req(gwas_data(), selected_gwas(), input$clumping_type)
      snp_assoc <- gd_read(gwas_data(), selected_gwas(), "snp_assoc")
      snp_assoc_lead <- NULL
      if (input$clumping_type == "ld_clumping") {
        snp_assoc_lead <- snp_assoc$clump
      }

      if (input$clumping_type == "cojo_analysis") {
        snp_assoc_lead <- snp_assoc$cojo
      }

      req(snp_assoc_lead)

      snp_assoc_lead <- snp_assoc_lead[, names(snp_assoc_lead) %in% c("CHR","BP","SNP","A1","A2","BETA","SE","P","NearestGene"), with = F]

      snp_assoc_lead$P <- as.numeric(snp_assoc_lead$P)

      snp_assoc_lead <- snp_assoc_lead[snp_assoc_lead$P < as.numeric(input$pvalue_threshold), ]

      snp_assoc_lead$BETA <- round(snp_assoc_lead$BETA, 2)
      snp_assoc_lead$SE <- round(snp_assoc_lead$SE, 2)

      snp_assoc_lead <- snp_assoc_lead[order(snp_assoc_lead$CHR, snp_assoc_lead$BP), ]

      return(snp_assoc_lead)
    })

    # COJO can fail per chromosome when the number of independent genome-wide-significant
    # signals exceeds the LD reference sample size (~503, 1000 Genomes EUR). package_results
    # carries snp_assoc$cojo_status so we can flag partial results, or a full failure.
    output$cojo_status_message <- renderUI({
      req(gwas_data(), selected_gwas(), input$clumping_type)
      if (input$clumping_type != "cojo_analysis") return(NULL)

      cojo_status <- gd_read(gwas_data(), selected_gwas(), "snp_assoc")$cojo_status
      if (is.null(cojo_status) || !isTRUE(cojo_status$any_failed)) return(NULL)

      failed_txt <- paste(cojo_status$failed_chrs, collapse = ", ")

      if (isTRUE(cojo_status$all_failed)) {
        div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 30px 20px; text-align: center; color: #6c757d;",
          icon("exclamation-triangle", style = "font-size: 1.5em;"),
          br(), br(),
          tags$strong("COJO could not be completed"),
          br(),
          paste0("The number of independent genome-wide-significant signals exceeded the LD ",
                 "reference panel size (~503 individuals, 1000 Genomes EUR) on every chromosome ",
                 "attempted (", failed_txt, "), so no joint model could be fitted. LD-based ",
                 "clumping is available as an alternative.")
        )
      } else {
        div(
          style = "background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 8px; padding: 15px 20px; color: #856404; margin-bottom: 15px;",
          icon("exclamation-triangle"),
          tags$strong(" Partial COJO results"),
          br(),
          paste0("Chromosome(s) ", failed_txt, " could not be analysed because the number of ",
                 "independent genome-wide-significant signals exceeded the LD reference panel size ",
                 "(~503 individuals, 1000 Genomes EUR). The table below shows COJO results for the ",
                 "chromosomes that completed successfully.")
        )
      }
    })

    output$snp_assoc_lead_table <- renderDataTable({
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "}"
      )

      datatable(snp_assoc_lead_data(), rownames = F, width = 7,
        selection = list(mode = 'single', target = 'row'),
        options = list(
          rowCallback = JS(js),
          columnDefs = list(list(className = 'dt-center', targets = 0:7),
                            list(width = '60px', targets = 7),
                            list(width = '600px', targets = 8))))
    })

    output$locus_plot_ui <- renderUI({
      req(gwas_data(), selected_gwas(), input$clumping_type)

      if (input$clumping_type != "ld_clumping") {
        return(div(
          style = "padding: 20px; text-align: center; color: #6c757d; font-style: italic;",
          "Locus plots are available for LD-based clumping results only."
        ))
      }

      locus_plots <- gd_read(gwas_data(), selected_gwas(), "snp_assoc")$locus_plots

      if (is.null(locus_plots) || length(locus_plots) == 0) {
        return(div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 30px 20px; text-align: center; color: #6c757d;",
          icon("chart-area", style = "font-size: 1.5em;"),
          br(), br(),
          tags$strong("No locus plots available"),
          br(),
          "Locus plots are generated for index variants with p < 1e-5 (top 50 by p-value)."
        ))
      }

      selected_row <- input$snp_assoc_lead_table_rows_selected
      if (is.null(selected_row) || length(selected_row) == 0) {
        return(div(
          style = "padding: 20px; text-align: center; color: #6c757d; font-style: italic;",
          "Select a variant from the table above to view its locus plot."
        ))
      }

      tbl <- snp_assoc_lead_data()
      if (selected_row > nrow(tbl)) return(NULL)
      snp_id <- as.character(tbl[selected_row, "SNP"][[1]])

      b64 <- locus_plots[[snp_id]]
      if (is.null(b64)) {
        return(div(
          style = "padding: 20px; text-align: center; color: #6c757d; font-style: italic;",
          paste0("No locus plot available for ", snp_id,
                 " (variant may be outside the top 50 loci by p-value).")
        ))
      }

      tags$img(
        src = paste0("data:image/png;base64,", b64),
        style = "max-width: 100%; height: auto;"
      )
    })

    snp_assoc_finemap_data <- reactive({
      req(gwas_data(), selected_gwas(), input$l_param)
      snp_assoc <- gd_read(gwas_data(), selected_gwas(), "snp_assoc")
      snp_assoc_finemap <- NULL
      if (input$l_param == "L1") {
        snp_assoc_finemap <- snp_assoc$susie$L1
      }

      if (input$l_param == "L10") {
        snp_assoc_finemap <- snp_assoc$susie$L10
      }

      req(snp_assoc_finemap)

      snp_assoc_finemap$cs_log10bf <- NULL
      snp_assoc_finemap$cs_avg_r2 <- NULL
      snp_assoc_finemap$cs_min_r2 <- NULL
      snp_assoc_finemap$TopPIP <- round(snp_assoc_finemap$TopPIP, 2)

      return(snp_assoc_finemap)
    })

    output$snp_assoc_finemap_table <- renderDataTable({
      datatable(snp_assoc_finemap_data(), rownames = F, width = 7, options = list(
        columnDefs = list(list(className = 'dt-center', targets = 0:5))))
    })
  })
}
