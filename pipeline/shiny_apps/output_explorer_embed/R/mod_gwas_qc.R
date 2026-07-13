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
        div(style = "max-width: 700px;", tableOutput(ns("qc_table")))
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
      ),
      tabPanel(
        title = "Bivariate LDSC",
        br(),
        uiOutput(ns("gencor_ui"))
      )
    )
  )
}

gwasQcServer <- function(id, gwas_data, selected_gwas, gwas_list, config_flags) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # Create a table showing key statistics
    qc_val <- reactive({
      req(gwas_data(), selected_gwas(), gwas_list(), config_flags())

      g <- gwas_data()[[selected_gwas()]]
      cf <- config_flags()

      # Build isn't always identified from CHR/BP (the cleaner can fall back to
      # matching by SNP ID instead, which never determines a build - see
      # extract_build() in package_results_functions.R). Coalesce to a scalar
      # so a single missing QC field can't collapse the whole table via
      # data.table()'s zero-length recycling.
      build_val <- g$gwas_qc$cleaner_dat$val$build$build
      if (length(build_val) == 0 || is.na(build_val)) build_val <- "Unknown"

      qc_val<-data.table(name=selected_gwas(),
                         label=gwas_list()$label[gwas_list()$name == selected_gwas()],
                         n_var_orig=g$gwas_qc$cleaner_dat$val$n_var_orig,
                         build=build_val,
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

    # Small, static 2-column summary - plain shiny::renderTable instead of a
    # DT widget, since this needs no sorting/searching/pagination (avoids the
    # DT/DataTables client-side rendering issues seen with this table).
    output$qc_table <- renderTable({
      dat <- qc_val()
      dat$Parameter <- paste0("<strong>", dat$Parameter, "</strong>")
      dat
    }, sanitize.text.function = function(x) x, colnames = FALSE, align = 'cc', rownames = FALSE)

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

    # Bivariate LDSC (genetic correlation) rendering
    gencor_table <- reactive({
      req(gwas_data(), selected_gwas())
      d <- gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_gencor_dat
      if (is.null(d) || is.null(d$table)) return(NULL)
      d$table
    })

    output$gencor_ui <- renderUI({
      req(gwas_data(), selected_gwas(), config_flags())
      cf  <- config_flags()
      tab <- gencor_table()

      if (!isTRUE(cf$ldsc_gencor) || is.null(tab) || nrow(tab) == 0) {
        return(div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          icon("link", style = "font-size: 2em;"),
          br(), br(),
          tags$strong("Bivariate LDSC not run"),
          br(),
          "Set gencor_gwas_list in the pipeline config to enable genetic-correlation analysis."
        ))
      }

      tagList(
        tags$p(
          style = "font-size: 0.85em; color: #6c757d;",
          "Note: rg estimated from externally-munged sumstats can be affected by ",
          "strand-ambiguous SNPs and by low SNP overlap (see N SNPs)."
        ),
        div(style = "max-width: 900px;", dataTableOutput(ns("gencor_table"))),
        br(),
        plotOutput(ns("gencor_forest"), height = "auto")
      )
    })

    output$gencor_table <- renderDataTable({
      tab <- gencor_table()
      req(tab)
      show <- tab[, .(label, rg, rg_se, rg_p, rg_p_fdr, n_snps)]
      # rowCallback replaces null cells (R's NA, after DT's JSON encoding) with
      # the literal string "NA". Runs after formatRound / formatSignif, so the
      # numeric formatting still applies to non-NA values.
      na_callback <- JS(
        "function(row, data, displayNum, displayIndex, dataIndex) {",
        "  for (var i = 0; i < data.length; i++) {",
        "    if (data[i] === null) { $('td:eq(' + i + ')', row).html('NA'); }",
        "  }",
        "}"
      )
      datatable(
        show, rownames = FALSE,
        colnames = c("Secondary trait", "rg", "SE", "p", "FDR p", "N SNPs"),
        options  = list(
          pageLength   = 25,
          order        = list(list(1, 'asc')),
          rowCallback  = na_callback
        )
      ) %>%
        formatRound(columns = c("rg", "rg_se"), digits = 3) %>%
        formatSignif(columns = c("rg_p", "rg_p_fdr"), digits = 3)
    })

    output$gencor_forest <- renderPlot({
      tab <- gencor_table()
      req(tab)
      dt <- tab[!is.na(rg) & !is.na(rg_se)]
      if (nrow(dt) == 0) return(NULL)
      dt[, ci_lo := rg - 1.96 * rg_se]
      dt[, ci_hi := rg + 1.96 * rg_se]
      dt[, tier := fifelse(!is.na(rg_p_fdr) & rg_p_fdr < 0.05, "FDR significant",
                    fifelse(!is.na(rg_p)     & rg_p     < 0.05, "Nominal (p < 0.05)",
                                                                "Non-significant"))]
      dt[, tier  := factor(tier, levels = c("FDR significant",
                                            "Nominal (p < 0.05)",
                                            "Non-significant"))]
      dt[, label := factor(label, levels = label[order(rg)])]

      tier_levels <- c("FDR significant",
                       "Nominal (p < 0.05)",
                       "Non-significant")
      tier_cols <- c("FDR significant"    = "#d62728",
                     "Nominal (p < 0.05)" = "#ff7f0e",
                     "Non-significant"    = "#7f7f7f")
      tier_shapes <- c("FDR significant"    = 15,  # filled square
                       "Nominal (p < 0.05)" = 17,  # filled triangle
                       "Non-significant"    = 16)  # filled circle

      # Ghost rows: one per tier, plotted with alpha = 0 so they're invisible
      # but make the colour + shape scales see every tier level. Without these,
      # `limits = tier_levels` shows legend text for empty tiers but no glyph,
      # because ggplot draws the legend key from the geom's data and tiers
      # with zero observations contribute no key glyph.
      ghost <- data.table(
        rg    = rep(dt$rg[1],    length(tier_levels)),
        ci_lo = rep(dt$rg[1],    length(tier_levels)),
        ci_hi = rep(dt$rg[1],    length(tier_levels)),
        label = rep(dt$label[1], length(tier_levels)),
        tier  = factor(tier_levels, levels = tier_levels)
      )

      ggplot(dt, aes(x = rg, y = label, colour = tier, shape = tier)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
        # Invisible ghost points seed all tier levels into the scales so the
        # legend renders a glyph for every tier, not just observed ones.
        geom_point(data = ghost, size = 3, alpha = 0) +
        # height = 0 removes the vertical end-cap whiskers, so the only lines
        # touching each point are the horizontal 95% CI segments.
        geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0) +
        geom_point(size = 3) +
        scale_colour_manual(values = tier_cols,   limits = tier_levels, drop = FALSE, name = NULL) +
        scale_shape_manual( values = tier_shapes, limits = tier_levels, drop = FALSE, name = NULL) +
        guides(
          colour = guide_legend(override.aes = list(size = 3.5, alpha = 1)),
          shape  = guide_legend(override.aes = list(size = 3.5, alpha = 1))
        ) +
        labs(x = expression("Genetic correlation ("*r[g]*")"), y = NULL) +
        cowplot::theme_cowplot() +
        theme(legend.position = "top")
    }, res = 96, height = function() {
      tab <- isolate(gencor_table())
      if (is.null(tab)) return(220)
      n <- sum(!is.na(tab$rg))
      max(220, 40 + 22 * n)
    })
  })
}
