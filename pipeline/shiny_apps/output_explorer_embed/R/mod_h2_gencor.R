#' Build the bivariate-LDSC forest plot as a ggplot object.
#'
#' Shared by the on-screen renderPlot and the download handler so saved
#' files match what's on screen. Returns `NULL` if there's no data.
#'
#' @param tab data.table with columns label, rg, rg_se, rg_p, rg_p_fdr.
#' @param sort_choice one of "rg", "significance", "alphabetical".
#' @param font_size,point_size,theme_fn,title standard plot-options params.
# Columns that are always in the gencor results table (not offered as
# facet options).
GENCOR_CORE_COLS <- c("name", "label", "rg", "rg_se", "rg_p", "rg_p_fdr",
                      "n_snps", "gcov_int")

#' Detect facet-able extra columns in a gencor table (anything beyond the
#' core LDSC output columns).
detect_gencor_facet_cols <- function(tab) {
  if (is.null(tab) || nrow(tab) == 0) return(character(0))
  setdiff(names(tab), GENCOR_CORE_COLS)
}

build_gencor_plot <- function(tab, sort_choice = "rg", facet_by = NULL,
                               font_size = 13, point_size = 3,
                               theme_fn = ggplot2::theme_bw, title = "") {
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  dt <- data.table::as.data.table(tab)[!is.na(rg) & !is.na(rg_se)]
  if (nrow(dt) == 0) return(NULL)
  dt[, ci_lo := rg - 1.96 * rg_se]
  dt[, ci_hi := rg + 1.96 * rg_se]
  dt[, tier := data.table::fifelse(
        !is.na(rg_p_fdr) & rg_p_fdr < 0.05, "FDR significant",
        data.table::fifelse(!is.na(rg_p) & rg_p < 0.05, "Nominal (p < 0.05)",
                                                        "Non-significant"))]
  tier_levels <- c("FDR significant", "Nominal (p < 0.05)", "Non-significant")
  dt[, tier := factor(tier, levels = tier_levels)]

  # y-axis order — first factor level renders at the BOTTOM.
  ord <- switch(sort_choice %||% "rg",
    alphabetical = order(as.character(dt$label), decreasing = TRUE),
    significance = order(-log10(pmax(dt$rg_p, 1e-300)), na.last = FALSE),
    order(dt$rg))
  dt[, label := factor(label, levels = as.character(label)[ord])]

  # Optional facet column: bucket NA / empty into "Unknown".
  use_facet <- !is.null(facet_by) && nzchar(facet_by) && facet_by %in% names(dt)
  if (use_facet) {
    vals <- as.character(dt[[facet_by]])
    vals[is.na(vals) | vals == ""] <- "Unknown"
    dt[[facet_by]] <- vals
  }

  tier_cols <- c("FDR significant"    = "#d62728",
                 "Nominal (p < 0.05)" = "#ff7f0e",
                 "Non-significant"    = "#7f7f7f")
  tier_shapes <- c("FDR significant"    = 15,
                   "Nominal (p < 0.05)" = 17,
                   "Non-significant"    = 16)

  # Ghost rows so the legend renders a glyph for every tier even when the
  # data doesn't cover all three levels.
  ghost <- data.table::data.table(
    rg    = rep(dt$rg[1],    length(tier_levels)),
    ci_lo = rep(dt$rg[1],    length(tier_levels)),
    ci_hi = rep(dt$rg[1],    length(tier_levels)),
    label = rep(dt$label[1], length(tier_levels)),
    tier  = factor(tier_levels, levels = tier_levels)
  )
  if (use_facet) ghost[[facet_by]] <- dt[[facet_by]][1]

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = rg, y = label, colour = tier, shape = tier)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_point(data = ghost, size = point_size, alpha = 0) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lo, xmax = ci_hi), height = 0) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_manual(values = tier_cols, limits = tier_levels,
                                  drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = tier_shapes, limits = tier_levels,
                                 drop = FALSE, name = NULL) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = point_size + 1, alpha = 1)),
      shape  = ggplot2::guide_legend(override.aes = list(size = point_size + 1, alpha = 1))
    ) +
    ggplot2::labs(x = expression("Genetic correlation ("*r[g]*")"), y = NULL,
                  title = if (nzchar(title)) title else NULL) +
    theme_fn(base_size = font_size) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # facet_grid rows with space="free_y" so each band is sized by its trait
  # count, keeping row height constant. Strips render on the right (ggplot's
  # default for row facets).
  if (use_facet) {
    p <- p +
      ggplot2::facet_grid(rows = as.formula(paste0("`", facet_by, "` ~ .")),
                          scales = "free_y", space = "free_y") +
      ggplot2::theme(
        strip.background.y = ggplot2::element_rect(fill = "grey93", colour = NA),
        strip.text.y       = ggplot2::element_text(angle = 0, face = "bold")
      )
  }
  p
}

# SNP-h² and rG tab
#
# A dedicated top-level tab that groups the LDSC heritability estimates and
# bivariate LDSC (genetic correlation) results. Previously both lived under
# GWAS QC — but they are analyses in their own right, not just QC metrics.
# The QC-relevant LDSC output (intercept, ratio) can be inferred from what's
# here plus what's on the GWAS QC tab (Lambda GC, Max chi²).
#
# Two sub-tabs:
#   - "SNP heritability" — DT table with h², h²_SE, LDSC intercept, one row
#      per GWAS. Same table for single- and multi-GWAS mode (single = one row).
#   - "Genetic correlation (rG)" — single-GWAS forest plot (moved verbatim
#      from mod_gwas_qc.R) or the cross-GWAS gencor matrix (compare mode).

h2GencorUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "SNP-h² & rG",
    br(),
    p("SNP heritability (LDSC) and cross-trait genetic correlations (bivariate LDSC)."),
    hr(),
    tabsetPanel(
      id = ns("h2_gencor_tabs"),
      tabPanel(
        title = "SNP heritability",
        value = "h2",
        br(),
        div(style = "max-width: 900px;",
          DT::DTOutput(ns("h2_table"))
        ),
        br(),
        uiOutput(ns("h2_legend"))
      ),
      tabPanel(
        title = "Genetic correlation (rG)",
        value = "gencor",
        br(),
        uiOutput(ns("gencor_ui"))
      )
    )
  )
}

h2GencorServer <- function(id, gwas_data, selected_gwas, gwas_list, config_flags,
                            selected_gwas_multi = NULL,
                            comparison_mode = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Cross-GWAS gencor sub-module (registered unconditionally; self-guards
    # on selected_gwas_multi being length >= 1).
    if (!is.null(selected_gwas_multi)) {
      gencor_compare_server("gencor_compare",
                              gwas_data, selected_gwas_multi)
    }

    # -------------------------------------------------------------------
    # SNP heritability table
    # -------------------------------------------------------------------

    h2_gwas_vec <- reactive({
      req(gwas_data())
      if (!is.null(comparison_mode) && isTRUE(comparison_mode())) {
        selected_gwas_multi()
      } else {
        selected_gwas()
      }
    })

    h2_data <- reactive({
      req(gwas_data(), config_flags())
      build_heritability_long(gwas_data(), h2_gwas_vec())
    })

    output$h2_table <- DT::renderDT({
      req(config_flags())
      cf <- config_flags()
      if (!isTRUE(cf$ldsc)) {
        return(NULL)
      }
      d <- h2_data()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      out <- data.frame(
        GWAS                     = d$gwas,
        `SNP-h² (SE)`            = ifelse(is.na(d$h2), "—",
                                            sprintf("%.3f (%.3f)", d$h2, d$h2_se)),
        `LDSC intercept (SE)`    = ifelse(is.na(d$int_est), "—",
                                            sprintf("%.3f (%.3f)", d$int_est, d$int_se)),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(
        out, rownames = FALSE, escape = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE)
      )
    })

    output$h2_legend <- renderUI({
      req(config_flags())
      cf <- config_flags()
      if (!isTRUE(cf$ldsc)) {
        return(div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          tags$strong("LDSC not run"),
          br(),
          "Enable ldsc in the pipeline config to populate this view."
        ))
      }
      gd_legend(list(
        "SNP-h² (SE; observed scale)" = paste0(
          "LDSC estimate of the proportion of trait variance explained by ",
          "common SNPs (observed scale)."),
        "LDSC intercept (SE)"         = paste0(
          "Attenuation of the LDSC regression at the y-intercept. Values ",
          "close to 1 imply the test statistics are inflated by polygenic ",
          "signal rather than confounding; values noticeably above 1 ",
          "suggest population stratification / cryptic relatedness ",
          "(contrast with Lambda GC on the GWAS QC tab).")
      ))
    })

    # -------------------------------------------------------------------
    # Bivariate LDSC (genetic correlation)
    #
    # Single-GWAS: forest plot of rG between the selected GWAS and each
    # secondary trait tested. Compare mode: cross-GWAS rG matrix via the
    # gencor_compare_ui swap. Below reactives / renderers are lifted from
    # the previous location in mod_gwas_qc.R.
    # -------------------------------------------------------------------

    gencor_table <- reactive({
      req(gwas_data(), selected_gwas())
      d <- gd_read(gwas_data(), selected_gwas(), "gwas_qc")$ldsc_gencor_dat
      if (is.null(d) || is.null(d$table)) return(NULL)
      d$table
    })

    gencor_table_filtered <- reactive({
      tab <- gencor_table()
      if (is.null(tab)) return(NULL)
      picked <- input$gencor_traits
      if (!is.null(picked) && length(picked) > 0) {
        tab <- tab[as.character(tab$label) %in% picked, ]
      }
      if (isTRUE(as.logical(input$gencor_conf_only))) {
        tab <- tab[!is.na(tab$rg_p_fdr) & tab$rg_p_fdr < 0.05, ]
      }
      tab
    })

    gencor_table_filtered_plot <- gencor_table_filtered %>% debounce(400)

    plot_dim_gencor <- reactive({
      tab <- gencor_table_filtered_plot()
      if (is.null(tab) || nrow(tab) == 0) return(list(width_px = 700L, height_px = 220L))
      fs <- input$gencor_font_size %||% 13
      label_pt <- max(strwidth_pt(as.character(tab$label), ps = fs)) + 10
      label_px <- label_pt * 96 / 72
      chrome_px <- 60 + 480 + 40
      fb <- input$gencor_facet
      has_facet <- !is.null(fb) && nzchar(fb) && fb %in% names(tab)
      if (has_facet) {
        facet_vals <- as.character(tab[[fb]])
        facet_vals[is.na(facet_vals) | facet_vals == ""] <- "Unknown"
        strip_pt   <- max(strwidth_pt(facet_vals, ps = fs)) + 20
        chrome_px  <- chrome_px + strip_pt * 96 / 72
      }
      total_px <- round(label_px + chrome_px)
      width_px <- as.integer(min(max(total_px, 700L), 1800L))

      n <- sum(!is.na(tab$rg))
      height_px <- as.integer(max(220, 40 + 22 * n))
      if (has_facet) {
        n_facets  <- length(unique(tab[[fb]]))
        height_px <- height_px + as.integer(30 * max(n_facets - 1, 0))
      }
      list(width_px = width_px, height_px = height_px)
    })

    output$gencor_ui <- renderUI({
      req(gwas_data(), selected_gwas(), config_flags())

      in_compare <- !is.null(comparison_mode) && isTRUE(comparison_mode())
      if (in_compare) {
        return(gencor_compare_ui(NS(ns("gencor_compare"))))
      }

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

      all_labels <- sort(unique(as.character(tab$label)))
      facet_choices <- c("None" = "",
                         stats::setNames(detect_gencor_facet_cols(tab),
                                         detect_gencor_facet_cols(tab)))

      tagList(
        tags$p(
          style = "font-size: 0.85em; color: var(--gd-text-mute);",
          "Note: rg estimated from externally-munged sumstats can be affected by ",
          "strand-ambiguous SNPs and by low SNP overlap (see N SNPs)."
        ),
        tags$details(class = "gd-details",
          tags$summary("Filter data"),
          tags$div(class = "gd-details-body",
            tags$p(class = "gd-details-intro",
              "Restrict which secondary traits appear in the plot and how they ",
              "are ordered."),
            fluidRow(
              column(4,
                selectInput(ns("gencor_traits"), "Include these traits:",
                            choices = all_labels, selected = all_labels,
                            multiple = TRUE),
                selectInput(ns("gencor_facet"), "Facet by:",
                            choices = facet_choices, selected = "")
              ),
              column(4,
                selectInput(ns("gencor_sort"), "Sort traits by:",
                            choices = c("Genetic correlation (rg)" = "rg",
                                        "Significance" = "significance",
                                        "Alphabetical" = "alphabetical"),
                            selected = "rg")
              ),
              column(4,
                radioButtons(ns("gencor_conf_only"),
                             "Show FDR-significant traits only :",
                             choices = c("True" = T, "False" = F),
                             selected = F)
              )
            )
          )
        ),
        tags$details(class = "gd-details",
          tags$summary("Plot options"),
          tags$div(class = "gd-details-body",
            tags$p(class = "gd-details-intro",
              "Customise how the plot looks (title, theme, font size, point size) ",
              "and download it as a PNG, PDF, or SVG at the size and resolution you choose."),
            fluidRow(
              column(4,
                textInput(ns("gencor_title"), "Plot title (optional):", value = ""),
                selectInput(ns("gencor_theme"), "Theme:",
                            choices = c("Black & white" = "bw", "Minimal" = "minimal",
                                        "Classic" = "classic", "Light" = "light"),
                            selected = "bw")
              ),
              column(4,
                sliderInput(ns("gencor_font_size"), "Font size (pt):",
                            min = 8, max = 20, value = 13, step = 1),
                sliderInput(ns("gencor_point_size"), "Point size:",
                            min = 1, max = 8, value = 3, step = 1)
              ),
              column(4,
                selectInput(ns("gencor_dl_format"), "Download format:",
                            choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                            selected = "png"),
                numericInput(ns("gencor_dl_width"), "Width (inches):",
                             value = 8, min = 2, max = 40, step = 0.5),
                numericInput(ns("gencor_dl_height"), "Height (inches):",
                             value = 9, min = 2, max = 40, step = 0.5),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'png'", ns("gencor_dl_format")),
                  numericInput(ns("gencor_dl_dpi"), "Resolution (DPI, PNG only):",
                               value = 300, min = 72, max = 600, step = 25)
                ),
                downloadButton(ns("gencor_download"), "Download plot")
              )
            )
          )
        ),
        br(),
        uiOutput(ns("gencor_plot_ui")),
        gd_legend(list(
          "rg" = "Genetic correlation between this GWAS and each secondary trait, ranging from -1 to 1.",
          "Whiskers" = "95% confidence interval (rg plus/minus 1.96 x SE).",
          "Dashed vertical line" = "rg = 0 (no genetic correlation).",
          "Colour and shape" = paste0(
            "FDR-significant (red square), nominally significant p < 0.05 (orange triangle), ",
            "or non-significant (grey circle).")
        ), heading = "How to read this plot"),
        br(),
        fluidRow(column(width = 8, dataTableOutput(ns("gencor_table"))))
      )
    })

    output$gencor_table <- renderDataTable({
      tab <- gencor_table_filtered()
      req(tab)
      show <- tab[, .(label, rg, rg_se, rg_p, rg_p_fdr, n_snps)]
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

    output$gencor_plot_ui <- renderUI({
      plotOutput(ns("gencor_forest"), height = "auto",
                 width = paste0(plot_dim_gencor()$width_px, "px"))
    })

    output$gencor_forest <- renderPlot({
      build_gencor_plot(
        tab = gencor_table_filtered_plot(),
        sort_choice = input$gencor_sort %||% "rg",
        facet_by = input$gencor_facet %||% "",
        font_size = input$gencor_font_size %||% 13,
        point_size = input$gencor_point_size %||% 3,
        theme_fn = mol_theme_fn(input$gencor_theme),
        title = input$gencor_title %||% ""
      )
    }, res = 96, height = function() plot_dim_gencor()$height_px)

    observeEvent(plot_dim_gencor(), {
      dims <- plot_dim_gencor()
      updateNumericInput(session, "gencor_dl_width",
                         value = round(dims$width_px  / 96, 1))
      updateNumericInput(session, "gencor_dl_height",
                         value = round(dims$height_px / 96, 1))
    })

    output$gencor_download <- downloadHandler(
      filename = function() {
        sprintf("bivariate_ldsc_%s.%s",
                format(Sys.time(), "%Y%m%d_%H%M%S"),
                input$gencor_dl_format %||% "png")
      },
      content = function(file) {
        p <- build_gencor_plot(
          tab = gencor_table_filtered(),
          sort_choice = input$gencor_sort %||% "rg",
          font_size = input$gencor_font_size %||% 13,
          point_size = input$gencor_point_size %||% 3,
          theme_fn = mol_theme_fn(input$gencor_theme),
          title = input$gencor_title %||% ""
        )
        if (is.null(p)) {
          grDevices::png(file, width = 4, height = 1, units = "in", res = 96)
          grid::grid.text("No genetic-correlation results to plot.")
          grDevices::dev.off()
          return(invisible())
        }
        w <- input$gencor_dl_width  %||% 8
        h <- input$gencor_dl_height %||% 9
        fmt <- input$gencor_dl_format %||% "png"
        switch(fmt,
          png = grDevices::png(file, width = w, height = h, units = "in",
                               res = input$gencor_dl_dpi %||% 300),
          pdf = grDevices::pdf(file, width = w, height = h),
          svg = grDevices::svg(file, width = w, height = h)
        )
        print(p)
        grDevices::dev.off()
      }
    )

    # Auto-jump between sub-tabs on compare-mode transition would confuse
    # users here — both sub-tabs work in both modes (h² shows 1 vs N rows;
    # rG swaps content in-place), so no showTab/hideTab / auto-focus is
    # needed.
  })
}
