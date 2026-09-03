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
        div(style = "max-width: 700px;", tableOutput(ns("qc_table"))),
        uiOutput(ns("qc_legend"))
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

gwasQcServer <- function(id, gwas_data, selected_gwas, gwas_list, config_flags,
                          selected_gwas_multi = NULL,
                          comparison_mode = NULL) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # Cross-GWAS gencor sub-module. Registered unconditionally; self-guards
    # on selected_gwas_multi being length >= 1.
    if (!is.null(selected_gwas_multi)) {
      gencor_compare_server("gencor_compare",
                              gwas_data, selected_gwas_multi)
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
      "LDSC SNP-heritability (SE; observed scale)" = paste0(
        "Proportion of trait variance explained by common SNPs, estimated by LD-score ",
        "regression, on the observed scale. Standard error in brackets."),
      "LDSC intercept (SE)" = paste0(
        "LD-score regression intercept. About 1.0 indicates little confounding or sample ",
        "overlap; values above 1 suggest residual confounding rather than polygenic signal ",
        "(contrast with Lambda GC). Standard error in brackets.")
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

      if (isTRUE(cf$ldsc)) {
        qc_val[, obs_h2 := paste0(round(gwas_qc$ldsc_dat$val$obs_h2_est,3), " (",round(gwas_qc$ldsc_dat$val$obs_h2_se,3),")")]
        qc_val[, int := paste0(round(gwas_qc$ldsc_dat$val$int_est,3), " (",round(gwas_qc$ldsc_dat$val$int_se,3),")")]
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

    # Always-visible legend under the QC table; LDSC rows only when LDSC was run.
    output$qc_legend <- renderUI({
      req(config_flags())
      cf <- config_flags()
      items <- list(
        "Lambda GC" = qc_help[["Lambda GC"]],
        "Max. chi^2" = qc_help[["Max. chi^2"]],
        "N genome-wide significant variants" = qc_help[["N genome-wide significant variants"]]
      )
      if (isTRUE(cf$ldsc)) {
        items[["LDSC SNP-heritability"]] <- qc_help[["LDSC SNP-heritability (SE; observed scale)"]]
        items[["LDSC intercept"]] <- qc_help[["LDSC intercept (SE)"]]
      }
      gd_legend(items)
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

    # Bivariate LDSC (genetic correlation) rendering
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

    # Debounced version used by the (expensive) plot render + dim/download
    # sync so quick successive filter tweaks (e.g. deselecting many traits
    # in the selectize) collapse into one plot rebuild. The table renderer
    # and download handler stay on the immediate reactive.
    gencor_table_filtered_plot <- gencor_table_filtered %>% debounce(400)

    # Plot width adapts to the longest visible trait label — long labels
    # (e.g. "Consumption (AUDIT-C)") were being cut off by the 900px cap.
    # When faceting, add a small extra allowance for the left-side facet strip.
    plot_dim_gencor <- reactive({
      tab <- gencor_table_filtered_plot()
      if (is.null(tab) || nrow(tab) == 0) return(list(width_px = 700L, height_px = 220L))
      fs <- input$gencor_font_size %||% 13
      label_pt <- max(strwidth_pt(as.character(tab$label), ps = fs)) + 10
      label_px <- label_pt * 96 / 72   # approx pt -> CSS px at 96dpi
      # Chrome: y-axis text padding + plot panel area + right margin.
      chrome_px <- 60 + 480 + 40
      fb <- input$gencor_facet
      has_facet <- !is.null(fb) && nzchar(fb) && fb %in% names(tab)
      if (has_facet) {
        # Left-side facet strip: measure widest facet value at the current font
        # so long category names ("Endocrine / metabolic") don't get clipped.
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

      # In compare mode, swap to the cross-GWAS gencor matrix view; the
      # single-GWAS forest-plot layout below is untouched.
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

    # Keep download width/height in sync with the on-screen plot so the
    # default download matches the current window. User overrides get
    # replaced whenever the on-screen dims change.
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
  })
}
