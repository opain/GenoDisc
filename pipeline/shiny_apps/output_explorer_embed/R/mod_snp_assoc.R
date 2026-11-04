#' Parse a user-typed P-value string into a positive number, falling back
#' to `default` if the value is missing, non-numeric, or non-positive.
#' Lets the UI use textInput (which round-trips 5e-8 correctly) instead of
#' numericInput (which shows 5e-8 as 0.00000005).
parse_p <- function(x, default) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1 || is.na(v) || v <= 0) default else v
}

#' Build the Manhattan plot as a ggplot object.
#'
#' Shared by the on-screen renderPlot and the download handler so saved
#' files match what's on screen. Returns NULL if `manhattan_data` is missing.
#'
#' @param md List with elements `data` (filtered variants with cum_bp,
#'   neglogP, chr_colour, category), `offsets` (per-chromosome axis labels
#'   with cum_offset, label_pos), `indexes` (clumping leads with cum_bp,
#'   neglogP, label = nearest gene), `p_threshold` (effective filter P).
#' @param sig_threshold Numeric P for the dashed significance line.
#' @param highlight_threshold Numeric P below which index variants get the
#'   enlarged diamond glyph.
#' @param label_variants Whether to draw ggrepel labels on index leads.
#' @param label_threshold Numeric P — only label indexes below this.
#' @param font_size,point_size,theme_fn,title Standard plot-options params.
build_manhattan_plot <- function(md,
                                  sig_threshold = 5e-8,
                                  suggestive_line = FALSE,
                                  suggestive_threshold = 1e-5,
                                  highlight_threshold = 5e-8,
                                  label_variants = FALSE,
                                  label_threshold = 5e-8,
                                  protein_coding_labels = FALSE,
                                  font_size = 14,
                                  point_size = 0.9,
                                  theme_fn = ggplot2::theme_bw,
                                  title = "",
                                  col_chr_odd  = "#8c8c8c",
                                  col_chr_even = "#4d4d4d",
                                  col_member   = "#2ca02c",
                                  col_index    = "#006400") {
  if (is.null(md) || is.null(md$data) || is.null(md$offsets)) return(NULL)
  ss <- data.table::copy(md$data)
  offsets <- md$offsets
  idx <- md$indexes

  # Recompute alternating-chromosome colour from CHR so it reacts to the
  # user's colour choices — ignore any baked-in chr_colour from old bundles.
  ss[, chr_bg := ifelse(CHR %% 2 == 1, col_chr_odd, col_chr_even)]

  # Only highlight clump members whose parent index variant passes the
  # user's highlight threshold. Members of indexes above the threshold
  # fall back to the alternating chromosome grey (drawn by layer 1 below).
  # Older bundles that predate the `index_snp` column highlight all members
  # (backwards-compatible fallback).
  hi_index_snps <- if (!is.null(idx) && nrow(idx) > 0) {
    idx[P < highlight_threshold, SNP]
  } else character(0)
  ss_members <- if ("index_snp" %in% names(ss)) {
    ss[category == "clump_member" & index_snp %in% hi_index_snps]
  } else {
    ss[category == "clump_member"]
  }

  # Layer 1: plot every non-index variant in the chromosome-alternating grey.
  # This is where "other" AND "clump_member" both live — so members whose
  # parent index sits above the highlight threshold appear as regular grey
  # points instead of vanishing. Layer 2 then overlays the highlighted
  # members subset in green on top.
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = ss[category != "index"],
      ggplot2::aes(x = cum_bp, y = neglogP, colour = chr_bg),
      size = point_size) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_point(
      data = ss_members,
      ggplot2::aes(x = cum_bp, y = neglogP),
      colour = col_member, size = point_size) +
    ggplot2::geom_hline(yintercept = -log10(sig_threshold),
                        colour = "red", linetype = "dashed") +
    (if (isTRUE(suggestive_line))
       ggplot2::geom_hline(yintercept = -log10(suggestive_threshold),
                           colour = "blue", linetype = "dotted")
     else NULL) +
    ggplot2::scale_x_continuous(breaks = offsets$label_pos,
                                labels = offsets$CHR,
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    # Clip the y-axis at -log10(p_threshold) so we don't leave a huge empty
    # band below the least-significant plotted variant. p_threshold is the
    # pipeline's Manhattan filter (default 1e-3); md$p_threshold is 1 in the
    # fallback case where < 10 variants passed and all were retained.
    ggplot2::coord_cartesian(
      ylim = c(-log10(md$p_threshold %||% 1e-3), NA), clip = "off") +
    ggplot2::labs(x = "Chromosome",
                  y = expression(-log[10](italic(P))),
                  title = if (nzchar(title)) title else NULL) +
    theme_fn(base_size = font_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Highlight independent leads below the user's threshold as diamonds.
  if (!is.null(idx) && nrow(idx) > 0) {
    idx_hi <- idx[P < highlight_threshold]
    if (nrow(idx_hi) > 0) {
      p <- p + ggplot2::geom_point(
        data = idx_hi, ggplot2::aes(x = cum_bp, y = neglogP),
        shape = 23, size = point_size * 3.5,
        fill = col_index, colour = "black")
    }
    if (isTRUE(label_variants)) {
      # Optionally label with the nearest PROTEIN-CODING gene (label_pc) when the
      # bundle carries it; otherwise fall back to the nearest gene of any biotype.
      idx <- data.table::copy(idx)
      if (isTRUE(protein_coding_labels) && "label_pc" %in% names(idx)) {
        idx[, label := label_pc]
      }
      idx_lab <- idx[!is.na(label) & P < label_threshold]
      if (nrow(idx_lab) > 0) {
        p <- p + ggrepel::geom_text_repel(
          data = idx_lab,
          ggplot2::aes(x = cum_bp, y = neglogP, label = label),
          ylim = c(-log10(sig_threshold), NA),
          direction = "both",
          force = 2,
          box.padding = 0.5,
          point.padding = 0.3,
          min.segment.length = 0,
          max.overlaps = Inf,
          segment.size = 0.3,
          segment.colour = "grey50",
          size = font_size * 0.28)
      }
    }
  }
  p
}

#' SNP Associations module UI
#'
#' @param id Module namespace id
#' @return UI elements for the SNP Associations tab
snpAssocUI <- function(id) {
  ns <- NS(id)
  # Per-tab GWAS picker — rendered via `uiOutput` so the choices are
  # baked in at widget-creation time (see mod_compare.R commit fa9d73e).
  # `updateSelectInput` sent to a selectize-backed selectInput before
  # the client-side binding is ready silently drops the choices update
  # here, leaving the widget with only its default one-item state.
  # No wrapping div here — the picker's outer div (with max-width and
  # margin) is emitted from the renderUI itself, so an empty return
  # in single-GWAS mode leaves ZERO residual vertical space.
  gwas_picker <- function(output_id) uiOutput(ns(output_id))

  tabPanel(
    title = "SNP Associations",
    br(),
    p("This tab shows SNP association results. Manhattan plot and per-GWAS ",
      "content include a GWAS selector at the top. The Lead variants and ",
      "Fine-mapping sub-tabs appear when those analyses (LD-based clumping / ",
      "COJO, and SuSiE fine-mapping) were run for the job."),
    hr(),
    tabsetPanel(
      id = ns("snp_assoc_tabs"),
      tabPanel(
        title = "Manhattan plot",
        value = "manhattan",
        br(),
        gwas_picker("manhattan_gwas_ui"),
        uiOutput(ns("manhattan_desc")),
        hr(),
        uiOutput(ns("manhattan_options_ui")),
        br(),
        uiOutput(ns("manhattan_plot_ui"))
      ),
      tabPanel(
        title = "Lead variants",
        value = "lead_variants",
        br(),
        tabsetPanel(
          id = ns("lead_sub_tabs"),
          tabPanel(
            title = "Single GWAS",
            value = "single",
            br(),
            gwas_picker("lead_gwas_ui"),
            uiOutput(ns("lead_single_body"))
          ),
          tabPanel(
            title = "Multi-GWAS",
            value = "multi",
            br(),
            locus_compare_ui(NS(ns("locus_compare")))
          )
        )
      ),
      tabPanel(
        title = "Fine-mapping",
        value = "finemapping",
        br(),
        gwas_picker("finemap_gwas_ui"),
        fluidPage(
          sidebarPanel(
            radioButtons(ns("l_param"), "Select L parameter:",
                         choices = c("L1" = "L1",
                                     "L10" = "L10"),
                         selected = "L1"),
            width = 3
          ),
          mainPanel(
            uiOutput(ns("finemap_status_message")),
            br(),
            dataTableOutput(ns("snp_assoc_finemap_table")),
            gd_legend(list(
              "L parameter" = paste0(
                "Maximum number of causal signals SuSiE may fit per locus ",
                "(L1 = single causal variant; L10 = up to ten)."),
              "cs" = "Credible-set identifier: a set of variants that jointly capture one causal signal.",
              "NSNP" = "Number of variants in the credible set (smaller = better resolved).",
              "TopPIP" = "Highest posterior inclusion probability in the set (closer to 1 = more confident).",
              "Gene" = "Gene(s) whose body contains all of the variants in the 95% credible set ('None' if no gene spans them).",
              "CALDERA top gene" = "CALDERA's highest-probability candidate causal gene for this credible set, at the selected L.",
              "CALDERA prob" = "That gene's CALDERA causal-gene probability (0–1, relative within the credible set)."
            ), heading = "Column guide"),
            width = 9
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
snpAssocServer <- function(id, gwas_data, selected_gwas, config_flags,
                            selected_gwas_multi = NULL,
                            comparison_mode = NULL,
                            comparison_long = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Cross-GWAS locus compare sub-module. Registered unconditionally; it
    # self-guards on its own reactives.
    if (!is.null(selected_gwas_multi) && !is.null(comparison_long)) {
      locus_compare_server("locus_compare",
                             gwas_data, selected_gwas_multi, comparison_long)
    }

    # Per-tab GWAS pickers, rendered via renderUI so the choice list is
    # baked in at widget-creation time. The wrapping tabPanel is created
    # eagerly (static UI), but the underlying selectize widget only
    # binds on first tab visit — updateSelectInput sent before that bind
    # is silently dropped, leaving the picker with an empty option set.
    .render_per_gwas_picker <- function(output_id, input_id, label = "GWAS:") {
      output[[output_id]] <- renderUI({
        gwas_vec <- selected_gwas_multi()
        req(length(gwas_vec) >= 1)
        # Suppress the picker for single-GWAS bundles.
        if (length(gwas_vec) < 2L) return(NULL)
        cur <- isolate(input[[input_id]])
        keep <- if (!is.null(cur) && cur %in% gwas_vec) cur else gwas_vec[1L]
        choices <- gwas_choices(gwas_data(), gwas_vec)
        # Width sized to the longest label so the widget shows the full
        # GWAS name and the dropdown arrow doesn't wrap to a new line.
        widest <- max(nchar(names(choices)), 12L, na.rm = TRUE)
        w_px   <- max(320L, as.integer(round(widest * 8.5)) + 60L)
        div(style = sprintf("width: %dpx; margin-bottom: 10px;", w_px),
          selectInput(session$ns(input_id), label,
                       choices = choices, selected = keep,
                       multiple = FALSE, width = "100%"))
      })
    }

    # Hide the Multi-GWAS sub-tab under Lead variants when the bundle
    # has only one GWAS — no cross-trait comparison makes sense with 1.
    # Also hide the sub-tab nav-tabs row entirely so the lone Single
    # GWAS tab-header doesn't clutter the UI (the .hide-inner-nav CSS
    # class defined in app.R hides `> ul.nav-tabs`).
    observe({
      is_multi <- length(selected_gwas_multi()) > 1L
      if (is_multi) {
        showTab("lead_sub_tabs", "multi", session = session)
        shinyjs::removeClass(id = "lead_sub_tabs", class = "hide-inner-nav")
      } else {
        hideTab("lead_sub_tabs", "multi", session = session)
        cur <- isolate(input$lead_sub_tabs)
        if (!is.null(cur) && cur == "multi") {
          updateTabsetPanel(session, "lead_sub_tabs", selected = "single")
        }
        shinyjs::addClass(id = "lead_sub_tabs", class = "hide-inner-nav")
      }
    })
    if (!is.null(selected_gwas_multi)) {
      .render_per_gwas_picker("manhattan_gwas_ui", "manhattan_gwas")
      .render_per_gwas_picker("lead_gwas_ui",      "lead_gwas")
      .render_per_gwas_picker("finemap_gwas_ui",   "finemap_gwas")
    }

    # Fallbacks used by outputs before their per-tab observer has fired.
    .pick_gwas <- function(input_id) {
      v <- input[[input_id]]
      if (is.null(v) || !nzchar(v)) selected_gwas() else v
    }

    # Lead variants — single-GWAS body. Kept in a renderUI so the
    # method-radio choices (LD-based clumping / COJO) reflect what the
    # pipeline actually ran; baked in at creation time to avoid the
    # race that used to leak COJO into the dropdown when COJO wasn't run.
    output$lead_single_body <- renderUI({
      cf <- config_flags()
      lead_choices <- c()
      if (isTRUE(cf$cojo))  lead_choices <- c(lead_choices, "COJO" = "cojo_analysis")
      if (isTRUE(cf$clump)) lead_choices <- c(lead_choices, "LD-based clumping" = "ld_clumping")
      fluidPage(
        sidebarPanel(
          radioButtons(ns("clumping_type"), "Select method:",
                       choices = lead_choices,
                       selected = if (length(lead_choices) > 0) unname(lead_choices[1]) else character(0)),
          radioButtons(ns("pvalue_threshold"), "Select P-value Threshold:",
                       choices = c("Genome-wide significance (p < 5e-8)" = 5e-8),
                       selected = 5e-8),
          uiOutput(ns("lead_pc_ui")),
          width = 3
        ),
        mainPanel(
          uiOutput(ns("cojo_status_message")),
          uiOutput(ns("lead_status_message")),
          dataTableOutput(ns("snp_assoc_lead_table")),
          uiOutput(ns("lead_legend")),
          tags$hr(),
          tags$h4("Locus plot"),
          uiOutput(ns("lp_status")),
          uiOutput(ns("lp_controls_ui")),
          tags$div(style = "max-width: 1100px;",
                   plotOutput(ns("lp_plot"), height = "auto")),
          uiOutput(ns("lp_legends")),
          tags$hr(),
          tags$h5("Gene association summary (genes in positional order)"),
          uiOutput(ns("lp_grid_controls_ui")),
          uiOutput(ns("lp_grid_ui")),
          uiOutput(ns("lp_grid_legend")),
          tags$div(style = "height: 350px;"),   # space so the figure can scroll to the bottom
          width = 9
        )
      )
    })

    # Tab visibility now depends only on which pipeline outputs exist; no
    # more compare-mode hiding since Single / Multi are sub-tabs of each
    # per-GWAS section.
    apply_snp_assoc_visibility <- function() {
      cf <- config_flags()
      if (is.null(cf)) return()
      # Lead variants tab appears when either clumping or COJO was run.
      if (any(cf$clump, cf$cojo)) {
        showTab("snp_assoc_tabs", "lead_variants", session = session)
      } else {
        hideTab("snp_assoc_tabs", "lead_variants", session = session)
      }
      # Fine-mapping tab appears when the pipeline ran SuSiE.
      if (isTRUE(cf$finemap)) {
        showTab("snp_assoc_tabs", "finemapping", session = session)
      } else {
        hideTab("snp_assoc_tabs", "finemapping", session = session)
      }
    }

    observeEvent(config_flags(), apply_snp_assoc_visibility())

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

    # Convenience: pull the whole gwas_qc block once per bundle/GWAS change.
    manhattan_qc <- reactive({
      req(gwas_data())
      g <- .pick_gwas("manhattan_gwas"); req(g)
      gd_read(gwas_data(), g, "gwas_qc")
    })
    manhattan_data <- reactive({ manhattan_qc()$manhattan_data })

    # Plot options collapsible — only shown when we have raw manhattan_data
    # (bundles produced before the pipeline change fall through to the
    # legacy PNG fallback in manhattan_plot_ui).
    # Per-plot description. Clump-derived features (green members, lead
    # diamonds, gene labels) are described only when clumping ran; otherwise
    # we note that those features require clumping.
    output$manhattan_desc <- renderUI({
      has_clump <- isTRUE(config_flags()$clump)
      base <- tagList(
        "Genome-wide summary of the association signal. Each point is a variant, ",
        "positioned by chromosome (x) and ", tags$code("-log10(P)"), " (y)."
      )
      if (has_clump) {
        p(base,
          " Green points are variants clumped with an index (lead) variant; ",
          "dark-green diamonds are the index variants themselves. Use ",
          tags$b("Plot options"),
          " to customise thresholds, gene labels, appearance, and download the figure.")
      } else {
        tagList(
          p(base, " Use ", tags$b("Plot options"),
            " to customise the significance threshold and appearance, and download the figure."),
          tags$p(style = "color: var(--gd-text-mute);",
            tags$em("Independent lead-variant highlighting and nearest-gene labels are ",
                    "only available when LD-based clumping was selected for this job."))
        )
      }
    })

    output$manhattan_options_ui <- renderUI({
      if (is.null(manhattan_data())) return(NULL)
      # Independent-lead highlighting and nearest-gene labels are derived from
      # LD-based clumping only, so those controls (and the member/lead colour
      # pickers) are shown only when clumping ran for this job.
      has_clump <- isTRUE(config_flags()$clump)
      intro_text <- if (has_clump) {
        paste0("Customise thresholds, gene labels, appearance (title, theme, ",
               "font size, point size), and download the figure as PNG / PDF / SVG.")
      } else {
        paste0("Customise the significance threshold and appearance (title, theme, ",
               "font size, point size), and download the figure as PNG / PDF / SVG.")
      }
      tags$details(class = "gd-details",
        tags$summary("Plot options"),
        tags$div(class = "gd-details-body",
          tags$p(class = "gd-details-intro", intro_text),
          fluidRow(
            column(4,
              textInput(ns("mh_title"), "Plot title (optional):", value = ""),
              selectInput(ns("mh_theme"), "Theme:",
                          choices = c("Black & white" = "bw", "Minimal" = "minimal",
                                      "Classic" = "classic", "Light" = "light"),
                          selected = "bw"),
              sliderInput(ns("mh_font_size"), "Font size (pt):",
                          min = 8, max = 20, value = 14, step = 1),
              sliderInput(ns("mh_point_size"), "Point size:",
                          min = 0.3, max = 3, value = 0.9, step = 0.1)
            ),
            column(4,
              textInput(ns("mh_sig_threshold"),
                        "Significance line at P <:", value = "5e-8"),
              checkboxInput(ns("mh_suggestive"),
                            "Show suggestive line", value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("mh_suggestive")),
                textInput(ns("mh_suggestive_threshold"),
                          "Suggestive line at P <:", value = "1e-5")
              ),
              if (has_clump) tagList(
                textInput(ns("mh_highlight_threshold"),
                          "Highlight independent leads at P < (max 1e-5):", value = "5e-8",
                          placeholder = "max 1e-5 (the clumping P threshold)"),
                checkboxInput(ns("mh_label"),
                              "Label nearest gene at lead variants", value = FALSE),
                conditionalPanel(
                  condition = sprintf("input['%s'] == true", ns("mh_label")),
                  textInput(ns("mh_label_threshold"),
                            "Label leads with P < (max 1e-5):", value = "5e-8",
                            placeholder = "max 1e-5 (the clumping P threshold)"),
                  checkboxInput(ns("mh_label_pc"),
                                "Protein-coding genes only", value = FALSE)
                )
              ) else NULL
            ),
            column(4,
              selectInput(ns("mh_dl_format"), "Download format:",
                          choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                          selected = "png"),
              selectInput(ns("mh_dl_units"), "Units:",
                          choices = .dl_units_choices, selected = "in"),
              numericInput(ns("mh_dl_width"), "Width:",
                           value = 12, min = 0.1, max = 20000, step = 0.5),
              numericInput(ns("mh_dl_height"), "Height:",
                           value = 5, min = 0.1, max = 20000, step = 0.5),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'png'", ns("mh_dl_format")),
                numericInput(ns("mh_dl_dpi"), "Resolution (DPI, PNG only):",
                             value = 300, min = 72, max = 600, step = 25)
              ),
              downloadButton(ns("mh_download"), "Download plot")
            )
          ),
          br(),
          fluidRow(
            column(3, colourpicker::colourInput(
              ns("mh_col_chr_odd"),  "Colour: odd chromosomes",  value = "#8c8c8c")),
            column(3, colourpicker::colourInput(
              ns("mh_col_chr_even"), "Colour: even chromosomes", value = "#4d4d4d")),
            if (has_clump) column(3, colourpicker::colourInput(
              ns("mh_col_member"),   "Colour: clumped members",  value = "#2ca02c")) else NULL,
            if (has_clump) column(3, colourpicker::colourInput(
              ns("mh_col_index"),    "Colour: independent leads", value = "#006400")) else NULL
          ),
          tags$div(style = "text-align: right; margin-top: 4px;",
            actionButton(ns("mh_reset_colours"), "Reset colours")
          )
        )
      )
    })

    output$manhattan_plot_ui <- renderUI({
      req(gwas_data())
      req(.pick_gwas("manhattan_gwas"))
      gwas_qc <- manhattan_qc()

      # Preferred: raw data → modifiable ggplot render.
      if (!is.null(gwas_qc$manhattan_data)) {
        p_thr <- gwas_qc$manhattan_data$p_threshold %||% 1e-3
        variant_filter_note <- if (p_thr < 1) {
          sprintf(
            "Only variants with P < %s are plotted (a pipeline-side filter that keeps the plot lightweight). The y-axis starts at -log10(%s) = %.2f so the shown range matches what's actually plotted.",
            format(p_thr, scientific = TRUE), format(p_thr, scientific = TRUE), -log10(p_thr))
        } else {
          "All variants are plotted (fewer than 10 passed the usual P < 1e-3 filter, so the pipeline fell back to showing everything)."
        }
        legend_items <- list(
          "Axes" = "Each point is a variant, positioned by chromosome (x) and -log10(P) (y); higher points are more strongly associated.",
          "Variant filter" = variant_filter_note,
          "Red dashed line" = "Significance line (default P < 5e-8; adjustable in Plot options)."
        )
        # Clump-derived marks only exist when clumping ran.
        if (isTRUE(config_flags()$clump)) {
          legend_items <- c(legend_items, list(
            "Green points" = "Variants clumped with an index variant (part of the same association signal).",
            "Dark-green diamonds" = "Index (lead) variants below the highlight threshold (adjustable in Plot options).",
            "Gene labels" = "Nearest gene to each labelled index variant (toggle and threshold in Plot options)."
          ))
        }
        return(tagList(
          tags$div(style = "max-width: 1200px;",
            plotOutput(ns("manhattan_plot"), height = "500px")
          ),
          gd_legend(legend_items, heading = "How to read this plot")
        ))
      }

      # Legacy fallback: pre-manhattan_data bundles still ship the base64 PNGs.
      b64 <- gwas_qc$manhattan_plot_unlabelled_base64
      if (is.null(b64)) {
        return(div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          icon("chart-bar", style = "font-size: 2em;"),
          br(), br(),
          tags$strong("Manhattan plot Unavailable"),
          br(),
          "This results package was produced before the Manhattan plot rule was added."
        ))
      }
      tagList(
        tags$img(
          src = paste0("data:image/png;base64,", b64),
          style = "max-width: 100%; height: auto;"
        ),
        tags$p(style = "font-size: 0.85em; color: var(--gd-text-mute); margin-top: 8px;",
               tags$em("This bundle predates the modifiable Manhattan plot; ",
                       "regenerate it with the current pipeline to unlock plot options."))
      )
    })

    output$manhattan_plot <- renderPlot({
      md <- manhattan_data()
      if (is.null(md)) return(NULL)
      build_manhattan_plot(
        md,
        sig_threshold        = parse_p(input$mh_sig_threshold,        5e-8),
        suggestive_line      = isTRUE(input$mh_suggestive),
        suggestive_threshold = parse_p(input$mh_suggestive_threshold, 1e-5),
        # Cap at 1e-5 — the clumping P threshold. No index variant with
        # P above that exists in the data, so a looser value is meaningless.
        highlight_threshold  = min(parse_p(input$mh_highlight_threshold, 5e-8), 1e-5),
        label_variants      = isTRUE(input$mh_label),
        # Cap at 1e-5 as above — no labelled index variant exists above it.
        label_threshold     = min(parse_p(input$mh_label_threshold, 5e-8), 1e-5),
        protein_coding_labels = isTRUE(input$mh_label_pc),
        font_size           = input$mh_font_size  %||% 14,
        point_size          = input$mh_point_size %||% 0.9,
        theme_fn            = mol_theme_fn(input$mh_theme),
        title               = input$mh_title      %||% "",
        col_chr_odd  = input$mh_col_chr_odd  %||% "#8c8c8c",
        col_chr_even = input$mh_col_chr_even %||% "#4d4d4d",
        col_member   = input$mh_col_member   %||% "#2ca02c",
        col_index    = input$mh_col_index    %||% "#006400"
      )
    })

    mh_dl_units_prev <- reactiveVal("in")
    dl_units_auto_convert(input, session,
                           units_id = "mh_dl_units",
                           width_id = "mh_dl_width",
                           height_id = "mh_dl_height",
                           prev_val = mh_dl_units_prev,
                           dpi_id   = "mh_dl_dpi",
                           default_dpi = 300)

    output$mh_download <- downloadHandler(
      filename = function() {
        sprintf("manhattan_%s.%s",
                format(Sys.time(), "%Y%m%d_%H%M%S"),
                input$mh_dl_format %||% "png")
      },
      content = function(file) {
        md <- manhattan_data()
        if (is.null(md)) {
          grDevices::png(file, width = 4, height = 1, units = "in", res = 96)
          grid::grid.text("No manhattan_data in this bundle.")
          grDevices::dev.off()
          return(invisible())
        }
        p <- build_manhattan_plot(
          md,
          sig_threshold       = parse_p(input$mh_sig_threshold,       5e-8),
          # Cap at 1e-5 — the clumping P threshold — so the saved figure
          # matches the on-screen plot (see renderPlot above).
          highlight_threshold = min(parse_p(input$mh_highlight_threshold, 5e-8), 1e-5),
          label_variants      = isTRUE(input$mh_label),
          label_threshold     = min(parse_p(input$mh_label_threshold, 5e-8), 1e-5),
          font_size           = input$mh_font_size  %||% 14,
          point_size          = input$mh_point_size %||% 0.9,
          theme_fn            = mol_theme_fn(input$mh_theme),
          title               = input$mh_title      %||% "",
          col_chr_odd  = input$mh_col_chr_odd  %||% "#8c8c8c",
          col_chr_even = input$mh_col_chr_even %||% "#4d4d4d",
          col_member   = input$mh_col_member   %||% "#2ca02c",
          col_index    = input$mh_col_index    %||% "#006400"
        )
        open_plot_device(file,
          fmt   = input$mh_dl_format %||% "png",
          w     = input$mh_dl_width  %||% 12,
          h     = input$mh_dl_height %||% 5,
          units = input$mh_dl_units  %||% "in",
          dpi   = input$mh_dl_dpi    %||% 300)
        print(p)
        grDevices::dev.off()
      }
    )

    # Reset the four Manhattan colour pickers to the defaults declared in
    # the UI / build_manhattan_plot(). Kept in sync manually — small enough
    # not to warrant a shared constant.
    observeEvent(input$mh_reset_colours, {
      colourpicker::updateColourInput(session, "mh_col_chr_odd",  value = "#8c8c8c")
      colourpicker::updateColourInput(session, "mh_col_chr_even", value = "#4d4d4d")
      colourpicker::updateColourInput(session, "mh_col_member",   value = "#2ca02c")
      colourpicker::updateColourInput(session, "mh_col_index",    value = "#006400")
    })

    snp_assoc_lead_data <- reactive({
      req(gwas_data(), input$clumping_type)
      g <- .pick_gwas("lead_gwas"); req(g)
      snp_assoc <- gd_read(gwas_data(), g, "snp_assoc")
      snp_assoc_lead <- NULL
      if (input$clumping_type == "ld_clumping") {
        snp_assoc_lead <- snp_assoc$clump
      }

      if (input$clumping_type == "cojo_analysis") {
        snp_assoc_lead <- snp_assoc$cojo
      }

      req(snp_assoc_lead)

      # FREQ = effect-allele frequency reported in the GWAS; REF.FREQ = allele
      # frequency in the ancestry-matched reference panel (both carried through
      # from the cleaned sumstats, if present). NearestGene kept last.
      snp_assoc_lead <- snp_assoc_lead[, names(snp_assoc_lead) %in% c("CHR","BP","SNP","A1","A2","BETA","SE","P","FREQ","REF.FREQ","NearestGene","NearestGeneProteinCoding"), with = F]

      # Optionally show the nearest PROTEIN-CODING gene instead (companion column
      # from process_clump.R; absent in older bundles and for COJO). Swap into the
      # single NearestGene column so the table layout is unchanged, then drop it.
      if ("NearestGeneProteinCoding" %in% names(snp_assoc_lead)) {
        if (isTRUE(input$lead_pc_only))
          snp_assoc_lead[, NearestGene := NearestGeneProteinCoding]
        snp_assoc_lead[, NearestGeneProteinCoding := NULL]
      }

      snp_assoc_lead$P <- as.numeric(snp_assoc_lead$P)

      snp_assoc_lead <- snp_assoc_lead[snp_assoc_lead$P < as.numeric(input$pvalue_threshold), ]

      snp_assoc_lead$BETA <- round(snp_assoc_lead$BETA, 2)
      snp_assoc_lead$SE <- round(snp_assoc_lead$SE, 2)
      for (cc in c("FREQ", "REF.FREQ"))
        if (cc %in% names(snp_assoc_lead))
          snp_assoc_lead[[cc]] <- round(as.numeric(snp_assoc_lead[[cc]]), 4)

      snp_assoc_lead <- snp_assoc_lead[order(snp_assoc_lead$CHR, snp_assoc_lead$BP), ]

      return(snp_assoc_lead)
    })

    # COJO can fail per chromosome when the number of independent genome-wide-significant
    # signals exceeds the LD reference sample size (~503, 1000 Genomes EUR). package_results
    # carries snp_assoc$cojo_status so we can flag partial results, or a full failure.
    output$cojo_status_message <- renderUI({
      req(gwas_data(), input$clumping_type)
      g <- .pick_gwas("lead_gwas"); req(g)
      if (input$clumping_type != "cojo_analysis") return(NULL)

      cojo_status <- gd_read(gwas_data(), g, "snp_assoc")$cojo_status
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

    # Distinct from the COJO reference-too-small banner above: this fires when the
    # analysis ran fine but simply found no variants passing the threshold (the
    # common underpowered-GWAS case). Suppressed when the COJO failure banner shows.
    output$lead_status_message <- renderUI({
      req(gwas_data(), input$clumping_type)
      g <- .pick_gwas("lead_gwas"); req(g)

      if (input$clumping_type == "cojo_analysis") {
        cojo_status <- gd_read(gwas_data(), g, "snp_assoc")$cojo_status
        if (!is.null(cojo_status) && isTRUE(cojo_status$any_failed)) return(NULL)
      }

      dat <- snp_assoc_lead_data()
      if (nrow(dat) > 0) return(NULL)

      thr <- if (isTRUE(as.numeric(input$pvalue_threshold) == 1e-5)) "1e-5" else "5e-8"
      method <- if (input$clumping_type == "cojo_analysis") "COJO" else "LD-based clumping"
      div(
        style = "background-color: #e9ecef; border-radius: 8px; padding: 30px 20px; text-align: center; color: #6c757d;",
        icon("info-circle", style = "font-size: 1.5em;"),
        br(), br(),
        tags$strong("No lead variants"),
        br(),
        paste0("No variants passed the p < ", thr, " threshold for ", method, ". ",
               "This is expected for an underpowered GWAS with no strong associations.")
      )
    })

    # Column guide for the lead-variant table, specific to the selected method.
    # Both tables show the same columns, but the meaning of a "lead variant" and
    # of BETA/SE/P differs between the joint COJO model and marginal LD clumping.
    output$lead_legend <- renderUI({
      req(input$clumping_type)
      common <- list(
        "CHR / BP" = "Chromosome and base-pair position of the lead variant (genome build as identified on the QC Summary tab).",
        "SNP" = "Variant identifier (rsID).",
        "A1 / A2" = "Effect allele (A1) and other allele (A2); BETA is expressed per copy of A1.",
        "FREQ" = "A1 (effect-allele) frequency as reported in the GWAS, where available.",
        "REF.FREQ" = "A1 frequency in the ancestry-matched reference panel used for QC and clumping, where available.",
        "NearestGene" = "Gene nearest to the lead variant (not necessarily the causal gene)."
      )
      if (input$clumping_type == "cojo_analysis") {
        items <- c(list(
          "Method (COJO)" = paste0(
            "COJO (Manc-COJO) runs a stepwise conditional & joint analysis, keeping variants that stay ",
            "genome-wide significant after conditioning on each other - i.e. approximately ",
            "independent association signals at a locus.")),
          common,
          list("BETA / SE" = "Joint effect size from the COJO model and its standard error.",
               "P" = "Joint conditional association p-value (already thresholded at p < 5e-8)."))
      } else {
        items <- c(list(
          "Method (LD-based clumping)" = paste0(
            "PLINK clumping groups variants correlated (in LD) with a top 'index' variant into ",
            "one signal; each row is an index variant representing an approximately independent locus.")),
          common,
          list("BETA / SE" = "Marginal (single-variant) GWAS effect size and standard error.",
               "P" = "Marginal GWAS association p-value of the index variant."))
      }
      gd_legend(items, heading = "Column guide")
    })

    # Protein-coding toggle for the lead-variant table — shown only for LD
    # clumping and only when the bundle carries the NearestGeneProteinCoding
    # column with data (older bundles / COJO have none).
    output$lead_pc_ui <- renderUI({
      req(gwas_data(), input$clumping_type)
      if (input$clumping_type != "ld_clumping") return(NULL)
      g <- .pick_gwas("lead_gwas"); req(g)
      cl <- tryCatch(gd_read(gwas_data(), g, "snp_assoc")$clump, error = function(e) NULL)
      if (is.null(cl) || !("NearestGeneProteinCoding" %in% names(cl))) return(NULL)
      vals <- as.character(cl$NearestGeneProteinCoding)
      if (all(is.na(vals) | vals %in% c("", "NA"))) return(NULL)
      checkboxInput(ns("lead_pc_only"),
                    "Restrict nearest gene to protein-coding genes", value = FALSE)
    })

    output$snp_assoc_lead_table <- renderDataTable({
      d <- snp_assoc_lead_data()
      # Column indices are resolved by name (0-indexed) because the column set
      # differs between clumping (carries FREQ / REF.FREQ) and COJO.
      p_idx  <- match("P", names(d)) - 1L
      ng_idx <- match("NearestGene", names(d)) - 1L
      last   <- ncol(d) - 1L
      js <- c(
        "function(row, data, displayNum, index){",
        sprintf("  var x = data[%d];", p_idx),
        sprintf("  if (x !== null && x !== '') $('td:eq(%d)', row).html((+x).toExponential(2));", p_idx),
        "}"
      )

      datatable(d, rownames = F, width = 7,
        selection = list(mode = 'single', target = 'row'),
        options = list(
          rowCallback = JS(js),
          columnDefs = c(
            list(list(className = 'dt-center', targets = 0:last),
                 list(width = '60px', targets = p_idx)),
            if (!is.na(ng_idx)) list(list(width = '600px', targets = ng_idx)))))
    })

    snp_assoc_finemap_data <- reactive({
      req(gwas_data(), input$l_param)
      g <- .pick_gwas("finemap_gwas"); req(g)
      snp_assoc <- gd_read(gwas_data(), g, "snp_assoc")
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

      # Append CALDERA's top gene + probability for the SAME L (feedback i.2).
      # CALDERA is keyed by (locus = lead variant, cs, L); the fine-mapping
      # table's SNP column is that lead variant, so join on (SNP, cs).
      cal <- tryCatch(gd_read(gwas_data(), g, "caldera"), error = function(e) NULL)
      if (!is.null(cal) && nrow(cal) > 0 && "SNP" %in% names(snp_assoc_finemap)) {
        cal <- as.data.frame(cal)
        Lval <- if (identical(input$l_param, "L10")) 10L else 1L
        if ("L" %in% names(cal)) cal <- cal[suppressWarnings(as.integer(cal$L)) == Lval, , drop = FALSE]
        cal$caldera <- suppressWarnings(as.numeric(cal$caldera))
        cal <- cal[is.finite(cal$caldera), , drop = FALSE]
        if (nrow(cal) > 0) {
          top <- do.call(rbind, lapply(
            split(cal, list(cal$locus, cal$cs), drop = TRUE), function(x) {
              r <- x[which.max(x$caldera), , drop = FALSE]
              # Canonicalise the CALDERA (GENCODE) symbol so it matches the
              # biomart-derived Gene column spelling in the same table.
              ck <- if ("ensgid" %in% names(r) && !is.na(r$ensgid) && r$ensgid != "")
                as.character(r$ensgid) else as.character(r$gene)
              ch <- harmonise_genes(ck)
              cdisp <- if (!is.na(ch$symbol) && !grepl("^ENSG", ch$symbol)) ch$symbol else as.character(r$gene)
              data.frame(SNP = as.character(r$locus), cs = as.integer(r$cs),
                         `CALDERA top gene` = cdisp,
                         `CALDERA prob` = round(r$caldera, 3),
                         check.names = FALSE, stringsAsFactors = FALSE)
            }))
          snp_assoc_finemap$cs <- as.integer(snp_assoc_finemap$cs)
          snp_assoc_finemap <- merge(snp_assoc_finemap, top,
                                     by = c("SNP", "cs"), all.x = TRUE, sort = FALSE)
        }
      }

      return(snp_assoc_finemap)
    })

    # When there are no credible sets, the pipeline emits a single all-NA placeholder
    # row rather than an empty table; treat that (and a genuinely empty table) as empty.
    finemap_is_empty <- function(df) {
      is.null(df) || nrow(df) == 0 || all(is.na(df$SNP))
    }

    output$finemap_status_message <- renderUI({
      req(gwas_data(), input$l_param)
      req(.pick_gwas("finemap_gwas"))
      if (!finemap_is_empty(snp_assoc_finemap_data())) return(NULL)
      div(
        style = "background-color: #e9ecef; border-radius: 8px; padding: 30px 20px; text-align: center; color: #6c757d;",
        icon("info-circle", style = "font-size: 1.5em;"),
        br(), br(),
        tags$strong("No credible sets identified"),
        br(),
        paste0("Fine-mapping runs only at loci containing a variant at p < 5e-8. ",
               "This GWAS has no such loci, so no credible sets were produced.")
      )
    })

    output$snp_assoc_finemap_table <- renderDataTable({
      df <- as.data.frame(snp_assoc_finemap_data())
      req(!finemap_is_empty(df))
      # merge() moves the join keys first; restore a readable column order.
      # (as.data.frame above: susie blocks are data.tables, where `dt[, c(cols)]`
      # would return the name vector rather than reordered columns.)
      ord <- c("CHR", "BP", "SNP", "cs", "NSNP", "TopPIP", "Gene",
               "CALDERA top gene", "CALDERA prob")
      df <- df[, c(intersect(ord, names(df)), setdiff(names(df), ord)), drop = FALSE]
      # Canonicalise each gene in the (possibly comma-separated) biomart Gene
      # column so its spelling matches the CALDERA top gene column.
      if ("Gene" %in% names(df)) {
        gl <- strsplit(as.character(df$Gene), ",\\s*")
        df$Gene <- vapply(gl, function(parts) {
          parts <- parts[!is.na(parts) & nzchar(parts)]
          if (!length(parts)) return(NA_character_)
          h <- harmonise_genes(parts)
          disp <- ifelse(!is.na(h$symbol) & !grepl("^ENSG", h$symbol), h$symbol, parts)
          paste(unique(disp), collapse = ", ")
        }, character(1))
        # Parallel biotype column, only when the reference carries biotype.
        bt <- vapply(gl, function(parts) {
          parts <- parts[!is.na(parts) & nzchar(parts)]
          if (!length(parts)) return(NA_character_)
          b <- harmonise_genes(parts)$biotype
          if (all(is.na(b))) NA_character_ else paste(unique(b[!is.na(b)]), collapse = ", ")
        }, character(1))
        if (any(!is.na(bt))) {
          bt[is.na(bt)] <- "—"
          pos <- match("Gene", names(df))
          df <- data.frame(df[, seq_len(pos), drop = FALSE],
                           `Gene biotype` = bt,
                           df[, setdiff(seq_len(ncol(df)), seq_len(pos)), drop = FALSE],
                           check.names = FALSE, stringsAsFactors = FALSE)
        }
      }
      datatable(df, rownames = F, width = 7, selection = "none", options = list(
        columnDefs = list(list(className = 'dt-center', targets = "_all"))))
    })

    # ======================================================================
    # Locus plot (enriched fine-mapping dashboard) — builders in mod_locus_plot.R
    # ======================================================================
    # Fine-mapped loci (SuSiE leads) at the currently-selected dashboard L.
    lp_loci <- reactive({
      req(gwas_data())
      g <- .pick_gwas("lead_gwas"); req(g)
      sa <- tryCatch(gd_read(gwas_data(), g, "snp_assoc"), error = function(e) NULL)
      st <- if (!is.null(sa) && !is.null(sa$susie)) sa$susie[[input$lp_l %||% "L1"]] else NULL
      if (is.null(st) || nrow(st) == 0) return(character(0))
      sort(unique(as.character(st$SNP[!is.na(st$SNP)])))
    })
    # Any fine-mapped locus across both L — governs whether the dashboard shows.
    lp_any_finemap <- reactive({
      req(gwas_data())
      g <- .pick_gwas("lead_gwas"); req(g)
      sa <- tryCatch(gd_read(gwas_data(), g, "snp_assoc"), error = function(e) NULL)
      su <- if (!is.null(sa)) sa$susie else NULL
      if (is.null(su)) return(FALSE)
      length(c(su$L1$SNP, su$L10$SNP)) > 0
    })
    # Locus = lead variant of the clicked row in the lead-variant table.
    lp_selected_locus <- reactive({
      sel <- input$snp_assoc_lead_table_rows_selected
      if (is.null(sel) || !length(sel)) return(NULL)
      tbl <- snp_assoc_lead_data()
      if (sel > nrow(tbl)) return(NULL)
      as.character(tbl[sel, "SNP"][[1]])
    })
    # Non-NULL only when the clicked locus actually has a credible set at this L.
    lp_valid_locus <- reactive({
      sel <- lp_selected_locus(); if (is.null(sel)) return(NULL)
      if (sel %in% lp_loci()) sel else NULL
    })

    # Status / prompt line above the plot.
    output$lp_status <- renderUI({
      if (!isTRUE(lp_any_finemap()))
        return(div(class = "gd-legend", em("No fine-mapped loci (credible sets) in this GWAS.")))
      sel <- lp_selected_locus()
      if (is.null(sel))
        return(div(class = "gd-legend", em("Select a lead-variant row above to view its locus plot.")))
      if (!(sel %in% lp_loci())) {
        l_now <- input$lp_l %||% "L1"; other <- if (l_now == "L1") "L10" else "L1"
        return(div(class = "gd-legend", em(sprintf(
          "%s has no credible set at %s — try %s, or select a fine-mapped lead variant.",
          sel, l_now, other))))
      }
      NULL
    })

    # Controls render once fine-mapping exists (not per-locus, so a row click
    # doesn't rebuild them and reset the user's panel/HC choices). The region
    # inputs are synced per locus by updateNumericInput, below.
    output$lp_controls_ui <- renderUI({
      if (!isTRUE(lp_any_finemap())) return(NULL)
      tagList(
        fluidRow(
          column(3, radioButtons(ns("lp_l"), "L parameter:",
                                 choices = c("L1" = "L1", "L10" = "L10"),
                                 selected = isolate(input$lp_l) %||% "L1", inline = TRUE)),
          column(3, radioButtons(ns("lp_pip_type"), "PIP shown:",
                                 choices = c("Marginal" = "marginal", "Per-cs (alpha)" = "alpha"),
                                 selected = "marginal", inline = TRUE)),
          column(3, checkboxInput(ns("lp_pc_only"), "Protein-coding genes only", FALSE))),
        fluidRow(
          column(4, numericInput(ns("lp_lo"), "Region start (Mb):", value = NA, step = 0.05)),
          column(4, numericInput(ns("lp_hi"), "Region end (Mb):", value = NA, step = 0.05)),
          column(4, tags$div(style = "margin-top:25px;",
                    actionButton(ns("lp_reset_region"), "Reset to locus", class = "btn btn-default btn-sm")))),
        fluidRow(column(12, checkboxGroupInput(ns("lp_hc"), "Mark high-confidence genes:",
                       choices = c("CALDERA ≥ 0.5" = "caldera", "Coding (alpha ≥ 0.5)" = "coding",
                                   "Nearest gene" = "nearest", "Colocalised (E/S/P)" = "coloc"),
                       selected = c("caldera", "coding", "coloc"), inline = TRUE))),
        fluidRow(column(12, checkboxGroupInput(ns("lp_panels"), "Panels to show:",
                       choices = c("Association (Manhattan)" = "A", "Fine-mapping" = "B", "Genes" = "C"),
                       selected = c("A", "B", "C"), inline = TRUE))),
        tags$details(class = "gd-details", style = "margin-top:8px;margin-bottom:10px;",
          tags$summary("Plot options"),
          tags$div(class = "gd-details-body",
            fluidRow(
              column(4, sliderInput(ns("lp_font"), "Font size (pt):", min = 8, max = 22, value = 14, step = 1),
                        selectInput(ns("lp_theme"), "Theme:",
                          choices = c("Black & white"="bw","Minimal"="minimal","Classic"="classic","Light"="light"), selected = "bw")),
              column(4, selectInput(ns("lp_dl_format"), "Download format:", choices = c("PNG"="png","PDF"="pdf","SVG"="svg")),
                        numericInput(ns("lp_dl_dpi"), "DPI (PNG):", value = 300, min = 72, max = 600, step = 1)),
              column(4, numericInput(ns("lp_dl_w"), "Width (in):", value = 10, min = 3, max = 30, step = 0.5),
                        numericInput(ns("lp_dl_h"), "Height (in):", value = 8, min = 2, max = 60, step = 0.5))),
            downloadButton(ns("lp_download"), "Download plot", class = "btn btn-default btn-sm")))
      )
    })

    # Natural window (no override) for defaulting / resetting the BP-limit inputs.
    lp_natural_window <- reactive({
      g <- .pick_gwas("lead_gwas"); req(gwas_data(), g, lp_valid_locus())
      d <- lp_build_locus_data(gwas_data(), g, lp_valid_locus(), input$lp_l %||% "L1", win_override = NULL)
      if (is.null(d)) NULL else d$window
    })
    # Debounce only the typing of the region bounds, so a redraw doesn't fire on
    # every keystroke.
    lp_win_raw <- shiny::debounce(reactive({
      lo <- suppressWarnings(as.numeric(input$lp_lo)); hi <- suppressWarnings(as.numeric(input$lp_hi))
      # length guards first: before the inputs are created as.numeric(NULL) is
      # numeric(0), and is.na(numeric(0)) || ... errors ("missing value").
      if (length(lo) != 1 || length(hi) != 1 || is.na(lo) || is.na(hi) || hi <= lo) NULL else c(lo, hi)
    }), 800)
    # The applied region override, kept in a de-duped reactiveVal so a no-op
    # change doesn't re-run lp_data (which would fade/redraw the plot for no
    # reason). It is NULL when there is no meaningful override — either the
    # region matches the locus's natural window (e.g. just after a locus switch
    # resets the inputs), or the region is stale from a previously-selected
    # locus (doesn't overlap the new window) — so switching locus renders the
    # natural window exactly once, with no second redraw.
    lp_override_rv <- reactiveVal(NULL)
    observe({
      v <- lp_win_raw()
      newval <- if (is.null(v)) NULL else {
        nw <- lp_natural_window()   # only evaluated once a region is set (avoids a startup req-halt)
        if (!is.null(nw) && (v[2] < nw[1] || v[1] > nw[2])) NULL
        else if (!is.null(nw) && isTRUE(all.equal(round(v, 3), round(nw, 3)))) NULL
        else v
      }
      if (!identical(isolate(lp_override_rv()), newval)) lp_override_rv(newval)
    })
    lp_data <- reactive({
      g <- .pick_gwas("lead_gwas"); req(gwas_data(), g, lp_valid_locus())
      lp_build_locus_data(gwas_data(), g, lp_valid_locus(), input$lp_l %||% "L1",
                          pip_type = input$lp_pip_type %||% "marginal",
                          pc_only = isTRUE(input$lp_pc_only),
                          hc_methods = input$lp_hc %||% character(0),
                          win_override = lp_override_rv())
    })
    .lp_set_region <- function() {
      nw <- lp_natural_window(); if (is.null(nw)) return()
      updateNumericInput(session, "lp_lo", value = round(nw[1], 3))
      updateNumericInput(session, "lp_hi", value = round(nw[2], 3))
    }
    # Reset the region whenever the locus (clicked row) or L changes.
    observeEvent(list(lp_valid_locus(), input$lp_l), .lp_set_region())
    observeEvent(input$lp_reset_region, .lp_set_region())

    lp_nrows <- reactive({
      if (is.null(lp_valid_locus())) return(1L)   # keep the plot-height fn from erroring when nothing is selected
      d <- lp_data(); if (is.null(d) || is.null(d$genes)) 1L else max(d$genes$row) })
    lp_panels_sel <- reactive({ s <- input$lp_panels; if (is.null(s)) c("A","B","C") else s })

    lp_h_px <- reactive({
      sel <- lp_panels_sel(); h <- 0
      if ("A" %in% sel) h <- h + 240
      if ("B" %in% sel) h <- h + 240
      if ("C" %in% sel) h <- h + 120 + 20 * lp_nrows()
      max(280, h + 30)
    })

    # Locus plot image = the continuous-x panels A/B/C (shared position axis via
    # cowplot). The gene x method summary grid (D) is a SEPARATE image below with
    # its own controls + download.
    lp_build <- function() {
      d <- lp_data(); if (is.null(d)) return(NULL)
      sel <- lp_panels_sel()
      tf <- mol_theme_fn(input$lp_theme %||% "bw"); fs <- input$lp_font %||% 14
      gene_h <- max(0.9, 0.3 + 0.16 * lp_nrows())
      cont <- list(); rh <- numeric(0)
      if ("A" %in% sel) { pA <- build_locus_assoc_panel(d, font_size = fs, theme_fn = tf)
                          if (!is.null(pA)) { cont <- c(cont, list(pA)); rh <- c(rh, 1.0) } }
      if ("B" %in% sel) { pB <- build_locus_finemap_panel(d, font_size = fs, theme_fn = tf)
                          if (!is.null(pB)) { cont <- c(cont, list(pB)); rh <- c(rh, 1.0) } }
      if ("C" %in% sel) { pC <- build_locus_gene_panel(d, font_size = fs, theme_fn = tf)
                          if (!is.null(pC)) { cont <- c(cont, list(pC)); rh <- c(rh, gene_h) } }
      stack <- assemble_locus_plot(cont, rel_heights = rh)
      if (is.null(stack) || !("A" %in% sel) || is.null(d$assoc)) return(stack)
      # Horizontal LD r² legend above the association (Manhattan) panel.
      cowplot::plot_grid(lp_ld_legend(fs, horizontal = TRUE), stack,
                         ncol = 1, rel_heights = c(0.05, 1))
    }

    output$lp_plot <- renderPlot({
      p <- lp_build()
      if (!is.null(p)) p else { grid::grid.newpage(); grid::grid.text("No data for this locus.") }
    }, height = function() lp_h_px())

    # ---- Gene association summary grid (panel D) — separate, own controls ----
    lp_grid_opts <- reactive({
      d <- lp_data(); if (is.null(d)) return(list(methods = character(0), expr = character(0), protein = character(0)))
      g <- .pick_gwas("lead_gwas"); req(g)
      tryCatch(lp_grid_options(gwas_data(), g, config_flags(), d),
               error = function(e) list(methods = character(0), expr = character(0), protein = character(0)))
    })
    lpg_panels_sel <- reactive(unique(c(input$lpg_expr, input$lpg_protein)))
    output$lp_grid_controls_ui <- renderUI({
      o <- lp_grid_opts()
      if (length(o$methods) == 0)
        return(div(class = "gd-legend", em("No molecular-association results with a gene in this window.")))
      m_labels <- gsub("\n", " ", o$methods)
      # Default to L=1 variants (L=10 available to add), like the mol-assoc summary.
      m_default <- o$methods[!o$methods %in% c("SuSiE (L=10)", "CALDERA (L=10)")]
      tagList(
        tags$details(class = "gd-details", style = "margin-bottom:6px;",
          tags$summary("Filter data"),
          tags$div(class = "gd-details-body",
            fluidRow(
              column(4, selectInput(ns("lpg_methods"), "Include results from these methods:",
                                    choices = setNames(o$methods, m_labels), selected = m_default, multiple = TRUE)),
              column(4, selectInput(ns("lpg_expr"), "Include expression / splicing panels:",
                                    choices = o$expr, selected = o$expr, multiple = TRUE)),
              column(4, selectInput(ns("lpg_protein"), "Include protein panels:",
                                    choices = o$protein, selected = o$protein, multiple = TRUE))))),
        tags$details(class = "gd-details", style = "margin-bottom:8px;",
          tags$summary("Plot options"),
          tags$div(class = "gd-details-body",
            fluidRow(
              column(4, sliderInput(ns("lpg_font"), "Font size (pt):", min = 8, max = 22, value = 14, step = 1),
                        selectInput(ns("lpg_theme"), "Theme:",
                          choices = c("Black & white"="bw","Minimal"="minimal","Classic"="classic","Light"="light"), selected = "bw")),
              column(4, selectInput(ns("lpg_dl_format"), "Download format:", choices = c("PNG"="png","PDF"="pdf","SVG"="svg")),
                        numericInput(ns("lpg_dl_dpi"), "DPI (PNG):", value = 300, min = 72, max = 600, step = 1)),
              column(4, numericInput(ns("lpg_dl_w"), "Width (in):", value = 12, min = 3, max = 40, step = 0.5),
                        numericInput(ns("lpg_dl_h"), "Height (in):", value = 5, min = 2, max = 60, step = 0.5))),
            downloadButton(ns("lpg_download"), "Download summary", class = "btn btn-default btn-sm")))
      )
    })
    lp_grid_build <- function() {
      d <- lp_data(); if (is.null(d)) return(NULL)
      g <- .pick_gwas("lead_gwas"); req(g)
      tryCatch(build_locus_assoc_grid(gwas_data(), g, config_flags(), d,
               methods = input$lpg_methods, panels = lpg_panels_sel(),
               font_size = input$lpg_font %||% 14, theme_fn = mol_theme_fn(input$lpg_theme %||% "bw")),
               error = function(e) NULL)
    }
    lpg_dims <- reactive({
      d <- lp_data(); if (is.null(d)) return(list(nrows = 6L, ngenes = 10L, nlegends = 2L))
      g <- .pick_gwas("lead_gwas")
      tryCatch(lp_grid_dims(gwas_data(), g, config_flags(), d, input$lpg_methods, lpg_panels_sel()),
               error = function(e) list(nrows = 6L, ngenes = 10L, nlegends = 2L))
    })
    # Size adapts to genes (width) + methods/panels (height), like the mol-assoc plot.
    # Height floors on the legend stack (each colourbar ~150px) so legends never crowd.
    output$lp_grid_ui <- renderUI({
      dm <- lpg_dims()
      w <- max(560, 240 + 26 * dm$ngenes)   # + legends/label margin
      h <- max(220, 90 + 24 * dm$nrows, 60 + 150 * (dm$nlegends %||% 0))
      tags$div(style = "overflow-x: auto;",
        plotOutput(ns("lp_grid"), width = paste0(w, "px"), height = paste0(h, "px")))
    })
    output$lp_grid <- renderPlot({
      p <- lp_grid_build()
      if (!is.null(p)) p else { grid::grid.newpage()
        grid::grid.text("No molecular-association results for the selected methods/panels in this window.") }
    })
    # Legend beneath the grid, mirroring the Molecular Associations summary plot.
    output$lp_grid_legend <- renderUI({
      if (length(lp_grid_opts()$methods) == 0) return(NULL)
      gd_legend(list(
        "Rows / columns" = "Each column is a gene (in genomic-position order); each row is a reference panel, grouped by method in the facet headers.",
        "Colour" = "Association Z-score: blue = the molecular feature is lower with the trait-increasing allele, red = higher, white ≈ no association.",
        "Black circle" = "The association is FDR-significant in that panel and method.",
        "Black square" = "FDR-significant AND supported by colocalisation.",
        "Green" = "Gene highlighted by SuSiE fine-mapping or as the nearest gene to a lead variant.",
        "CALDERA (purple)" = "Effector-gene prioritiser, not an association test (no Z-score): each gene is shaded by its maximum CALDERA causal probability; a black square behind it marks genes above the high-confidence threshold."
      ), heading = "How to read this plot")
    })
    output$lpg_download <- downloadHandler(
      filename = function() sprintf("locus_gene_summary_%s_%s.%s",
                 lp_valid_locus() %||% "locus", input$lp_l %||% "L", input$lpg_dl_format %||% "png"),
      content = function(file) {
        p <- tryCatch(lp_grid_build(), error = function(e) NULL)
        open_plot_device(file, fmt = input$lpg_dl_format %||% "png",
                         w = input$lpg_dl_w %||% 11, h = input$lpg_dl_h %||% 5, units = "in", dpi = input$lpg_dl_dpi %||% 300)
        if (!is.null(p)) print(p) else grid::grid.text("No data to plot.")
        grDevices::dev.off()
      })

    output$lp_legends <- renderUI({
      d <- lp_data(); if (is.null(d)) return(NULL)
      tagList(
        if (!isTRUE(d$has_baseline)) div(class = "gd-legend",
          em("Full-locus baseline unavailable for this bundle - showing credible-set members only (re-run with the current pipeline to add every SNP).")),
        gd_legend(list(
          "Association (top panel)" = "Each SNP's -log10(P); colour = LD r² to the lead (blue low → red high); purple diamond = lead; red line p < 5e-8, blue dashed p < 1e-5; light-blue line (right axis) = recombination rate (cM/Mb). Shaded bands = credible-set spans.",
          "Fine-mapping (middle panel)" = "Per-SNP posterior (PIP); points coloured by credible set; circle = non-coding, triangle = protein-altering variant. Per-set text gives SNP count, purity and coding share.",
          "Genes (bottom panel)" = paste0(
            "Gene bodies coloured by biotype (protein_coding dark; lncRNA / antisense blue / teal; ",
            "pseudogene grey), with a strand-direction arrow. A gold box marks a high-confidence gene; ",
            "the letter badges beneath it are C = CALDERA ≥ 0.5, M = coding variant (alpha ≥ 0.5), ",
            "N = nearest gene, and E / S / P = colocalised (FDR-significant and colocalised) expression / ",
            "splicing / protein.")
        ), heading = "Reading the locus plot"))
    })

    output$lp_download <- downloadHandler(
      filename = function() sprintf("locus_plot_%s_%s.%s",
                 lp_valid_locus() %||% "locus", input$lp_l %||% "L", input$lp_dl_format %||% "png"),
      content = function(file) {
        p <- tryCatch(lp_build(), error = function(e) NULL)
        open_plot_device(file, fmt = input$lp_dl_format %||% "png",
                         w = input$lp_dl_w %||% 10, h = input$lp_dl_h %||% 8, units = "in", dpi = input$lp_dl_dpi %||% 300)
        if (!is.null(p)) print(p) else grid::grid.text("No data to plot.")
        grDevices::dev.off()
      })
  })
}
