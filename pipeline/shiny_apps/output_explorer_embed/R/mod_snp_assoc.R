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
                                  font_size = 12,
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
  gwas_picker <- function(output_id) {
    div(style = "max-width: 260px; margin-bottom: 10px;",
      uiOutput(ns(output_id))
    )
  }

  tabPanel(
    title = "SNP Associations",
    br(),
    p("This tab shows SNP association results. Manhattan plot and per-GWAS ",
      "content include a GWAS selector at the top; Lead variants and ",
      "Fine-mapping are split into Single- and Multi-GWAS sub-tabs so both ",
      "the per-trait and the cross-trait views are always accessible."),
    hr(),
    tabsetPanel(
      id = ns("snp_assoc_tabs"),
      tabPanel(
        title = "Manhattan plot",
        value = "manhattan",
        br(),
        gwas_picker("manhattan_gwas_ui"),
        p("Genome-wide summary of the association signal. Each point is a variant, positioned by chromosome (x) and ", tags$code("-log10(P)"), " (y). Green points are variants clumped with an index (lead) variant; dark-green diamonds are the index variants themselves. Use ", tags$b("Plot options"), " to customise thresholds, gene labels, appearance, and download the figure."),
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
            dataTableOutput(ns("snp_assoc_finemap_table")),
            gd_legend(list(
              "L parameter" = paste0(
                "Maximum number of causal signals SuSiE may fit per locus ",
                "(L1 = single causal variant; L10 = up to ten)."),
              "cs" = "Credible-set identifier: a set of variants that jointly capture one causal signal.",
              "NSNP" = "Number of variants in the credible set (smaller = better resolved).",
              "TopPIP" = "Highest posterior inclusion probability in the set (closer to 1 = more confident).",
              "Gene" = "Nearest gene to the top variant ('None' if unavailable)."
            ), heading = "Column guide"),
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
        choices <- selected_gwas_multi()
        req(length(choices) >= 1)
        # Suppress the picker for single-GWAS bundles.
        if (length(choices) < 2L) return(NULL)
        cur <- isolate(input[[input_id]])
        keep <- if (!is.null(cur) && cur %in% choices) cur else choices[1L]
        selectInput(session$ns(input_id), label,
                     choices = choices, selected = keep, multiple = FALSE)
      })
    }

    # Hide the Multi-GWAS sub-tab under Lead variants when the bundle
    # has only one GWAS — no cross-trait comparison makes sense with 1.
    observe({
      is_multi <- length(selected_gwas_multi()) > 1L
      if (is_multi) {
        showTab("lead_sub_tabs", "multi", session = session)
      } else {
        hideTab("lead_sub_tabs", "multi", session = session)
        cur <- isolate(input$lead_sub_tabs)
        if (!is.null(cur) && cur == "multi") {
          updateTabsetPanel(session, "lead_sub_tabs", selected = "single")
        }
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
          width = 3
        ),
        mainPanel(
          uiOutput(ns("cojo_status_message")),
          uiOutput(ns("lead_status_message")),
          dataTableOutput(ns("snp_assoc_lead_table")),
          uiOutput(ns("lead_legend")),
          br(),
          uiOutput(ns("locus_plot_ui")),
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
    output$manhattan_options_ui <- renderUI({
      if (is.null(manhattan_data())) return(NULL)
      tags$details(class = "gd-details",
        tags$summary("Plot options"),
        tags$div(class = "gd-details-body",
          tags$p(class = "gd-details-intro",
            "Customise thresholds, gene labels, appearance (title, theme, ",
            "font size, point size), and download the figure as PNG / PDF / SVG."),
          fluidRow(
            column(4,
              textInput(ns("mh_title"), "Plot title (optional):", value = ""),
              selectInput(ns("mh_theme"), "Theme:",
                          choices = c("Black & white" = "bw", "Minimal" = "minimal",
                                      "Classic" = "classic", "Light" = "light"),
                          selected = "bw"),
              sliderInput(ns("mh_font_size"), "Font size (pt):",
                          min = 8, max = 20, value = 12, step = 1),
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
              textInput(ns("mh_highlight_threshold"),
                        "Highlight independent leads at P <:", value = "5e-8",
                        placeholder = "max 1e-5 (the clumping P threshold)"),
              checkboxInput(ns("mh_label"),
                            "Label nearest gene at lead variants", value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("mh_label")),
                textInput(ns("mh_label_threshold"),
                          "Label leads with P <:", value = "5e-8")
              )
            ),
            column(4,
              selectInput(ns("mh_dl_format"), "Download format:",
                          choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                          selected = "png"),
              numericInput(ns("mh_dl_width"), "Width (inches):",
                           value = 12, min = 2, max = 40, step = 0.5),
              numericInput(ns("mh_dl_height"), "Height (inches):",
                           value = 5, min = 2, max = 40, step = 0.5),
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
            column(3, colourpicker::colourInput(
              ns("mh_col_member"),   "Colour: clumped members",  value = "#2ca02c")),
            column(3, colourpicker::colourInput(
              ns("mh_col_index"),    "Colour: independent leads", value = "#006400"))
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
        return(tagList(
          tags$div(style = "max-width: 1200px;",
            plotOutput(ns("manhattan_plot"), height = "500px")
          ),
          gd_legend(list(
            "Axes" = "Each point is a variant, positioned by chromosome (x) and -log10(P) (y); higher points are more strongly associated.",
            "Variant filter" = variant_filter_note,
            "Red dashed line" = "Significance line (default P < 5e-8; adjustable in Plot options).",
            "Green points" = "Variants clumped with an index variant (part of the same association signal).",
            "Dark-green diamonds" = "Index (lead) variants below the highlight threshold (adjustable in Plot options).",
            "Gene labels" = "Nearest gene to each labelled index variant (toggle and threshold in Plot options)."
          ), heading = "How to read this plot")
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
        label_threshold     = parse_p(input$mh_label_threshold,     5e-8),
        font_size           = input$mh_font_size  %||% 12,
        point_size          = input$mh_point_size %||% 0.9,
        theme_fn            = mol_theme_fn(input$mh_theme),
        title               = input$mh_title      %||% "",
        col_chr_odd  = input$mh_col_chr_odd  %||% "#8c8c8c",
        col_chr_even = input$mh_col_chr_even %||% "#4d4d4d",
        col_member   = input$mh_col_member   %||% "#2ca02c",
        col_index    = input$mh_col_index    %||% "#006400"
      )
    })

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
          highlight_threshold = parse_p(input$mh_highlight_threshold, 5e-8),
          label_variants      = isTRUE(input$mh_label),
          label_threshold     = parse_p(input$mh_label_threshold,     5e-8),
          font_size           = input$mh_font_size  %||% 12,
          point_size          = input$mh_point_size %||% 0.9,
          theme_fn            = mol_theme_fn(input$mh_theme),
          title               = input$mh_title      %||% "",
          col_chr_odd  = input$mh_col_chr_odd  %||% "#8c8c8c",
          col_chr_even = input$mh_col_chr_even %||% "#4d4d4d",
          col_member   = input$mh_col_member   %||% "#2ca02c",
          col_index    = input$mh_col_index    %||% "#006400"
        )
        w   <- input$mh_dl_width  %||% 12
        h   <- input$mh_dl_height %||% 5
        fmt <- input$mh_dl_format %||% "png"
        switch(fmt,
          png = grDevices::png(file, width = w, height = h, units = "in",
                               res = input$mh_dl_dpi %||% 300),
          pdf = grDevices::pdf(file, width = w, height = h),
          svg = grDevices::svg(file, width = w, height = h)
        )
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
        "NearestGene" = "Gene nearest to the lead variant (not necessarily the causal gene)."
      )
      if (input$clumping_type == "cojo_analysis") {
        items <- c(list(
          "Method (COJO)" = paste0(
            "GCTA-COJO runs a stepwise conditional & joint analysis, keeping variants that stay ",
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
      req(gwas_data(), input$clumping_type)
      g <- .pick_gwas("lead_gwas"); req(g)

      if (input$clumping_type != "ld_clumping") {
        return(div(
          style = "padding: 20px; text-align: center; color: #6c757d; font-style: italic;",
          "Locus plots are available for LD-based clumping results only."
        ))
      }

      locus_plots <- gd_read(gwas_data(), g, "snp_assoc")$locus_plots

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
      df <- snp_assoc_finemap_data()
      req(!finemap_is_empty(df))
      datatable(df, rownames = F, width = 7, options = list(
        columnDefs = list(list(className = 'dt-center', targets = 0:5))))
    })
  })
}
