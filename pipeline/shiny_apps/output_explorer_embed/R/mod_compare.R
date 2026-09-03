# Cross-GWAS comparison sub-modules
#
# Compare-mode content for the Tissue and ATC sub-tabs in the Enrichment
# module. Each entity type has:
#   - a *_compare_ui(ns) function returning the tabPanel body
#   - a *_compare_server(...) function registering outputs, download handlers
#     and reactives
#
# The Enrichment module calls these when comparison_mode() is TRUE. When only
# one GWAS is selected the existing single-GWAS UI is used unchanged.
#
# All data flows through the shared comparison_long reactive; no compare code
# reads the raw package directly.

# Null-coalesce (R 4.4+ has it natively; guard for older versions).
if (!exists("%||%", mode = "function", envir = baseenv(), inherits = FALSE)) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

.gd_hex_alpha <- function(hex, alpha) {
  a <- pmin(1, pmax(0, alpha))
  paste0(hex, sprintf("%02X", round(a * 255)))
}

# Colour palette used by both compare views. Tissue uses a single-hue teal
# ramp with three intensity tiers; ATC-MAGMA uses the same teal ramp;
# ATC-TWAS-GSEA uses signed red/blue with intensity.
.gd_teal <- "#0f766e"
.gd_red  <- "#dc2626"
.gd_blue <- "#2563eb"
.gd_grey <- "#c9d1e1"     # not significant
.gd_hatch_bg <- "#f6f7f9" # not tested background
.gd_hatch_fg <- "#c4cad5" # not tested stroke

# Plot dimensions for a compare heatmap given its row and column counts. Mimics
# the sizing used by calc_plot_dims() for the single-GWAS heatmaps but does not
# require a facet column. Returns pixel dimensions suitable for plotOutput().
.compare_plot_dims <- function(n_rows, n_cols, font_size = 11,
                                 y_labels = NULL,
                                 min_height = 350, min_width = 500) {
  fs_ratio <- font_size / 11
  y_label_px <- if (!is.null(y_labels) && length(y_labels) > 0) {
    max(nchar(as.character(y_labels)), na.rm = TRUE) * font_size * 0.55 + 20
  } else 160
  panel_w <- max(30, 40 * fs_ratio) * max(1L, n_cols)
  panel_h <- max(18, 22 * fs_ratio) * max(1L, n_rows)
  x_label_h <- 6 * font_size + 30
  legend_h <- 70
  height <- max(min_height, x_label_h + panel_h + legend_h + 40)
  width  <- max(min_width, y_label_px + panel_w + 180)
  list(height = round(height), width = round(width))
}

# Add the "outline for nominal-sig" (larger solid black circle behind the
# data circle) and "black square for FDR-sig" (larger solid black square)
# overlays, then re-draw the data circles on top so the fill shows through.
# Matches build_tx_atc_gtable / build_tx_drug_gtable in the single-GWAS
# code path. The base data layer must already be added to `gg`.
.compare_add_sig_overlays <- function(gg, data, base_layer,
                                        nominal_flag = "nom_sig",
                                        fdr_flag = "fdr_sig",
                                        point_size = 4) {
  nom <- data[data[[nominal_flag]] %in% TRUE, ]
  fdr <- data[data[[fdr_flag]]     %in% TRUE, ]
  if (nrow(nom) > 0) {
    gg <- gg + ggplot2::geom_point(data = nom, colour = "black", fill = NA,
                                     size = point_size + 1)
  }
  if (nrow(fdr) > 0) {
    gg <- gg + ggplot2::geom_point(data = fdr, colour = "black", fill = NA,
                                     size = point_size + 2, shape = 15)
  }
  # Re-draw the data layer so the fill sits on top of the overlay markers.
  gg + base_layer
}

########################################
# TISSUE COMPARE
########################################

tissue_compare_ui <- function(ns) {
  tagList(
    br(),
    p(
      "Cross-GWAS view of MAGMA tissue-specific enrichment across GTEx v8 tissues. ",
      "Cells are coloured in three tiers: retained after conditional analysis, ",
      "FDR-significant but not retained, and nominal-significant only. Rows are ",
      "ordered by recurrence (how many GWAS reach the chosen significance basis) ",
      "then by minimum p-value."
    ),
    hr(),
    tags$details(class = "gd-details",
      tags$summary("Filter data"),
      tags$div(class = "gd-details-body",
        fluidRow(
          column(3,
            radioButtons(ns("sig_basis"), "Significance basis:",
                          choices = c("FDR" = "fdr", "P" = "p"),
                          selected = "fdr", inline = TRUE),
            numericInput(ns("sig_threshold"), "Significance threshold:",
                          value = 0.05, min = 1e-12, max = 1, step = 0.01)
          ),
          column(3,
            checkboxInput(ns("only_recurrent"),
                          "Only show tissues significant in ≥ k GWAS",
                          value = FALSE),
            conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("only_recurrent")),
              sliderInput(ns("k_min"), "k:",
                           min = 1, max = 9, value = 2, step = 1)
            )
          ),
          column(3,
            radioButtons(ns("cell_metric"), "Cell colour:",
                         choices = c("-log10(FDR)" = "fdr",
                                      "-log10(P)"   = "p"),
                         selected = "fdr", inline = TRUE),
            selectInput(ns("gwas_sort"), "Order GWAS by:",
                         choices = c("As selected"   = "as_selected",
                                      "Alphabetical" = "alphabetical",
                                      "Sample size"  = "n",
                                      "N sig SNPs"   = "n_sig_snp",
                                      "SNP-h²"       = "h2"),
                         selected = "as_selected"),
            sliderInput(ns("plot_font_size"), "Font size (pt):",
                         min = 8, max = 20, value = 12, step = 1),
            sliderInput(ns("plot_point_size"), "Point size:",
                         min = 2, max = 10, value = 4, step = 1)
          ),
          column(3,
            selectInput(ns("dl_format"), "Download format:",
                        choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                        selected = "png"),
            numericInput(ns("dl_width"), "Width (in):",
                         value = 10, min = 4, max = 24, step = 0.5),
            numericInput(ns("dl_height"), "Height (in):",
                         value = 10, min = 4, max = 24, step = 0.5),
            downloadButton(ns("download_plot"), "Download plot"),
            downloadButton(ns("download_csv"), "Download matrix CSV")
          )
        )
      )
    ),
    br(),
    tags$div(style = "max-width: 1100px; overflow-x: auto;",
      uiOutput(ns("tissue_compare_plot_ui"))
    ),
    gd_legend(list(
      "Retained (dark teal)"   = "FDR-significant and kept after the conditional analysis.",
      "FDR-sig (medium teal)"  = "FDR-significant but not retained after conditioning.",
      "Nominal-sig (pale teal)"= "Nominal (P < 0.05) but not FDR-significant.",
      "Not significant (grey)" = "Tissue was tested but did not reach nominal significance.",
      "Not tested (hatched)"   = "Tissue was not present in this GWAS's results."
    ), heading = "How to read this heatmap"),
    br(),
    tags$div(style = "max-width: 1100px;",
      h4("Underlying data"),
      DT::DTOutput(ns("tissue_compare_tbl"))
    )
  )
}

# Long -> tidy plotting frame with tier categorisation
.tissue_compare_frame <- function(long, gwas_vec, only_recurrent, k_min,
                                    sig_basis = "fdr", sig_threshold = 0.05) {
  # Slice to tissue rows and to the selected GWAS
  slice <- long[method == "MAGMA-tissue" & gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)

  # Recurrence uses the user-chosen basis + threshold; row order is by
  # recurrence desc then min p asc.
  basis_col <- if (identical(sig_basis, "p")) "p" else "fdr"
  rec <- slice[, .(
    k         = sum(!is.na(.SD[[1L]]) & .SD[[1L]] < sig_threshold),
    n_fdr_sig = sum(!is.na(fdr) & fdr < 0.05),
    n_nom_sig = sum(!is.na(p)   & p   < 0.05),
    min_p     = suppressWarnings(min(p, na.rm = TRUE))
  ), by = entity_id, .SDcols = basis_col]
  rec <- rec[order(-k, min_p)]

  if (isTRUE(only_recurrent)) {
    rec <- rec[k >= k_min]
  }

  if (nrow(rec) == 0) return(NULL)
  slice <- slice[entity_id %in% rec$entity_id]
  slice[, entity_id := factor(entity_id, levels = rev(rec$entity_id))]
  slice[, gwas := factor(gwas, levels = gwas_vec)]

  # Tier fill uses standard fixed thresholds (P<0.05, FDR<0.05) so the
  # legend is unambiguous regardless of the user's recurrence basis.
  slice[, tier := "Not significant"]
  slice[!is.na(p) & p < 0.05,        tier := "Nominal-sig"]
  slice[!is.na(fdr) & fdr < 0.05,    tier := "FDR-sig"]
  slice[!is.na(evidence) & evidence, tier := "Retained"]
  slice[, tier := factor(tier, levels = c("Not significant", "Nominal-sig",
                                            "FDR-sig", "Retained"))]
  list(slice = slice, rec = rec)
}

.tissue_compare_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                    metric = "fdr",
                                    sig_basis = "fdr",
                                    sig_threshold = 0.05,
                                    font_size = 12,
                                    point_size = 4) {
  packed <- .tissue_compare_frame(long, gwas_vec, only_recurrent, k_min,
                                    sig_basis, sig_threshold)
  if (is.null(packed)) return(NULL)
  slice <- packed$slice

  # Point-style heatmap: filled coloured circle by -log10(FDR) (or -log10(P)
  # when metric = "p"), with a solid-black-circle overlay for nominal-sig,
  # a solid-black-square overlay for FDR-sig, and the data circle re-drawn
  # on top so the fill sits inside the overlay markers. Retained tissues
  # get an inner black dot (shape 19) so they are distinguishable from
  # FDR-sig-but-not-retained.
  slice[, minus_log10 := pmin(-log10(if (identical(metric, "p")) p else fdr), 12)]
  slice[, nom_sig := !is.na(p) & p < 0.05]
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]
  slice[, retained := !is.na(evidence) & evidence]

  scale_lab <- if (identical(metric, "p")) "-log10(P)" else "-log10(FDR)"

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = minus_log10), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(slice, ggplot2::aes(x = gwas, y = entity_id)) +
    base_layer +
    ggplot2::scale_fill_gradient(low = "#e6f4f1", high = .gd_teal,
                                   name = scale_lab,
                                   na.value = .gd_grey,
                                   limits = c(0, 12))

  gg <- .compare_add_sig_overlays(gg, slice, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)
  # Retained tissues get a small solid black dot inside the circle.
  ret <- slice[slice$retained, ]
  if (nrow(ret) > 0) {
    gg <- gg + ggplot2::geom_point(data = ret, shape = 19,
                                     colour = "black", size = point_size / 2)
  }

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = c(0.5, 0.5)) +
    ggplot2::scale_y_discrete(expand = c(0.5, 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL,
                   caption = "Ring = nominal-sig (P<0.05). Black square = FDR-sig. Inner dot = retained after conditional.") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

tissue_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                    comparison_long) {
  moduleServer(id, function(input, output, session) {

    # Cap the k slider at the number of selected GWAS.
    observeEvent(selected_gwas_multi(), {
      n_sel <- length(selected_gwas_multi())
      if (n_sel < 1) return()
      cur <- if (is.null(input$k_min)) 2L else as.integer(input$k_min)
      updateSliderInput(session, "k_min",
                        max = max(n_sel, 1L),
                        value = min(cur, n_sel))
    })

    gwas_vec_r <- reactive({
      order_gwas(selected_gwas_multi(),
                  if (is.null(input$gwas_sort)) "as_selected" else input$gwas_sort,
                  NULL)
    })

    plot_obj <- reactive({
      req(comparison_long())
      .tissue_compare_ggplot(
        long           = comparison_long(),
        gwas_vec       = gwas_vec_r(),
        only_recurrent = isTRUE(input$only_recurrent),
        k_min          = if (is.null(input$k_min)) 2L else as.integer(input$k_min),
        metric         = if (is.null(input$cell_metric)) "fdr" else input$cell_metric,
        sig_basis      = if (is.null(input$sig_basis)) "fdr" else input$sig_basis,
        sig_threshold  = if (is.null(input$sig_threshold)) 0.05 else as.numeric(input$sig_threshold),
        font_size      = if (is.null(input$plot_font_size)) 12 else as.numeric(input$plot_font_size),
        point_size     = if (is.null(input$plot_point_size)) 4 else as.numeric(input$plot_point_size)
      )
    })

    plot_dims <- reactive({
      req(comparison_long())
      slice <- comparison_long()[method == "MAGMA-tissue" & gwas %in% gwas_vec_r()]
      .compare_plot_dims(
        n_rows    = length(unique(slice$entity_id)),
        n_cols    = length(gwas_vec_r()),
        font_size = if (is.null(input$plot_font_size)) 12 else as.numeric(input$plot_font_size),
        y_labels  = unique(slice$entity_id)
      )
    })

    output$tissue_compare_plot_ui <- renderUI({
      dims <- plot_dims()
      plotOutput(session$ns("tissue_compare_plot"),
                  height = dims$height, width = dims$width)
    })

    output$tissue_compare_plot <- renderPlot({
      p <- plot_obj()
      if (is.null(p)) {
        plot.new(); title("No tissues meet the current filter")
        return(invisible())
      }
      print(p)
    })

    output$tissue_compare_tbl <- DT::renderDT({
      req(comparison_long())
      gwas_vec <- gwas_vec_r()
      slice <- comparison_long()[method == "MAGMA-tissue" & gwas %in% gwas_vec]
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS     = factor(slice$gwas, levels = gwas_vec),
        Tissue   = slice$entity_id,
        BETA     = signif(slice$statistic, 3),
        SE       = signif(slice$se, 3),
        P        = signif(slice$p, 3),
        `P.FDR`  = signif(slice$fdr, 3),
        Retained = ifelse(is.na(slice$evidence), "—",
                          ifelse(slice$evidence, "Yes", "No")),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE,
                    filter = "top",
                    options = list(pageLength = 20, server = TRUE,
                                    order = list(list(4, "asc"))))
    }, server = TRUE)

    output$download_plot <- downloadHandler(
      filename = function() sprintf("tissue_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$dl_format),
      content = function(file) {
        p <- plot_obj()
        if (is.null(p)) {
          grDevices::png(file, width = 400, height = 200); dev.off(); return()
        }
        fmt <- input$dl_format
        w <- input$dl_width; h <- input$dl_height
        if (fmt == "png") {
          grDevices::png(file, width = w, height = h, units = "in", res = 300)
        } else if (fmt == "pdf") {
          grDevices::pdf(file, width = w, height = h)
        } else {
          grDevices::svg(file, width = w, height = h)
        }
        print(p)
        grDevices::dev.off()
      }
    )

    output$download_csv <- downloadHandler(
      filename = function() sprintf("tissue_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- gwas_vec_r()
        long <- comparison_long()[method == "MAGMA-tissue" & gwas %in% gwas_vec]
        wide <- pivot_matrix(long, "fdr", gwas_vec)
        ret  <- pivot_matrix(long, "evidence", gwas_vec)
        display <- ifelse(
          is.na(wide), "NT",
          ifelse(!is.na(ret) & ret == 1,
                 sprintf("%.2e (R)", wide),
                 sprintf("%.2e", wide))
        )
        display <- matrix(display, nrow = nrow(wide),
                          dimnames = dimnames(wide))
        header <- sprintf("# GenoDisc tissue compare CSV | %s | sig_basis=%s threshold=%g k_min=%s",
                          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                          input$sig_basis %||% "fdr",
                          input$sig_threshold %||% 0.05,
                          input$k_min %||% 2L)
        con <- file(file, "w")
        on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: P.FDR (R = retained after conditional). NT = not tested.", con)
        out <- data.frame(Tissue = rownames(display), display,
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )
  })
}

########################################
# LOCUS COMPARE
########################################

.locus_methods <- c("clump", "COJO")

locus_compare_ui <- function(ns) {
  tagList(
    br(),
    p(
      "Cross-GWAS locus-level associations. Loci are keyed by their nearest ",
      "gene (from the clump / COJO table). Each cell shows the smallest p-value ",
      "at that locus for the selected GWAS. Blank cell = the locus has no ",
      "clumped / independent signal in that GWAS. Default significance threshold ",
      "is P < 5×10⁻⁸ (genome-wide significance)."
    ),
    hr(),
    tags$details(class = "gd-details",
      tags$summary("Filter data"),
      tags$div(class = "gd-details-body",
        fluidRow(
          column(3,
            selectInput(ns("method"), "Method:",
                        choices = .locus_methods,
                        selected = "clump")
          ),
          column(3,
            radioButtons(ns("sig_basis"), "Significance basis:",
                          choices = c("P" = "p"),
                          selected = "p", inline = TRUE),
            numericInput(ns("sig_threshold"), "Significance threshold:",
                          value = 5e-8, min = 1e-20, max = 1, step = 1e-8)
          ),
          column(3,
            checkboxInput(ns("only_recurrent"),
                          "Only show loci significant in ≥ k GWAS",
                          value = TRUE),
            conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("only_recurrent")),
              sliderInput(ns("k_min"), "k:",
                           min = 1, max = 9, value = 2, step = 1)
            ),
            numericInput(ns("row_cap"), "Rows shown in heatmap:",
                          value = 50, min = 5, max = 500, step = 5)
          ),
          column(3,
            selectInput(ns("gwas_sort"), "Order GWAS by:",
                         choices = .compare_gwas_sort_choices,
                         selected = "as_selected"),
            sliderInput(ns("plot_font_size"), "Font size (pt):",
                         min = 8, max = 20, value = 11, step = 1),
            sliderInput(ns("plot_point_size"), "Point size:",
                         min = 2, max = 10, value = 4, step = 1)
          )
        ),
        fluidRow(
          .dl_and_download_column(ns, "locus", default_h = 12)
        )
      )
    ),
    br(),
    tags$div(style = "max-width: 1100px; overflow-x: auto;",
      uiOutput(ns("locus_compare_plot_ui"))
    ),
    br(),
    tags$div(style = "max-width: 1100px;",
      h4("Underlying data"),
      DT::DTOutput(ns("locus_compare_tbl"))
    )
  )
}

.locus_compare_ggplot <- function(long, gwas_vec, method_pick,
                                     only_recurrent, k_min,
                                     sig_threshold = 5e-8,
                                     row_cap = 50, font_size = 11,
                                     point_size = 4) {
  slice <- long[method == method_pick & gwas %in% gwas_vec &
                    entity_type == "locus"]
  if (nrow(slice) == 0) return(NULL)
  slice <- pick_best_per_cell(slice, c("gwas", "entity_id"))

  rec <- slice[, .(
    k       = sum(!is.na(p) & p < sig_threshold),
    min_val = suppressWarnings(min(p, na.rm = TRUE))
  ), by = entity_id]
  rec <- rec[order(-k, min_val)]

  if (isTRUE(only_recurrent)) rec <- rec[k >= k_min]
  if (nrow(rec) == 0) return(NULL)

  cap <- max(1L, min(as.integer(row_cap), nrow(rec)))
  rec <- rec[seq_len(cap)]
  slice <- slice[entity_id %in% rec$entity_id]

  # Loci use a wide p-value range; clamp at 40 (P ~ 1e-40) which is generous
  # for most GWAS. Anything beyond that saturates at the top of the ramp.
  slice[, minus_log10 := pmin(-log10(p), 40)]
  slice[, nom_sig := !is.na(p) & p < 5e-8]   # GWS threshold as the "hit" mark
  slice[, fdr_sig := !is.na(p) & p < sig_threshold]
  slice[, entity_id := factor(entity_id, levels = rev(rec$entity_id))]
  slice[, gwas := factor(gwas, levels = gwas_vec)]

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = minus_log10), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(slice, ggplot2::aes(x = gwas, y = entity_id)) +
    base_layer +
    ggplot2::scale_fill_gradient(low = "#e6f4f1", high = .gd_teal,
                                   name = "-log10(P)",
                                   na.value = .gd_grey, limits = c(0, 40))

  gg <- .compare_add_sig_overlays(gg, slice, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = c(0.5, 0.5)) +
    ggplot2::scale_y_discrete(expand = c(0.5, 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL,
                   caption = sprintf(
                     "Ring = GWS (P<5x10^-8). Black square = P<%.0e. Blank cell = no clumped / independent signal for this GWAS.",
                     sig_threshold)) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

locus_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                    comparison_long) {
  moduleServer(id, function(input, output, session) {

    observeEvent(selected_gwas_multi(), {
      n_sel <- length(selected_gwas_multi())
      if (n_sel < 1) return()
      cur <- if (is.null(input$k_min)) 2L else as.integer(input$k_min)
      updateSliderInput(session, "k_min",
                        max = max(n_sel, 1L),
                        value = min(cur, n_sel))
    })

    gwas_vec_r <- reactive({
      order_gwas(selected_gwas_multi(),
                  if (is.null(input$gwas_sort)) "as_selected" else input$gwas_sort,
                  NULL)
    })

    plot_obj <- reactive({
      req(comparison_long())
      .locus_compare_ggplot(
        long           = comparison_long(),
        gwas_vec       = gwas_vec_r(),
        method_pick    = input$method %||% "clump",
        only_recurrent = isTRUE(input$only_recurrent),
        k_min          = if (is.null(input$k_min)) 2L else as.integer(input$k_min),
        sig_threshold  = if (is.null(input$sig_threshold)) 5e-8 else as.numeric(input$sig_threshold),
        row_cap        = if (is.null(input$row_cap)) 50L else as.integer(input$row_cap),
        font_size      = if (is.null(input$plot_font_size)) 11 else as.numeric(input$plot_font_size),
        point_size     = if (is.null(input$plot_point_size)) 4 else as.numeric(input$plot_point_size)
      )
    })

    plot_dims <- reactive({
      p <- plot_obj()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_id))
      .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = length(gwas_vec_r()),
        font_size = if (is.null(input$plot_font_size)) 11 else as.numeric(input$plot_font_size),
        y_labels  = if (is.null(p)) NULL else as.character(unique(p$data$entity_id))
      )
    })

    output$locus_compare_plot_ui <- renderUI({
      dims <- plot_dims()
      plotOutput(session$ns("locus_compare_plot"),
                  height = dims$height, width = dims$width)
    })

    output$locus_compare_plot <- renderPlot({
      p <- plot_obj()
      if (is.null(p)) {
        plot.new(); title("No loci meet the current filter")
        return(invisible())
      }
      print(p)
    })

    output$locus_compare_tbl <- DT::renderDT({
      req(comparison_long())
      gwas_vec <- gwas_vec_r()
      picked_method <- input$method %||% "clump"
      slice <- comparison_long()[method == picked_method & gwas %in% gwas_vec &
                                    entity_type == "locus"]
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS       = factor(slice$gwas, levels = gwas_vec),
        Locus      = slice$entity_id,
        BETA       = signif(slice$statistic, 3),
        SE         = signif(slice$se, 3),
        P          = signif(slice$p, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(4, "asc"))))
    }, server = TRUE)

    output$locus_download_plot <- downloadHandler(
      filename = function() sprintf("locus_compare_%s_%s.%s",
                                     input$method %||% "clump",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$locus_dl_format),
      content = function(file) {
        p <- plot_obj(); if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$locus_dl_format
        w <- input$locus_dl_width; h <- input$locus_dl_height
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$locus_download_csv <- downloadHandler(
      filename = function() sprintf("locus_compare_%s_matrix_%s.csv",
                                     input$method %||% "clump",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- gwas_vec_r()
        picked_method <- input$method %||% "clump"
        long <- comparison_long()[method == picked_method & gwas %in% gwas_vec &
                                    entity_type == "locus"]
        long <- pick_best_per_cell(long, c("gwas", "entity_id"))
        wide <- pivot_matrix(long, "p", gwas_vec)
        display <- ifelse(is.na(wide), "",
                          format(wide, scientific = TRUE, digits = 3))
        display <- matrix(display, nrow = nrow(wide), dimnames = dimnames(wide))
        header <- sprintf("# GenoDisc locus compare CSV | %s | method=%s threshold=%g k_min=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           picked_method,
                           input$sig_threshold %||% 5e-8,
                           input$k_min %||% 2L)
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: P (best per locus x GWAS). Empty = no clumped / independent signal for that GWAS.", con)
        out <- data.frame(Locus = rownames(display), display,
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )
  })
}

########################################
# ATC COMPARE
########################################

.compare_gwas_sort_choices <- c(
  "As selected"  = "as_selected",
  "Alphabetical" = "alphabetical",
  "Sample size"  = "n",
  "N sig SNPs"   = "n_sig_snp",
  "SNP-h²"       = "h2"
)

.dl_and_download_column <- function(ns, prefix, default_w = 10, default_h = 12) {
  column(3,
    selectInput(ns(paste0(prefix, "_dl_format")), "Download format:",
                choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                selected = "png"),
    numericInput(ns(paste0(prefix, "_dl_width")), "Width (in):",
                 value = default_w, min = 4, max = 24, step = 0.5),
    numericInput(ns(paste0(prefix, "_dl_height")), "Height (in):",
                 value = default_h, min = 4, max = 32, step = 0.5),
    downloadButton(ns(paste0(prefix, "_download_plot")), "Download plot"),
    downloadButton(ns(paste0(prefix, "_download_csv")), "Download matrix CSV")
  )
}

########################################
# GENE COMPARE
########################################

.gene_methods <- c("MAGMA-gene", "TWAS-FUSION", "SMR-expression",
                    "PWAS-FUSION", "SMR-protein")

# Methods where "evidence" is meaningful (colocalisation / HEIDI).
.gene_evidence_methods <- c("TWAS-FUSION", "SMR-expression",
                             "PWAS-FUSION", "SMR-protein")

gene_compare_ui <- function(ns) {
  tagList(
    br(),
    p(
      "Cross-GWAS gene-level associations. Pick a method to compare a ",
      "gene x GWAS matrix; where a method has multiple panels, the default ",
      "reduces to the smallest p-value per (gene, GWAS) cell. Rows are ",
      "ordered by recurrence (how many GWAS meet the chosen significance ",
      "basis) then by minimum p-value; the top rows are shown in the ",
      "heatmap and the full slice appears in the table below."
    ),
    hr(),
    tags$details(class = "gd-details",
      tags$summary("Filter data"),
      tags$div(class = "gd-details-body",
        fluidRow(
          column(3,
            selectInput(ns("method"), "Method:",
                        choices = .gene_methods,
                        selected = "MAGMA-gene"),
            selectInput(ns("panel"), "Panel:",
                        choices = c("Best-per-cell (min P)" = "__best__"),
                        selected = "__best__"),
            checkboxInput(ns("evidence_required"),
                          "Require colocalisation / HEIDI evidence",
                          value = FALSE)
          ),
          column(3,
            radioButtons(ns("sig_basis"), "Significance basis:",
                          choices = c("FDR" = "fdr", "P" = "p"),
                          selected = "fdr", inline = TRUE),
            numericInput(ns("sig_threshold"), "Significance threshold:",
                          value = 0.05, min = 1e-12, max = 1, step = 0.01)
          ),
          column(3,
            checkboxInput(ns("only_recurrent"),
                          "Only show genes significant in ≥ k GWAS",
                          value = TRUE),
            conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("only_recurrent")),
              sliderInput(ns("k_min"), "k:",
                           min = 1, max = 9, value = 2, step = 1)
            ),
            numericInput(ns("row_cap"), "Rows shown in heatmap:",
                          value = 50, min = 5, max = 500, step = 5)
          ),
          column(3,
            selectInput(ns("gwas_sort"), "Order GWAS by:",
                         choices = .compare_gwas_sort_choices,
                         selected = "as_selected"),
            sliderInput(ns("plot_font_size"), "Font size (pt):",
                         min = 8, max = 20, value = 11, step = 1),
            sliderInput(ns("plot_point_size"), "Point size:",
                         min = 2, max = 10, value = 4, step = 1)
          )
        ),
        fluidRow(
          .dl_and_download_column(ns, "gene", default_h = 12)
        )
      )
    ),
    br(),
    tags$div(style = "max-width: 1100px; overflow-x: auto;",
      uiOutput(ns("gene_compare_plot_ui"))
    ),
    br(),
    tags$div(style = "max-width: 1100px;",
      h4("Underlying data"),
      DT::DTOutput(ns("gene_compare_tbl"))
    )
  )
}

.gene_compare_ggplot <- function(long, gwas_vec, method_pick, panel_pick,
                                    only_recurrent, k_min,
                                    sig_basis = "fdr", sig_threshold = 0.05,
                                    evidence_required = FALSE,
                                    row_cap = 50, font_size = 11,
                                    point_size = 4) {
  slice <- long[method == method_pick & gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)

  # panel reduction (best-per-cell or single panel)
  if (!identical(panel_pick, "__best__") && !is.null(panel_pick) &&
      nzchar(panel_pick)) {
    slice <- slice[panel == panel_pick]
  } else {
    slice <- pick_best_per_cell(slice, c("gwas", "entity_id"))
  }
  if (nrow(slice) == 0) return(NULL)

  # evidence-required only meaningful for methods that have evidence
  if (isTRUE(evidence_required) && method_pick %in% .gene_evidence_methods) {
    slice <- slice[!is.na(evidence) & evidence == TRUE]
  }
  if (nrow(slice) == 0) return(NULL)

  # recurrence ordering
  basis_col <- if (identical(sig_basis, "p")) "p" else "fdr"
  rec <- slice[, .(
    k       = sum(!is.na(.SD[[1L]]) & .SD[[1L]] < sig_threshold),
    min_val = suppressWarnings(min(.SD[[1L]], na.rm = TRUE))
  ), by = entity_id, .SDcols = basis_col]
  rec <- rec[order(-k, min_val)]

  if (isTRUE(only_recurrent)) rec <- rec[k >= k_min]
  if (nrow(rec) == 0) return(NULL)

  # Cap heatmap rows
  cap <- max(1L, min(as.integer(row_cap), nrow(rec)))
  rec <- rec[seq_len(cap)]
  slice <- slice[entity_id %in% rec$entity_id]

  # Fill: -log10(FDR) or -log10(P) depending on user's basis; MAGMA/
  # unsigned methods -> single-hue teal. Direction is shown in the DT
  # table for methods that carry it.
  slice[, minus_log10 := pmin(-log10(.SD[[1L]]), 12), .SDcols = basis_col]
  slice[, nom_sig := !is.na(p)   & p   < 0.05]
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]
  slice[, entity_id := factor(entity_id, levels = rev(rec$entity_id))]
  slice[, gwas := factor(gwas, levels = gwas_vec)]

  scale_lab <- if (identical(sig_basis, "p")) "-log10(P)" else "-log10(FDR)"

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = minus_log10), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(slice, ggplot2::aes(x = gwas, y = entity_id)) +
    base_layer +
    ggplot2::scale_fill_gradient(low = "#e6f4f1", high = .gd_teal,
                                   name = scale_lab,
                                   na.value = .gd_grey, limits = c(0, 12))

  gg <- .compare_add_sig_overlays(gg, slice, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = c(0.5, 0.5)) +
    ggplot2::scale_y_discrete(expand = c(0.5, 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL,
                   caption = "Ring = nominal-sig. Black square = FDR-sig.") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

gene_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                  comparison_long) {
  moduleServer(id, function(input, output, session) {

    # Cap k slider at number of selected GWAS
    observeEvent(selected_gwas_multi(), {
      n_sel <- length(selected_gwas_multi())
      if (n_sel < 1) return()
      cur <- if (is.null(input$k_min)) 2L else as.integer(input$k_min)
      updateSliderInput(session, "k_min",
                        max = max(n_sel, 1L),
                        value = min(cur, n_sel))
    })

    # When method changes, refresh panel dropdown to only include panels
    # present for that method. MAGMA-gene has no panels.
    observeEvent(list(comparison_long(), input$method), {
      long <- comparison_long()
      picked_method <- input$method %||% "MAGMA-gene"
      panels <- unique(long[method == picked_method, panel])
      panels <- panels[!is.na(panels)]
      if (length(panels) == 0) {
        choices <- c("N/A (single-panel method)" = "__best__")
      } else {
        choices <- c("Best-per-cell (min P)" = "__best__",
                      setNames(panels, panels))
      }
      updateSelectInput(session, "panel", choices = choices, selected = "__best__")
    })

    gwas_vec_r <- reactive({
      order_gwas(selected_gwas_multi(),
                  if (is.null(input$gwas_sort)) "as_selected" else input$gwas_sort,
                  NULL)
    })

    plot_obj <- reactive({
      req(comparison_long())
      .gene_compare_ggplot(
        long              = comparison_long(),
        gwas_vec          = gwas_vec_r(),
        method_pick       = input$method %||% "MAGMA-gene",
        panel_pick        = input$panel %||% "__best__",
        only_recurrent    = isTRUE(input$only_recurrent),
        k_min             = if (is.null(input$k_min)) 2L else as.integer(input$k_min),
        sig_basis         = if (is.null(input$sig_basis)) "fdr" else input$sig_basis,
        sig_threshold     = if (is.null(input$sig_threshold)) 0.05 else as.numeric(input$sig_threshold),
        evidence_required = isTRUE(input$evidence_required),
        row_cap           = if (is.null(input$row_cap)) 50L else as.integer(input$row_cap),
        font_size         = if (is.null(input$plot_font_size)) 11 else as.numeric(input$plot_font_size),
        point_size        = if (is.null(input$plot_point_size)) 4 else as.numeric(input$plot_point_size)
      )
    })

    plot_dims <- reactive({
      p <- plot_obj()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_id))
      .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = length(gwas_vec_r()),
        font_size = if (is.null(input$plot_font_size)) 11 else as.numeric(input$plot_font_size),
        y_labels  = if (is.null(p)) NULL else as.character(unique(p$data$entity_id))
      )
    })

    output$gene_compare_plot_ui <- renderUI({
      dims <- plot_dims()
      plotOutput(session$ns("gene_compare_plot"),
                  height = dims$height, width = dims$width)
    })

    output$gene_compare_plot <- renderPlot({
      p <- plot_obj()
      if (is.null(p)) {
        plot.new(); title("No genes meet the current filter")
        return(invisible())
      }
      print(p)
    })

    output$gene_compare_tbl <- DT::renderDT({
      req(comparison_long())
      gwas_vec <- gwas_vec_r()
      picked_method <- input$method %||% "MAGMA-gene"
      slice <- comparison_long()[method == picked_method & gwas %in% gwas_vec]
      picked_panel <- input$panel %||% "__best__"
      if (!identical(picked_panel, "__best__") && nzchar(picked_panel)) {
        slice <- slice[panel == picked_panel]
      }
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS      = factor(slice$gwas, levels = gwas_vec),
        Gene      = slice$entity_id,
        Panel     = slice$panel,
        Statistic = signif(slice$statistic, 3),
        SE        = signif(slice$se, 3),
        P         = signif(slice$p, 3),
        `P.FDR`   = signif(slice$fdr, 3),
        Evidence  = ifelse(is.na(slice$evidence), "—",
                            ifelse(slice$evidence, "Yes", "No")),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(6, "asc"))))
    }, server = TRUE)

    output$gene_download_plot <- downloadHandler(
      filename = function() sprintf("gene_compare_%s_%s.%s",
                                     gsub("[^A-Za-z0-9]+", "_", input$method %||% "MAGMA-gene"),
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$gene_dl_format),
      content = function(file) {
        p <- plot_obj(); if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$gene_dl_format
        w <- input$gene_dl_width; h <- input$gene_dl_height
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$gene_download_csv <- downloadHandler(
      filename = function() sprintf("gene_compare_%s_matrix_%s.csv",
                                     gsub("[^A-Za-z0-9]+", "_", input$method %||% "MAGMA-gene"),
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- gwas_vec_r()
        picked_method <- input$method %||% "MAGMA-gene"
        long <- comparison_long()[method == picked_method & gwas %in% gwas_vec]
        picked_panel <- input$panel %||% "__best__"
        if (!identical(picked_panel, "__best__") && nzchar(picked_panel)) {
          long <- long[panel == picked_panel]
        } else {
          long <- pick_best_per_cell(long, c("gwas", "entity_id"))
        }
        wide <- pivot_matrix(long, "fdr", gwas_vec)
        display <- ifelse(is.na(wide), "",
                          format(wide, scientific = TRUE, digits = 3))
        display <- matrix(display, nrow = nrow(wide), dimnames = dimnames(wide))
        header <- sprintf("# GenoDisc gene compare CSV | %s | method=%s panel=%s sig_basis=%s threshold=%g k_min=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           picked_method, picked_panel,
                           input$sig_basis %||% "fdr",
                           input$sig_threshold %||% 0.05,
                           input$k_min %||% 2L)
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: P.FDR. Empty cell = gene not present in that GWAS's method output.", con)
        out <- data.frame(Gene = rownames(display), display,
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )
  })
}

atc_compare_ui <- function(ns) {
  tabsetPanel(
    tabPanel("MAGMA", br(),
      p("Cross-GWAS MAGMA drug-class enrichment. Cell colour intensity is ",
        "-log10(FDR); cells with FDR < 0.05 are outlined in black."),
      hr(),
      tags$details(class = "gd-details",
        tags$summary("Filter data"),
        tags$div(class = "gd-details-body",
          fluidRow(
            column(3,
              radioButtons(ns("magma_sig_basis"), "Significance basis:",
                            choices = c("FDR" = "fdr", "P" = "p"),
                            selected = "fdr", inline = TRUE),
              numericInput(ns("magma_sig_threshold"), "Significance threshold:",
                            value = 0.05, min = 1e-12, max = 1, step = 0.01)
            ),
            column(3,
              checkboxInput(ns("magma_only_recurrent"),
                            "Only show classes significant in ≥ k GWAS",
                            value = TRUE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("magma_only_recurrent")),
                sliderInput(ns("magma_k_min"), "k:",
                             min = 1, max = 9, value = 2, step = 1)
              )
            ),
            column(3,
              selectInput(ns("magma_gwas_sort"), "Order GWAS by:",
                           choices = .compare_gwas_sort_choices,
                           selected = "as_selected"),
              sliderInput(ns("magma_plot_font_size"), "Font size (pt):",
                           min = 8, max = 20, value = 11, step = 1),
              sliderInput(ns("magma_plot_point_size"), "Point size:",
                           min = 2, max = 10, value = 4, step = 1)
            ),
            .dl_and_download_column(ns, "magma", default_h = 12)
          )
        )
      ),
      br(),
      tags$div(style = "max-width: 1100px; overflow-x: auto;",
        uiOutput(ns("atc_magma_plot_ui"))
      ),
      br(),
      tags$div(style = "max-width: 1100px;",
        h4("Underlying data"),
        DT::DTOutput(ns("atc_magma_tbl"))
      )
    ),
    tabPanel("TWAS-GSEA", br(),
      p("Cross-GWAS TWAS-GSEA drug-class enrichment. Cell colour is signed by ",
        "direction of effect (blue = matches disease signature; red = opposes ",
        "disease signature). Intensity is -log10(FDR). Hatched cells = the ",
        "class was not tested in that GWAS (not the same as 'not significant')."),
      hr(),
      tags$details(class = "gd-details",
        tags$summary("Filter data"),
        tags$div(class = "gd-details-body",
          fluidRow(
            column(3,
              radioButtons(ns("gsea_sig_basis"), "Significance basis:",
                            choices = c("FDR" = "fdr", "P" = "p"),
                            selected = "fdr", inline = TRUE),
              numericInput(ns("gsea_sig_threshold"), "Significance threshold:",
                            value = 0.05, min = 1e-12, max = 1, step = 0.01)
            ),
            column(3,
              checkboxInput(ns("gsea_only_recurrent"),
                            "Only show classes significant in ≥ k GWAS",
                            value = TRUE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("gsea_only_recurrent")),
                sliderInput(ns("gsea_k_min"), "k:",
                             min = 1, max = 9, value = 2, step = 1)
              )
            ),
            column(3,
              selectInput(ns("gsea_panel"), "Panel:",
                          choices = c("Best-per-cell (min P)" = "__best__"),
                          selected = "__best__"),
              selectInput(ns("gsea_gwas_sort"), "Order GWAS by:",
                           choices = .compare_gwas_sort_choices,
                           selected = "as_selected"),
              sliderInput(ns("gsea_plot_font_size"), "Font size (pt):",
                           min = 8, max = 20, value = 11, step = 1),
              sliderInput(ns("gsea_plot_point_size"), "Point size:",
                           min = 2, max = 10, value = 4, step = 1)
            ),
            .dl_and_download_column(ns, "gsea", default_h = 12)
          )
        )
      ),
      br(),
      tags$div(style = "max-width: 1100px; overflow-x: auto;",
        uiOutput(ns("atc_gsea_plot_ui"))
      ),
      gd_legend(list(
        "Blue" = "Direction 'Matches disease' — drug class shares the trait's TWAS signature.",
        "Red"  = "Direction 'Opposes disease' — drug class counteracts the trait's TWAS signature (therapeutic hypothesis).",
        "Grey" = "Nominal or FDR sig with unassigned direction.",
        "Hatched" = "Not tested in this GWAS."
      ), heading = "Colour key"),
      br(),
      tags$div(style = "max-width: 1100px;",
        h4("Underlying data"),
        DT::DTOutput(ns("atc_gsea_tbl"))
      )
    )
  )
}

.atc_row_order <- function(slice, basis = "fdr", thr = 0.05) {
  rec <- slice[, .(
    n_sig = sum(!is.na(.SD[[1L]]) & .SD[[1L]] < thr),
    min_val = suppressWarnings(min(.SD[[1L]], na.rm = TRUE))
  ), by = entity_id, .SDcols = basis]
  rec[order(-n_sig, min_val)]
}

.atc_magma_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                 sig_basis = "fdr", sig_threshold = 0.05,
                                 font_size = 11, point_size = 4) {
  slice <- long[method == "MAGMA-ATC" & gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)
  rec <- .atc_row_order(slice, sig_basis, sig_threshold)
  if (isTRUE(only_recurrent)) rec <- rec[n_sig >= k_min]
  if (nrow(rec) == 0) return(NULL)
  slice <- slice[entity_id %in% rec$entity_id]
  # Attach a description-suffixed label from the first non-NA entity_label
  label_map <- slice[!duplicated(entity_id), .(entity_id, entity_label)]
  slice <- merge(slice, label_map, by = "entity_id", suffixes = c("", ".y"))
  slice[, entity_label := entity_label.y]
  slice[, entity_label.y := NULL]

  slice[, entity_label := factor(entity_label,
    levels = rev(label_map$entity_label[match(rec$entity_id, label_map$entity_id)]))]
  slice[, gwas := factor(gwas, levels = gwas_vec)]
  slice[, minus_log10_fdr := pmin(-log10(fdr), 12)]
  slice[, nom_sig := !is.na(p)   & p   < 0.05]
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = minus_log10_fdr), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(slice, ggplot2::aes(x = gwas, y = entity_label)) +
    base_layer +
    ggplot2::scale_fill_gradient(low = "#e6f4f1", high = .gd_teal,
                                   name = "-log10(FDR)",
                                   na.value = .gd_grey, limits = c(0, 12))

  gg <- .compare_add_sig_overlays(gg, slice, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = c(0.5, 0.5)) +
    ggplot2::scale_y_discrete(expand = c(0.5, 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL,
                   caption = "Ring = nominal-sig (P<0.05). Black square = FDR-sig.") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

.atc_gsea_frame <- function(long, gwas_vec, only_recurrent, k_min,
                              panel_pick = "__best__",
                              sig_basis = "fdr", sig_threshold = 0.05) {
  slice <- long[method == "TWAS-GSEA-ATC" & gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)

  # panel reduction
  if (identical(panel_pick, "__best__") || is.null(panel_pick) ||
      !nzchar(panel_pick)) {
    slice <- pick_best_per_cell(slice, c("gwas", "entity_id"))
  } else {
    slice <- slice[panel == panel_pick]
  }
  if (nrow(slice) == 0) return(NULL)

  # order rows
  rec <- .atc_row_order(slice, sig_basis, sig_threshold)
  if (isTRUE(only_recurrent)) rec <- rec[n_sig >= k_min]
  if (nrow(rec) == 0) return(NULL)
  slice <- slice[entity_id %in% rec$entity_id]

  # Union of entities across selected GWAS -> add "not tested" rows so the
  # heatmap shows the missing cells rather than dropping the row.
  present <- unique(slice[, .(entity_id, gwas)])
  full <- data.table::CJ(entity_id = rec$entity_id, gwas = gwas_vec)
  full <- merge(full, present, by = c("entity_id", "gwas"), all.x = TRUE,
                sort = FALSE)
  full[, tested := paste(entity_id, gwas) %in% paste(present$entity_id, present$gwas)]
  filled <- merge(full, slice, by = c("entity_id", "gwas"), all.x = TRUE,
                  sort = FALSE)

  # labels
  label_map <- slice[!duplicated(entity_id), .(entity_id, entity_label)]
  filled <- merge(filled, label_map, by = "entity_id", all.x = TRUE,
                  suffixes = c("", ".lab"), sort = FALSE)
  filled[!is.na(entity_label.lab), entity_label := entity_label.lab]
  filled[, entity_label.lab := NULL]

  filled[, entity_label := factor(entity_label,
    levels = rev(label_map$entity_label[match(rec$entity_id, label_map$entity_id)]))]
  filled[, gwas := factor(gwas, levels = gwas_vec)]
  filled[, fdr_sig := !is.na(fdr) & fdr < 0.05]
  # Signed -log10(FDR): positive when the class matches the disease
  # signature, negative when it opposes it. Untested cells and cells with
  # NA direction get NA and are handled by na.value / a hatch overlay.
  neg_log_fdr <- pmin(-log10(filled$fdr), 12)
  filled[, signed_score := ifelse(
    is.na(direction) | !tested, NA_real_,
    ifelse(direction == "Matches disease",  neg_log_fdr,
    ifelse(direction == "Opposes disease", -neg_log_fdr, NA_real_)))]
  filled
}

.atc_gsea_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                panel_pick, sig_basis, sig_threshold,
                                font_size = 11, point_size = 4) {
  filled <- .atc_gsea_frame(long, gwas_vec, only_recurrent, k_min,
                             panel_pick, sig_basis, sig_threshold)
  if (is.null(filled)) return(NULL)
  filled[, nom_sig := !is.na(p)   & p   < 0.05]
  filled[, fdr_sig := !is.na(fdr) & fdr < 0.05]
  # Points only drawn for tested cells; blank cells = "not tested".
  tested_dat <- filled[filled$tested, ]

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = signed_score), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(tested_dat, ggplot2::aes(x = gwas, y = entity_label)) +
    base_layer +
    ggplot2::scale_fill_gradient2(
      low = .gd_red, mid = .gd_grey, high = .gd_blue, midpoint = 0,
      limits = c(-12, 12), na.value = .gd_grey,
      name = "signed -log10(FDR)  (+ matches / − opposes)"
    )

  gg <- .compare_add_sig_overlays(gg, tested_dat, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = c(0.5, 0.5)) +
    ggplot2::scale_y_discrete(expand = c(0.5, 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL,
                   caption = "Ring = nominal-sig. Black square = FDR-sig. Blank cell = not tested in that GWAS.") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

atc_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                 comparison_long) {
  moduleServer(id, function(input, output, session) {

    # Update the TWAS-GSEA panel dropdown with the panels actually present.
    observeEvent(comparison_long(), {
      long <- comparison_long()
      panels <- unique(long[method == "TWAS-GSEA-ATC", panel])
      panels <- panels[!is.na(panels)]
      choices <- c("Best-per-cell (min P)" = "__best__",
                    setNames(panels, panels))
      updateSelectInput(session, "gsea_panel",
                         choices = choices, selected = "__best__")
    })

    # Cap the k sliders (one per sub-tab) at the number of selected GWAS.
    observeEvent(selected_gwas_multi(), {
      n_sel <- length(selected_gwas_multi())
      if (n_sel < 1) return()
      for (which_slider in c("magma_k_min", "gsea_k_min")) {
        cur <- input[[which_slider]]
        cur <- if (is.null(cur)) 2L else as.integer(cur)
        updateSliderInput(session, which_slider,
                          max = max(n_sel, 1L),
                          value = min(cur, n_sel))
      }
    })

    magma_gwas_vec <- reactive({
      order_gwas(selected_gwas_multi(),
                  if (is.null(input$magma_gwas_sort)) "as_selected" else input$magma_gwas_sort,
                  NULL)
    })
    gsea_gwas_vec <- reactive({
      order_gwas(selected_gwas_multi(),
                  if (is.null(input$gsea_gwas_sort)) "as_selected" else input$gsea_gwas_sort,
                  NULL)
    })

    # ---------- MAGMA plot / table ----------

    magma_plot <- reactive({
      req(comparison_long())
      .atc_magma_ggplot(
        long           = comparison_long(),
        gwas_vec       = magma_gwas_vec(),
        only_recurrent = isTRUE(input$magma_only_recurrent),
        k_min          = if (is.null(input$magma_k_min)) 2L else as.integer(input$magma_k_min),
        sig_basis      = if (is.null(input$magma_sig_basis)) "fdr" else input$magma_sig_basis,
        sig_threshold  = if (is.null(input$magma_sig_threshold)) 0.05 else as.numeric(input$magma_sig_threshold),
        font_size      = if (is.null(input$magma_plot_font_size)) 11 else as.numeric(input$magma_plot_font_size),
        point_size     = if (is.null(input$magma_plot_point_size)) 4 else as.numeric(input$magma_plot_point_size)
      )
    })

    magma_dims <- reactive({
      p <- magma_plot()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_label))
      .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = length(magma_gwas_vec()),
        font_size = if (is.null(input$magma_plot_font_size)) 11 else as.numeric(input$magma_plot_font_size),
        y_labels  = if (is.null(p)) NULL else as.character(unique(p$data$entity_label))
      )
    })

    output$atc_magma_plot_ui <- renderUI({
      dims <- magma_dims()
      plotOutput(session$ns("atc_magma_plot"),
                  height = dims$height, width = dims$width)
    })

    output$atc_magma_plot <- renderPlot({
      p <- magma_plot()
      if (is.null(p)) {
        plot.new(); title("No ATC classes meet the current filter")
        return(invisible())
      }
      print(p)
    })

    output$atc_magma_tbl <- DT::renderDT({
      req(comparison_long())
      gwas_vec <- magma_gwas_vec()
      slice <- comparison_long()[method == "MAGMA-ATC" & gwas %in% gwas_vec]
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS   = factor(slice$gwas, levels = gwas_vec),
        `ATC Code`  = slice$entity_id,
        Class       = slice$entity_label,
        `N Drugs`   = slice$n_units,
        P           = signif(slice$p, 3),
        `P.FDR`     = signif(slice$fdr, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(5, "asc"))))
    }, server = TRUE)

    output$magma_download_plot <- downloadHandler(
      filename = function() sprintf("atc_magma_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$magma_dl_format),
      content = function(file) {
        p <- magma_plot(); if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$magma_dl_format
        w <- input$magma_dl_width; h <- input$magma_dl_height
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$magma_download_csv <- downloadHandler(
      filename = function() sprintf("atc_magma_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- magma_gwas_vec()
        long <- comparison_long()[method == "MAGMA-ATC" & gwas %in% gwas_vec]
        wide <- pivot_matrix(long, "fdr", gwas_vec)
        header <- sprintf("# GenoDisc ATC MAGMA compare CSV | %s | sig_basis=%s threshold=%g k_min=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           input$magma_sig_basis %||% "fdr",
                           input$magma_sig_threshold %||% 0.05,
                           input$magma_k_min %||% 2L)
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: P.FDR", con)
        out <- data.frame(`ATC Code` = rownames(wide),
                          format(wide, scientific = TRUE, digits = 3),
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )

    # ---------- TWAS-GSEA plot / table ----------

    gsea_plot <- reactive({
      req(comparison_long())
      .atc_gsea_ggplot(
        long           = comparison_long(),
        gwas_vec       = gsea_gwas_vec(),
        only_recurrent = isTRUE(input$gsea_only_recurrent),
        k_min          = if (is.null(input$gsea_k_min)) 2L else as.integer(input$gsea_k_min),
        panel_pick     = input$gsea_panel,
        sig_basis      = if (is.null(input$gsea_sig_basis)) "fdr" else input$gsea_sig_basis,
        sig_threshold  = if (is.null(input$gsea_sig_threshold)) 0.05 else as.numeric(input$gsea_sig_threshold),
        font_size      = if (is.null(input$gsea_plot_font_size)) 11 else as.numeric(input$gsea_plot_font_size),
        point_size     = if (is.null(input$gsea_plot_point_size)) 4 else as.numeric(input$gsea_plot_point_size)
      )
    })

    gsea_dims <- reactive({
      p <- gsea_plot()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_label))
      .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = length(gsea_gwas_vec()),
        font_size = if (is.null(input$gsea_plot_font_size)) 11 else as.numeric(input$gsea_plot_font_size),
        y_labels  = if (is.null(p)) NULL else as.character(unique(p$data$entity_label))
      )
    })

    output$atc_gsea_plot_ui <- renderUI({
      dims <- gsea_dims()
      plotOutput(session$ns("atc_gsea_plot"),
                  height = dims$height, width = dims$width)
    })

    output$atc_gsea_plot <- renderPlot({
      p <- gsea_plot()
      if (is.null(p)) {
        plot.new(); title("No ATC classes meet the current filter")
        return(invisible())
      }
      print(p)
    })

    output$atc_gsea_tbl <- DT::renderDT({
      req(comparison_long())
      gwas_vec <- gsea_gwas_vec()
      slice <- comparison_long()[method == "TWAS-GSEA-ATC" & gwas %in% gwas_vec]
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS       = factor(slice$gwas, levels = gwas_vec),
        `ATC Code` = slice$entity_id,
        Class      = slice$entity_label,
        Panel      = slice$panel,
        Estimate   = signif(slice$statistic, 3),
        Direction  = slice$direction,
        Reversal_Z = signif(slice$reversal_z, 3),
        P          = signif(slice$p, 3),
        `P.FDR`    = signif(slice$fdr, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(8, "asc"))))
    }, server = TRUE)

    output$gsea_download_plot <- downloadHandler(
      filename = function() sprintf("atc_gsea_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$gsea_dl_format),
      content = function(file) {
        p <- gsea_plot(); if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$gsea_dl_format
        w <- input$gsea_dl_width; h <- input$gsea_dl_height
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$gsea_download_csv <- downloadHandler(
      filename = function() sprintf("atc_gsea_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- gsea_gwas_vec()
        long <- comparison_long()[method == "TWAS-GSEA-ATC" & gwas %in% gwas_vec]
        best <- pick_best_per_cell(long, c("gwas", "entity_id"))
        best[, signed := ifelse(!is.na(direction) & direction == "Opposes disease",
                                 -(-log10(fdr)), -log10(fdr))]
        wide <- pivot_matrix(best, "signed", gwas_vec)
        display <- ifelse(is.na(wide), "NT", sprintf("%+.2f", wide))
        display <- matrix(display, nrow = nrow(wide), dimnames = dimnames(wide))
        header <- sprintf("# GenoDisc ATC TWAS-GSEA compare CSV | %s | sig_basis=%s threshold=%g k_min=%s panel=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           input$gsea_sig_basis %||% "fdr",
                           input$gsea_sig_threshold %||% 0.05,
                           input$gsea_k_min %||% 2L,
                           input$gsea_panel %||% "__best__")
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: signed -log10(FDR): + = matches, - = opposes. NT = not tested.", con)
        out <- data.frame(`ATC Code` = rownames(display), display,
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )
  })
}

