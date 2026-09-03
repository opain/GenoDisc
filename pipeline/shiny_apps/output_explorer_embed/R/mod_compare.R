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
            sliderInput(ns("k_min"), "k:",
                         min = 1, max = 9, value = 2, step = 1)
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
                         min = 8, max = 20, value = 12, step = 1)
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
      plotOutput(ns("tissue_compare_plot"), height = "720px")
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

  if (isTRUE(only_recurrent) && k_min > 1L) {
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
                                    font_size = 12) {
  packed <- .tissue_compare_frame(long, gwas_vec, only_recurrent, k_min,
                                    sig_basis, sig_threshold)
  if (is.null(packed)) return(NULL)
  slice <- packed$slice

  fill_map <- c(
    "Not significant" = .gd_grey,
    "Nominal-sig"     = .gd_hex_alpha(.gd_teal, 0.35),
    "FDR-sig"         = .gd_hex_alpha(.gd_teal, 0.7),
    "Retained"        = .gd_teal
  )

  intensity <- if (identical(metric, "p")) -log10(slice$p) else -log10(slice$fdr)
  slice[, label := ifelse(!is.na(intensity) & intensity >= -log10(0.05),
                          sprintf("%.1f", pmin(intensity, 20)), "")]

  ggplot2::ggplot(slice, ggplot2::aes(x = gwas, y = entity_id, fill = tier)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3, colour = "black") +
    ggplot2::scale_fill_manual(values = fill_map, drop = FALSE,
                                name = "Cell tier") +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = c(0, 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL,
                   caption = "Numbers show -log10(FDR) where cell is FDR- or retained-significant.") +
    ggplot2::theme_minimal(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(t = 70, r = 40, b = 10, l = 10, unit = "pt")
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
        font_size      = if (is.null(input$plot_font_size)) 12 else as.numeric(input$plot_font_size)
      )
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
              sliderInput(ns("magma_k_min"), "k:",
                           min = 1, max = 9, value = 2, step = 1)
            ),
            column(3,
              selectInput(ns("magma_gwas_sort"), "Order GWAS by:",
                           choices = .compare_gwas_sort_choices,
                           selected = "as_selected"),
              sliderInput(ns("magma_plot_font_size"), "Font size (pt):",
                           min = 8, max = 20, value = 11, step = 1)
            ),
            .dl_and_download_column(ns, "magma", default_h = 12)
          )
        )
      ),
      br(),
      tags$div(style = "max-width: 1100px; overflow-x: auto;",
        plotOutput(ns("atc_magma_plot"), height = "800px")
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
              sliderInput(ns("gsea_k_min"), "k:",
                           min = 1, max = 9, value = 2, step = 1)
            ),
            column(3,
              selectInput(ns("gsea_panel"), "Panel:",
                          choices = c("Best-per-cell (min P)" = "__best__"),
                          selected = "__best__"),
              selectInput(ns("gsea_gwas_sort"), "Order GWAS by:",
                           choices = .compare_gwas_sort_choices,
                           selected = "as_selected"),
              sliderInput(ns("gsea_plot_font_size"), "Font size (pt):",
                           min = 8, max = 20, value = 11, step = 1)
            ),
            .dl_and_download_column(ns, "gsea", default_h = 12)
          )
        )
      ),
      br(),
      tags$div(style = "max-width: 1100px; overflow-x: auto;",
        plotOutput(ns("atc_gsea_plot"), height = "800px")
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

# Colour-map helper: signed -log10(FDR) into blue (matches) / red (opposes)
.atc_gsea_fill <- function(direction, fdr) {
  intensity <- pmin(-log10(fdr), 12) / 12
  intensity[is.na(intensity)] <- 0
  # start from grey and shift toward direction colour
  out <- rep(.gd_grey, length(direction))
  m <- !is.na(direction) & direction == "Matches disease"
  o <- !is.na(direction) & direction == "Opposes disease"
  ramp_b <- grDevices::colorRamp(c(.gd_grey, .gd_blue))
  ramp_r <- grDevices::colorRamp(c(.gd_grey, .gd_red))
  if (any(m)) {
    rgb_m <- ramp_b(intensity[m])
    out[m] <- grDevices::rgb(rgb_m[,1]/255, rgb_m[,2]/255, rgb_m[,3]/255)
  }
  if (any(o)) {
    rgb_o <- ramp_r(intensity[o])
    out[o] <- grDevices::rgb(rgb_o[,1]/255, rgb_o[,2]/255, rgb_o[,3]/255)
  }
  out
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
                                 font_size = 11) {
  slice <- long[method == "MAGMA-ATC" & gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)
  rec <- .atc_row_order(slice, sig_basis, sig_threshold)
  if (isTRUE(only_recurrent) && k_min > 1L) rec <- rec[n_sig >= k_min]
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
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]

  ggplot2::ggplot(slice, ggplot2::aes(x = gwas, y = entity_label,
                                        fill = minus_log10_fdr)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_tile(data = slice[fdr_sig == TRUE],
                        fill = NA, colour = "black", linewidth = 0.7) +
    ggplot2::scale_fill_gradient(low = "#e6f4f1", high = .gd_teal,
                                   name = "-log10(FDR)",
                                   na.value = .gd_grey, limits = c(0, 12)) +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = c(0, 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(t = 70, r = 40, b = 10, l = 10, unit = "pt")
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
  if (isTRUE(only_recurrent) && k_min > 1L) rec <- rec[n_sig >= k_min]
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
  filled[, fill_hex := .atc_gsea_fill(direction, fdr)]
  filled$fill_hex[!filled$tested] <- NA
  filled
}

.atc_gsea_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                panel_pick, sig_basis, sig_threshold,
                                font_size = 11) {
  filled <- .atc_gsea_frame(long, gwas_vec, only_recurrent, k_min,
                             panel_pick, sig_basis, sig_threshold)
  if (is.null(filled)) return(NULL)
  p <- ggplot2::ggplot(filled, ggplot2::aes(x = gwas, y = entity_label)) +
    ggplot2::geom_tile(ggplot2::aes(fill = I(fill_hex)),
                        colour = "white", linewidth = 0.3) +
    ggplot2::geom_tile(data = filled[tested == FALSE],
                        fill = .gd_hatch_bg, colour = .gd_hatch_fg,
                        linetype = "dotted", linewidth = 0.6) +
    ggplot2::geom_tile(data = filled[fdr_sig == TRUE & tested == TRUE],
                        fill = NA, colour = "black", linewidth = 0.7) +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = c(0, 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 70, r = 40, b = 10, l = 10, unit = "pt")
    )
  p
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
        font_size      = if (is.null(input$magma_plot_font_size)) 11 else as.numeric(input$magma_plot_font_size)
      )
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
        font_size      = if (is.null(input$gsea_plot_font_size)) 11 else as.numeric(input$gsea_plot_font_size)
      )
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

