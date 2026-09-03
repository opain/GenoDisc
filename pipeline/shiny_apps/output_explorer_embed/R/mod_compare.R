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
    fluidRow(
      column(3,
        checkboxInput(ns("only_recurrent_tissue"),
                      "Only show tissues significant in ≥ k GWAS", value = FALSE),
        helpText("k is set at the top of the page.")
      ),
      column(3,
        radioButtons(ns("cell_metric_tissue"), "Cell colour:",
                     choices = c("-log10(FDR)" = "fdr", "-log10(P)" = "p"),
                     selected = "fdr", inline = TRUE)
      ),
      column(3,
        numericInput(ns("dl_width_tissue_cmp"), "Download width (in):",
                     value = 10, min = 4, max = 24, step = 0.5),
        numericInput(ns("dl_height_tissue_cmp"), "Download height (in):",
                     value = 10, min = 4, max = 24, step = 0.5)
      ),
      column(3,
        selectInput(ns("dl_format_tissue_cmp"), "Download format:",
                    choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                    selected = "png"),
        downloadButton(ns("download_plot_tissue_cmp"), "Download plot"),
        br(), br(),
        downloadButton(ns("download_csv_tissue_cmp"), "Download matrix CSV")
      )
    ),
    br(),
    tags$div(style = "overflow-x: auto;",
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
    h4("Underlying data"),
    DT::DTOutput(ns("tissue_compare_tbl"))
  )
}

# Long -> tidy plotting frame with tier categorisation
.tissue_compare_frame <- function(long, gwas_vec, only_recurrent, k_min) {
  # Slice to tissue rows and to the selected GWAS
  slice <- long[method == "MAGMA-tissue" & gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)

  # Determine row order: recurrence (FDR<0.05) desc, then min p asc.
  rec <- slice[, .(
    n_fdr_sig = sum(!is.na(fdr) & fdr < 0.05),
    n_nom_sig = sum(!is.na(p) & p < 0.05),
    min_p     = suppressWarnings(min(p, na.rm = TRUE))
  ), by = entity_id]
  rec <- rec[order(-n_fdr_sig, -n_nom_sig, min_p)]

  if (isTRUE(only_recurrent) && k_min > 1L) {
    rec <- rec[n_fdr_sig >= k_min | n_nom_sig >= k_min]
  }

  if (nrow(rec) == 0) return(NULL)
  slice <- slice[entity_id %in% rec$entity_id]
  slice[, entity_id := factor(entity_id, levels = rev(rec$entity_id))]
  slice[, gwas := factor(gwas, levels = gwas_vec)]

  # Three-tier: retained > FDR-sig > nominal > not sig
  slice[, tier := "Not significant"]
  slice[!is.na(p) & p < 0.05,        tier := "Nominal-sig"]
  slice[!is.na(fdr) & fdr < 0.05,    tier := "FDR-sig"]
  slice[!is.na(evidence) & evidence, tier := "Retained"]
  slice[, tier := factor(tier, levels = c("Not significant", "Nominal-sig",
                                            "FDR-sig", "Retained"))]
  list(slice = slice, rec = rec)
}

.tissue_compare_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                    metric = "fdr") {
  packed <- .tissue_compare_frame(long, gwas_vec, only_recurrent, k_min)
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
    ggplot2::labs(x = NULL, y = NULL,
                   caption = "Numbers show -log10(FDR) where cell is FDR- or retained-significant.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(35, 25, 5, 5, "pt")
    )
}

tissue_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                    shared_filters, comparison_long) {
  moduleServer(id, function(input, output, session) {

    plot_obj <- reactive({
      req(comparison_long())
      sf <- shared_filters()
      gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
      .tissue_compare_ggplot(
        long           = comparison_long(),
        gwas_vec       = gwas_vec,
        only_recurrent = isTRUE(input$only_recurrent_tissue),
        k_min          = as.integer(sf$k_min %||% 2L),
        metric         = if (is.null(input$cell_metric_tissue)) "fdr" else input$cell_metric_tissue
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
      sf <- shared_filters()
      gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
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

    output$download_plot_tissue_cmp <- downloadHandler(
      filename = function() sprintf("tissue_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$dl_format_tissue_cmp),
      content = function(file) {
        p <- plot_obj()
        if (is.null(p)) {
          # write a placeholder so the download completes
          grDevices::png(file, width = 400, height = 200); dev.off(); return()
        }
        fmt <- input$dl_format_tissue_cmp
        w <- input$dl_width_tissue_cmp; h <- input$dl_height_tissue_cmp
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

    output$download_csv_tissue_cmp <- downloadHandler(
      filename = function() sprintf("tissue_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        sf <- shared_filters()
        gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
        long <- comparison_long()[method == "MAGMA-tissue" & gwas %in% gwas_vec]
        # matrix P.FDR + retained flag concatenated as "0.012 (R)"
        wide <- pivot_matrix(long, "fdr", gwas_vec)
        ret  <- pivot_matrix(long, "evidence", gwas_vec)   # 1/0 encoding
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
                          sf$sig_basis, sf$sig_threshold, sf$k_min %||% "NA")
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

atc_compare_ui <- function(ns) {
  tabsetPanel(
    tabPanel("MAGMA", br(),
      p("Cross-GWAS MAGMA drug-class enrichment. Cell colour intensity is ",
        "-log10(FDR); cells with FDR < 0.05 are outlined in black."),
      hr(),
      fluidRow(
        column(3,
          checkboxInput(ns("only_recurrent_atc_magma"),
                        "Only show classes significant in ≥ k GWAS",
                        value = TRUE)
        ),
        column(3,
          numericInput(ns("dl_width_atc_magma"), "Download width (in):",
                       value = 10, min = 4, max = 24, step = 0.5),
          numericInput(ns("dl_height_atc_magma"), "Download height (in):",
                       value = 12, min = 4, max = 32, step = 0.5)
        ),
        column(3,
          selectInput(ns("dl_format_atc_magma"), "Download format:",
                      choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                      selected = "png"),
          downloadButton(ns("download_plot_atc_magma"), "Download plot"),
          br(), br(),
          downloadButton(ns("download_csv_atc_magma"), "Download matrix CSV")
        )
      ),
      br(),
      tags$div(style = "overflow-x: auto;",
        plotOutput(ns("atc_magma_plot"), height = "800px")
      ),
      br(),
      h4("Underlying data"),
      DT::DTOutput(ns("atc_magma_tbl"))
    ),
    tabPanel("TWAS-GSEA", br(),
      p("Cross-GWAS TWAS-GSEA drug-class enrichment. Cell colour is signed by ",
        "direction of effect (blue = matches disease signature; red = opposes ",
        "disease signature). Intensity is -log10(FDR). Hatched cells = the ",
        "class was not tested in that GWAS (not the same as 'not significant')."),
      hr(),
      fluidRow(
        column(3,
          selectInput(ns("atc_gsea_panel"), "Panel:",
                      choices = c("Best-per-cell (min P)" = "__best__"),
                      selected = "__best__")
        ),
        column(3,
          checkboxInput(ns("only_recurrent_atc_gsea"),
                        "Only show classes significant in ≥ k GWAS",
                        value = TRUE)
        ),
        column(3,
          numericInput(ns("dl_width_atc_gsea"), "Download width (in):",
                       value = 10, min = 4, max = 24, step = 0.5),
          numericInput(ns("dl_height_atc_gsea"), "Download height (in):",
                       value = 12, min = 4, max = 32, step = 0.5)
        ),
        column(3,
          selectInput(ns("dl_format_atc_gsea"), "Download format:",
                      choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                      selected = "png"),
          downloadButton(ns("download_plot_atc_gsea"), "Download plot"),
          br(), br(),
          downloadButton(ns("download_csv_atc_gsea"), "Download matrix CSV")
        )
      ),
      br(),
      tags$div(style = "overflow-x: auto;",
        plotOutput(ns("atc_gsea_plot"), height = "800px")
      ),
      gd_legend(list(
        "Blue" = "Direction 'Matches disease' — drug class shares the trait's TWAS signature.",
        "Red"  = "Direction 'Opposes disease' — drug class counteracts the trait's TWAS signature (therapeutic hypothesis).",
        "Grey" = "Nominal or FDR sig with unassigned direction.",
        "Hatched" = "Not tested in this GWAS."
      ), heading = "Colour key"),
      br(),
      h4("Underlying data"),
      DT::DTOutput(ns("atc_gsea_tbl"))
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
                                 sig_basis = "fdr", sig_threshold = 0.05) {
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
    ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
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
                                panel_pick, sig_basis, sig_threshold) {
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
    ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )
  p
}

atc_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                 shared_filters, comparison_long) {
  moduleServer(id, function(input, output, session) {

    # Update panel dropdown with the panels actually present in the data.
    observeEvent(comparison_long(), {
      long <- comparison_long()
      panels <- unique(long[method == "TWAS-GSEA-ATC", panel])
      panels <- panels[!is.na(panels)]
      choices <- c("Best-per-cell (min P)" = "__best__",
                    setNames(panels, panels))
      updateSelectInput(session, "atc_gsea_panel",
                         choices = choices, selected = "__best__")
    })

    # ---------- MAGMA plot / table ----------

    magma_plot <- reactive({
      req(comparison_long())
      sf <- shared_filters()
      gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
      .atc_magma_ggplot(
        long           = comparison_long(),
        gwas_vec       = gwas_vec,
        only_recurrent = isTRUE(input$only_recurrent_atc_magma),
        k_min          = as.integer(sf$k_min %||% 2L),
        sig_basis      = sf$sig_basis,
        sig_threshold  = as.numeric(sf$sig_threshold %||% 0.05)
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
      sf <- shared_filters()
      gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
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

    output$download_plot_atc_magma <- downloadHandler(
      filename = function() sprintf("atc_magma_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$dl_format_atc_magma),
      content = function(file) {
        p <- magma_plot(); if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$dl_format_atc_magma
        w <- input$dl_width_atc_magma; h <- input$dl_height_atc_magma
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$download_csv_atc_magma <- downloadHandler(
      filename = function() sprintf("atc_magma_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        sf <- shared_filters()
        gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
        long <- comparison_long()[method == "MAGMA-ATC" & gwas %in% gwas_vec]
        wide <- pivot_matrix(long, "fdr", gwas_vec)
        header <- sprintf("# GenoDisc ATC MAGMA compare CSV | %s | sig_basis=%s threshold=%g k_min=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           sf$sig_basis, sf$sig_threshold, sf$k_min %||% "NA")
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
      sf <- shared_filters()
      gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
      .atc_gsea_ggplot(
        long           = comparison_long(),
        gwas_vec       = gwas_vec,
        only_recurrent = isTRUE(input$only_recurrent_atc_gsea),
        k_min          = as.integer(sf$k_min %||% 2L),
        panel_pick     = input$atc_gsea_panel,
        sig_basis      = sf$sig_basis,
        sig_threshold  = as.numeric(sf$sig_threshold %||% 0.05)
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
      sf <- shared_filters()
      gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
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

    output$download_plot_atc_gsea <- downloadHandler(
      filename = function() sprintf("atc_gsea_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$dl_format_atc_gsea),
      content = function(file) {
        p <- gsea_plot(); if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$dl_format_atc_gsea
        w <- input$dl_width_atc_gsea; h <- input$dl_height_atc_gsea
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$download_csv_atc_gsea <- downloadHandler(
      filename = function() sprintf("atc_gsea_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        sf <- shared_filters()
        gwas_vec <- order_gwas(selected_gwas_multi(), sf$gwas_sort, NULL)
        long <- comparison_long()[method == "TWAS-GSEA-ATC" & gwas %in% gwas_vec]
        best <- pick_best_per_cell(long, c("gwas", "entity_id"))
        # Signed -log10(FDR) or NT tag
        best[, signed := ifelse(!is.na(direction) & direction == "Opposes disease",
                                 -(-log10(fdr)), -log10(fdr))]
        wide <- pivot_matrix(best, "signed", gwas_vec)
        display <- ifelse(is.na(wide), "NT", sprintf("%+.2f", wide))
        display <- matrix(display, nrow = nrow(wide), dimnames = dimnames(wide))
        header <- sprintf("# GenoDisc ATC TWAS-GSEA compare CSV | %s | sig_basis=%s threshold=%g k_min=%s panel=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           sf$sig_basis, sf$sig_threshold, sf$k_min %||% "NA",
                           input$atc_gsea_panel %||% "__best__")
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

