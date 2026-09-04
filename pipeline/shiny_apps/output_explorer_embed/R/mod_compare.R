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

# NA-aware numeric formatter for DT tables. `signif(NA, 3)` returns
# NA_real_ which DT renders as an empty cell; users can't distinguish
# "value missing" from "cell not populated". Convert to a character
# vector so NA becomes the literal string "NA".
.fmt_sig <- function(v, digits = 3) {
  ifelse(is.na(v), "NA",
         formatC(signif(as.numeric(v), digits), format = "g", digits = digits))
}
.fmt_int <- function(v) {
  ifelse(is.na(v), "NA", formatC(as.integer(v), format = "d", big.mark = ","))
}

# Plot dimensions for a compare heatmap given its row and column counts. Mimics
# the sizing used by calc_plot_dims() for the single-GWAS heatmaps but does not
# require a facet column. Returns pixel dimensions suitable for plotOutput().
.compare_plot_dims <- function(n_rows, n_cols, font_size = 11,
                                 y_labels = NULL,
                                 min_height = 300, min_width = 500,
                                 min_panel_h_pt = 120) {
  fs_ratio <- font_size / 11
  y_label_px <- if (!is.null(y_labels) && length(y_labels) > 0) {
    max(nchar(as.character(y_labels)), na.rm = TRUE) * font_size * 0.55 + 20
  } else 160
  panel_w <- max(30, 40 * fs_ratio) * max(1L, n_cols)
  # Panel height scales with row count but is floored at min_panel_h_pt so
  # short heatmaps (1-3 rows) still leave room for the bottom legend rather
  # than being crushed.
  panel_h <- max(min_panel_h_pt, max(18, 22 * fs_ratio) * max(1L, n_rows))
  x_label_h <- 6 * font_size + 30
  legend_h  <- 60   # bottom colour-bar legend
  height <- max(min_height, x_label_h + panel_h + legend_h + 40)
  # Right margin is small now that the legend sits underneath the panel.
  width  <- max(min_width, y_label_px + panel_w + 60)
  list(height = round(height), width = round(width))
}

# Placeholder shown in a compare view's plot slot when the filter combination
# leaves no rows. A helpful text message with concrete suggestions is more
# actionable than an empty ggplot with a small caption. `entity_label` is the
# entity noun for the calling view (e.g. "loci", "tissues", "genes").
.compare_empty_state <- function(entity_label = "entities") {
  tags$div(
    class = "well",
    style = "max-width: 720px; margin-top: 12px; color: var(--gd-text-mute);",
    tags$p(
      tags$b(sprintf("No %s match the current filters.", entity_label)),
      " Try changing the ", tags$b("Significance threshold"),
      ", lowering the ", tags$b("k"),
      " value (or unticking ", tags$b("Only show ... significant in ≥ k GWAS"),
      "), or expanding the ", tags$b("Method"),
      " / ", tags$b("Panel"), " selection."
    )
  )
}

# Modal-dialog helper that shows every comparison_long row for a given
# (entity_type, entity_id) pair across all methods, panels and selected
# GWAS. Called from the row-click observers on each compare view's DT
# table. Answers "is this entity real across methods / traits, or one
# LD block / one method?".
.show_entity_detail <- function(ent_id, ent_type, long, gwas_vec) {
  if (is.null(ent_id) || is.na(ent_id) || !nzchar(ent_id)) return(invisible())
  detail <- long[entity_type == ent_type &
                    entity_id   == ent_id &
                    gwas %in% gwas_vec]
  if (nrow(detail) == 0) {
    showModal(modalDialog(
      title = paste0(ent_id, " — no rows in current selection"),
      easyClose = TRUE, size = "m"
    ))
    return(invisible())
  }
  detail <- detail[order(match(gwas, gwas_vec), method, panel)]
  # Build the display columns using NA-aware formatters so numeric NAs
  # render as "NA" (not blank / error) and character columns use "—" for
  # missing values.
  fmt_char <- function(v) ifelse(is.na(v) | !nzchar(as.character(v)),
                                    "—", as.character(v))
  fmt_bool <- function(v) ifelse(is.na(v), "—", ifelse(v, "Yes", "No"))
  out <- data.frame(
    GWAS       = factor(detail$gwas, levels = gwas_vec),
    Method     = fmt_char(detail$method),
    Panel      = fmt_char(detail$panel),
    Statistic  = .fmt_sig(detail$statistic, 3),
    SE         = .fmt_sig(detail$se, 3),
    P          = .fmt_sig(detail$p, 3),
    `P.FDR`    = .fmt_sig(detail$fdr, 3),
    Direction  = fmt_char(detail$direction),
    Evidence   = fmt_bool(detail$evidence),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  # Drop columns whose values are all placeholders (NA / "—" / "NA")
  # so they don't clutter the modal. isTRUE() on `any(...)` guards
  # against NA propagating into the selector — the old vapply pattern
  # returned NA for columns that were entirely NA, and passing a
  # logical vector containing NA to `[.data.frame` raises "undefined
  # columns selected", which is what was crashing the app on row-click.
  keep <- vapply(names(out), function(nm) {
    v <- as.character(out[[nm]])
    isTRUE(any(!is.na(v) & v != "—" & v != "NA"))
  }, logical(1))
  # GWAS column is always kept regardless.
  keep["GWAS"] <- TRUE
  out <- out[, keep, drop = FALSE]
  first_lab <- detail$entity_label[1L]
  title_txt <- if (nzchar(first_lab) && !identical(first_lab, ent_id)) {
    sprintf("%s — %s", ent_id, first_lab)
  } else {
    ent_id
  }
  showModal(modalDialog(
    title = title_txt,
    size  = "l",
    easyClose = TRUE,
    footer = modalButton("Close"),
    tags$p(
      style = "color: var(--gd-text-mute); margin-bottom: 8px;",
      sprintf("All %s rows for %s across the %d selected GWAS.",
              ent_type, ent_id, length(gwas_vec))),
    DT::datatable(
      out, rownames = FALSE,
      options = list(pageLength = 25, dom = "tip",
                      order = list(list(0, "asc")))
    )
  ))
}

# Sync the download Width/Height numericInputs to the current on-screen plot
# dimensions. Whenever `dims_reactive()` changes (font size, row count, GWAS
# selection etc.), the numericInputs update to match. User edits get
# overwritten on the next dims change - that's acceptable since the correct
# default is more useful than persistent overrides.
.sync_dl_dims <- function(session, dims_reactive, prefix = "") {
  w_key <- if (nzchar(prefix)) paste0(prefix, "_dl_width")  else "dl_width"
  h_key <- if (nzchar(prefix)) paste0(prefix, "_dl_height") else "dl_height"
  observe({
    dims <- dims_reactive()
    w_in <- max(2, round(dims$width  / 96 * 2) / 2)
    h_in <- max(2, round(dims$height / 96 * 2) / 2)
    updateNumericInput(session, w_key, value = w_in)
    updateNumericInput(session, h_key, value = h_in)
  })
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
      "Default plot lays out each GWAS in its own column as a lollipop of ",
      "-log10(P), with tissues shared across facets and ordered by recurrence. ",
      "Switch to the heatmap when comparing many GWAS."
    ),
    hr(),
    tags$details(class = "gd-details",
      tags$summary("Filter data"),
      tags$div(class = "gd-details-body",
        fluidRow(
          column(3,
            radioButtons(ns("plot_type"), "Plot type:",
                          choices = c("Facet by GWAS" = "facet",
                                       "Heatmap" = "heatmap"),
                          selected = "facet", inline = TRUE),
            numericInput(ns("sig_threshold"), "FDR significance threshold:",
                          value = 0.05, min = 1e-12, max = 1, step = 0.01),
            # `gwas_pick` / `tissue_pick` are rendered server-side via
            # renderUI so the choices are baked in at widget-creation
            # time (see gencor commit fa9d73e for why updateSelectInput
            # is unreliable here — tissue_compare_ui also lives inside
            # a renderUI swap in mod_enrichment).
            uiOutput(ns("gwas_pick_ui"))
          ),
          column(3,
            uiOutput(ns("tissue_pick_ui")),
            checkboxInput(ns("only_recurrent"),
                          "Only show tissues significant in ≥ k GWAS",
                          value = FALSE),
            conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("only_recurrent")),
              uiOutput(ns("k_slider_ui"))
            )
          ),
          column(3,
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
      "Facet mode"                   = "One panel per GWAS; each row is a tissue and the horizontal bar shows -log10(P). Dashed line = nominal significance (P = 0.05); dotted line = Bonferroni across shown tissues.",
      "Filled point (facet)"         = "Green fill = FDR-significant (P.FDR < 0.05); white fill = not FDR-significant.",
      "Inner white dot (black outline)" = "Retained after the conditional analysis (shown in both facet and heatmap modes).",
      "Cell colour (heatmap)"        = "-log10(P) on a teal ramp.",
      "Ring around a circle"         = "Heatmap: nominal-significant (P < 0.05).",
      "Black square around a circle" = "Heatmap: FDR-significant (P.FDR < 0.05)."
    ), heading = "How to read this plot"),
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
  basis_col <- "fdr"
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

#' Faceted per-GWAS lollipop for the Tissue compare view.
#'
#' Uses the same visual language as the single-GWAS build_tissue_plot()
#' but places each selected GWAS in its own column facet. Tissue order is
#' shared across facets so recurrent hits stack together.
#'
#' @param long comparison_long tibble (`method == "MAGMA-tissue"` slice).
#' @param gwas_vec GWAS order for facets.
#' @param only_recurrent,k_min,sig_threshold Passed through the same recurrence
#'   filter as the heatmap so both views react identically to the sidebar.
#' @param sort_choice "significance" (default) or "alphabetical" for tissue order.
.tissue_compare_facet_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                          sig_threshold = 0.05,
                                          sort_choice = "significance",
                                          font_size = 12,
                                          point_size = 3) {
  packed <- .tissue_compare_frame(long, gwas_vec, only_recurrent, k_min,
                                    "fdr", sig_threshold)
  if (is.null(packed)) return(NULL)
  slice <- data.table::copy(packed$slice)
  slice[, negLog10P := -log10(pmax(p, 1e-300))]
  slice[, FDR_Sig   := !is.na(fdr) & fdr < 0.05]
  slice[, Retained  := !is.na(evidence) & evidence]

  # Shared tissue ordering across facets. Recurrence-based row order from
  # `.tissue_compare_frame` puts recurrent hits first (via `rec$entity_id`),
  # so re-use that; `sort_choice = "alphabetical"` overrides to trait name.
  tissue_levels <- if (identical(sort_choice, "alphabetical")) {
    sort(as.character(unique(slice$entity_id)), decreasing = TRUE)
  } else {
    # rec$entity_id is high-recurrence-first — we want that at the top of
    # the facet (so reverse for factor levels — first level renders at bottom).
    rev(as.character(packed$rec$entity_id))
  }
  slice[, Label := factor(as.character(entity_id), levels = tissue_levels)]
  slice[, gwas  := factor(gwas, levels = gwas_vec)]

  # Bonferroni line: total tissues in the CURRENT slice (per-facet population).
  n_total  <- length(unique(slice$entity_id))
  nom_line  <- -log10(0.05)
  bonf_line <- -log10(0.05 / max(n_total, 1))

  gg <- ggplot2::ggplot(slice, ggplot2::aes(x = negLog10P, y = Label)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = negLog10P, yend = Label),
                          colour = "grey78", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = nom_line,  linetype = "dashed", colour = "grey55") +
    ggplot2::geom_vline(xintercept = bonf_line, linetype = "dotted", colour = "grey35") +
    ggplot2::geom_point(ggplot2::aes(fill = FDR_Sig),
                        shape = 21, size = point_size, colour = "black", stroke = 0.4) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "white", `TRUE` = "#0f766e"),
      labels = c(`FALSE` = "Not FDR-significant", `TRUE` = "FDR-significant"),
      name = NULL)

  # Retained-in-conditional marker: small inner WHITE dot with a thin
  # black outline (shape 21). White-on-dark-green (FDR-sig fill) and
  # white-on-white (not-sig fill) both read clearly because of the
  # black stroke.
  ret <- slice[Retained == TRUE]
  if (nrow(ret) > 0) {
    gg <- gg + ggplot2::geom_point(data = ret, shape = 21, fill = "white",
                                     colour = "black", stroke = 0.6,
                                     size = point_size * 0.45)
  }

  gg +
    ggplot2::facet_grid(cols = ggplot2::vars(gwas), scales = "free_x") +
    ggplot2::labs(x = expression(-log[10](italic(P))), y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey93", colour = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

.tissue_compare_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                    sig_basis = "fdr",
                                    sig_threshold = 0.05,
                                    font_size = 12,
                                    point_size = 4) {
  packed <- .tissue_compare_frame(long, gwas_vec, only_recurrent, k_min,
                                    sig_basis, sig_threshold)
  if (is.null(packed)) return(NULL)
  slice <- packed$slice

  # Point-style heatmap: filled coloured circle by -log10(P). Overlays:
  # solid black ring for nominal-sig (P < 0.05); solid black square for
  # FDR-sig; inner solid dot for tissues retained after the conditional
  # analysis. Data circle is redrawn last so the fill sits inside the
  # overlay markers.
  slice[, minus_log10 := pmin(-log10(p), 12)]
  slice[, nom_sig := !is.na(p) & p < 0.05]
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]
  slice[, retained := !is.na(evidence) & evidence]

  scale_lab <- "-log10(P)"

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
  # Retained tissues get a small inner white dot with black outline
  # (matching the facet plot) so it reads on both dark-teal FDR-sig
  # fills and pale not-sig fills.
  ret <- slice[slice$retained, ]
  if (nrow(ret) > 0) {
    gg <- gg + ggplot2::geom_point(data = ret, shape = 21, fill = "white",
                                     colour = "black", stroke = 0.5,
                                     size = point_size * 0.45)
  }

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

tissue_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                    comparison_long) {
  moduleServer(id, function(input, output, session) {

    # Server-side k slider so max always tracks the current selection.
    output$k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
    })

    # Include-GWAS picker — restrict which bundle GWAS appear in the plot.
    # Rendered via renderUI (choices baked in at creation) because this
    # whole UI lives inside a renderUI swap in mod_enrichment.
    output$gwas_pick_ui <- renderUI({
      choices <- selected_gwas_multi()
      req(length(choices) >= 1)
      cur <- isolate(input$gwas_pick)
      keep <- if (is.null(cur)) choices else intersect(cur, choices)
      if (length(keep) == 0) keep <- choices
      selectInput(session$ns("gwas_pick"), "Include GWAS:",
                   choices = choices, selected = keep, multiple = TRUE)
    })

    # Include-tissues picker — restrict the tissue set shown in the plot.
    # Defaults to all tissues present in the current MAGMA-tissue slice.
    output$tissue_pick_ui <- renderUI({
      long <- comparison_long()
      req(!is.null(long))
      tissues <- sort(unique(as.character(
        long[method == "MAGMA-tissue", entity_id])))
      req(length(tissues) >= 1)
      cur <- isolate(input$tissue_pick)
      keep <- if (is.null(cur)) tissues else intersect(cur, tissues)
      if (length(keep) == 0) keep <- tissues
      selectInput(session$ns("tissue_pick"), "Include tissues:",
                   choices = tissues, selected = keep, multiple = TRUE)
    })

    # Effective GWAS vector = intersection of sidebar's selected_gwas_multi()
    # and the per-view `gwas_pick`, then sorted per `gwas_sort`.
    gwas_vec_r <- reactive({
      base <- selected_gwas_multi()
      pick <- input$gwas_pick
      chosen <- if (is.null(pick) || length(pick) == 0) base else intersect(base, pick)
      if (length(chosen) == 0) chosen <- base
      order_gwas(chosen,
                  if (is.null(input$gwas_sort)) "as_selected" else input$gwas_sort,
                  NULL)
    })

    # Long slice pre-filtered by the tissue picker. Applies to both
    # facet and heatmap modes.
    tissue_long_filt <- reactive({
      long <- comparison_long()
      req(!is.null(long))
      picked <- input$tissue_pick
      if (is.null(picked) || length(picked) == 0) return(long)
      # Only restrict the tissue-method rows; other entity types
      # (drug/gene/atc/etc.) are not consumed by this view anyway.
      long[method != "MAGMA-tissue" | entity_id %in% picked]
    })

    plot_obj <- reactive({
      long <- tissue_long_filt()
      req(!is.null(long))
      if (identical(input$plot_type %||% "facet", "facet")) {
        .tissue_compare_facet_ggplot(
          long           = long,
          gwas_vec       = gwas_vec_r(),
          only_recurrent = isTRUE(input$only_recurrent),
          k_min          = if (is.null(input$k_min)) 2L else as.integer(input$k_min),
          sig_threshold  = if (is.null(input$sig_threshold)) 0.05 else as.numeric(input$sig_threshold),
          font_size      = if (is.null(input$plot_font_size)) 12 else as.numeric(input$plot_font_size),
          point_size     = if (is.null(input$plot_point_size)) 3 else as.numeric(input$plot_point_size)
        )
      } else {
        .tissue_compare_ggplot(
          long           = long,
          gwas_vec       = gwas_vec_r(),
          only_recurrent = isTRUE(input$only_recurrent),
          k_min          = if (is.null(input$k_min)) 2L else as.integer(input$k_min),
          sig_basis      = "fdr",
          sig_threshold  = if (is.null(input$sig_threshold)) 0.05 else as.numeric(input$sig_threshold),
          font_size      = if (is.null(input$plot_font_size)) 12 else as.numeric(input$plot_font_size),
          point_size     = if (is.null(input$plot_point_size)) 4 else as.numeric(input$plot_point_size)
        )
      }
    })

    plot_dims <- reactive({
      # Derive from the actual plot data so filtered / capped rows shrink
      # the panel height. Facet mode: n_cols = 1 (all facets share vertical
      # tissue axis, and widths are set by faceting rather than cell-per-GWAS).
      p <- plot_obj()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_id))
      y_lab  <- if (is.null(p)) NULL else as.character(unique(p$data$entity_id))
      is_facet <- identical(input$plot_type %||% "facet", "facet")
      n_cols_dim <- if (is_facet) max(2L, length(gwas_vec_r())) else length(gwas_vec_r())
      dims <- .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = n_cols_dim,
        font_size = if (is.null(input$plot_font_size)) 12 else as.numeric(input$plot_font_size),
        y_labels  = y_lab
      )
      # Facet mode needs extra width per facet (each panel is a full mini-plot,
      # not a single-cell heatmap column).
      if (is_facet) dims$width <- dims$width * 1.4
      dims
    })

    .sync_dl_dims(session, plot_dims)

    output$tissue_compare_plot_ui <- renderUI({
      if (is.null(plot_obj())) return(.compare_empty_state("tissues"))
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

    tissue_tbl_slice <- reactive({
      long <- tissue_long_filt()
      req(!is.null(long))
      long[method == "MAGMA-tissue" & gwas %in% gwas_vec_r()]
    })

    output$tissue_compare_tbl <- DT::renderDT({
      slice <- tissue_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS     = factor(slice$gwas, levels = gwas_vec_r()),
        Tissue   = slice$entity_id,
        BETA     = .fmt_sig(slice$statistic, 3),
        SE       = .fmt_sig(slice$se, 3),
        P        = .fmt_sig(slice$p, 3),
        `P.FDR`  = .fmt_sig(slice$fdr, 3),
        Retained = ifelse(is.na(slice$evidence), "—",
                          ifelse(slice$evidence, "Yes", "No")),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE,
                    filter = "top", selection = "none",
                    options = list(pageLength = 20, server = TRUE,
                                    order = list(list(4, "asc"))))
    }, server = TRUE)

    observeEvent(input$tissue_compare_tbl_rows_selected, {
      sel <- input$tissue_compare_tbl_rows_selected
      if (length(sel) != 1) return()
      slice <- tissue_tbl_slice()
      if (nrow(slice) < sel) return()
      .show_entity_detail(as.character(slice$entity_id[sel]),
                            "tissue", comparison_long(), gwas_vec_r())
    }, ignoreInit = TRUE)

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
        long <- tissue_long_filt()[method == "MAGMA-tissue" & gwas %in% gwas_vec]
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
                          "fdr" %||% "fdr",
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
      "Cross-GWAS locus-level associations. ",
      tags$b("Each row is one locus, labelled by the nearest gene to that trait's lead SNP."),
      " Cells show the smallest p-value at that locus for the selected ",
      "GWAS (from LD-based clumping or conditional / joint COJO, depending ",
      "on the Method selection). Blank cell = the locus has no clumped or ",
      "independent COJO signal in that GWAS. Default significance threshold ",
      "is P < 5×10⁻⁸ (genome-wide significance)."
    ),
    hr(),
    tags$details(class = "gd-details",
      tags$summary("Filter data"),
      tags$div(class = "gd-details-body",
        fluidRow(
          column(3,
            # Method dropdown is server-rendered so it reflects the
            # methods actually present in the bundle at the moment
            # the compare UI is shown ("COJO" only offered when the
            # pipeline ran COJO; "clump" whenever it has rows).
            uiOutput(ns("method_ui"))
          ),
          column(3,
            numericInput(ns("sig_threshold"), "Significance threshold (P):",
                          value = 5e-8, min = 1e-20, max = 1, step = 1e-8)
          ),
          column(3,
            checkboxInput(ns("only_recurrent"),
                          "Only show loci significant in ≥ k GWAS",
                          value = TRUE),
            conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("only_recurrent")),
              uiOutput(ns("k_slider_ui"))
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
    gd_legend(list(
      "Ring around a circle"         = "Genome-wide significant (P < 5×10⁻⁸).",
      "Black square around a circle" = "P < chosen significance threshold.",
      "Blank cell"                   = "Locus has no clumped / independent signal in that GWAS."
    ), heading = "How to read this heatmap"),
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
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

locus_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                    comparison_long) {
  moduleServer(id, function(input, output, session) {

    # Method dropdown built from methods actually present in the bundle.
    # Uses renderUI (not selectInput + updateSelectInput) so the
    # dropdown is guaranteed to reflect the current data at the moment
    # the compare UI first renders — earlier we saw COJO leak into the
    # dropdown when the update message raced against the widget being
    # created.
    output$method_ui <- renderUI({
      long <- comparison_long()
      present <- unique(long[entity_type == "locus", method])
      choices <- intersect(.locus_methods, present)
      if (length(choices) == 0) choices <- .locus_methods
      cur <- isolate(input$method)
      sel <- if (!is.null(cur) && cur %in% choices) cur else choices[1L]
      selectInput(session$ns("method"), "Method:",
                   choices = setNames(choices, choices),
                   selected = sel)
    })

    output$k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
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

    .sync_dl_dims(session, plot_dims, "locus")

    output$locus_compare_plot_ui <- renderUI({
      if (is.null(plot_obj())) return(.compare_empty_state("loci"))
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

    locus_tbl_slice <- reactive({
      req(comparison_long())
      picked_method <- input$method %||% "clump"
      comparison_long()[method == picked_method & entity_type == "locus" &
                          gwas %in% gwas_vec_r()]
    })

    output$locus_compare_tbl <- DT::renderDT({
      slice <- locus_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS       = factor(slice$gwas, levels = gwas_vec_r()),
        Locus      = slice$entity_id,
        BETA       = .fmt_sig(slice$statistic, 3),
        SE         = .fmt_sig(slice$se, 3),
        P          = .fmt_sig(slice$p, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(4, "asc"))))
    }, server = TRUE)

    observeEvent(input$locus_compare_tbl_rows_selected, {
      sel <- input$locus_compare_tbl_rows_selected
      if (length(sel) != 1) return()
      slice <- locus_tbl_slice()
      if (nrow(slice) < sel) return()
      .show_entity_detail(as.character(slice$entity_id[sel]),
                            "locus", comparison_long(), gwas_vec_r())
    }, ignoreInit = TRUE)

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
            # Rendered server-side so choices reflect the methods
            # actually present in comparison_long (skipping methods the
            # pipeline didn't run — e.g. no SMR-protein option if no
            # SMR-protein panels were selected).
            uiOutput(ns("method_ui")),
            selectInput(ns("panel"), "Panel:",
                        choices = c("Best-per-cell (min P)" = "__best__"),
                        selected = "__best__"),
            checkboxInput(ns("evidence_required"),
                          "Require colocalisation / HEIDI evidence",
                          value = FALSE)
          ),
          column(3,
            numericInput(ns("sig_threshold"), "Significance threshold:",
                          value = 0.05, min = 1e-12, max = 1, step = 0.01)
          ),
          column(3,
            checkboxInput(ns("only_recurrent"),
                          "Only show genes significant in ≥ k GWAS",
                          value = TRUE),
            conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("only_recurrent")),
              uiOutput(ns("k_slider_ui"))
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
    gd_legend(list(
      "Ring around a circle"         = "Nominal-significant (P < 0.05).",
      "Black square around a circle" = "FDR-significant (P.FDR < 0.05).",
      "Blank cell"                   = "Gene has no result for this method in that GWAS."
    ), heading = "How to read this heatmap"),
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
  basis_col <- "fdr"
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

  scale_lab <- "-log10(FDR)"

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
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

gene_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                  comparison_long) {
  moduleServer(id, function(input, output, session) {

    output$k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
    })

    # Method dropdown rendered from methods actually present in
    # comparison_long — build_comparison_long only emits rows for
    # analyses that ran, so this filters to what the pipeline produced.
    output$method_ui <- renderUI({
      long <- comparison_long()
      req(long)
      have <- intersect(.gene_methods,
                          unique(as.character(long[entity_type == "gene", method])))
      req(length(have) >= 1)
      cur <- isolate(input$method)
      keep <- if (!is.null(cur) && cur %in% have) cur else have[1L]
      selectInput(session$ns("method"), "Method:",
                   choices = have, selected = keep)
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
        sig_basis         = "fdr",
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

    .sync_dl_dims(session, plot_dims, "gene")

    output$gene_compare_plot_ui <- renderUI({
      if (is.null(plot_obj())) return(.compare_empty_state("genes"))
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

    gene_tbl_slice <- reactive({
      req(comparison_long())
      picked_method <- input$method %||% "MAGMA-gene"
      slice <- comparison_long()[method == picked_method & gwas %in% gwas_vec_r()]
      picked_panel <- input$panel %||% "__best__"
      if (!identical(picked_panel, "__best__") && nzchar(picked_panel)) {
        slice <- slice[panel == picked_panel]
      }
      slice
    })

    output$gene_compare_tbl <- DT::renderDT({
      slice <- gene_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      method_pick <- input$method %||% "MAGMA-gene"
      # Method-specific column set — MAGMA has no panel / Z / SE /
      # colocalisation, so those columns are irrelevant and would fill
      # with NAs. FUSION methods report Z / P / P.FDR / COLOC; SMR
      # methods report BETA / SE / P / P.FDR / HEIDI-passed. Column
      # names and the "evidence" flag are renamed to reflect what the
      # column actually means for the chosen method.
      is_fusion <- method_pick %in% c("TWAS-FUSION", "PWAS-FUSION")
      is_smr    <- method_pick %in% c("SMR-expression", "SMR-protein")
      is_magma  <- method_pick == "MAGMA-gene"

      out <- data.frame(
        GWAS  = factor(slice$gwas, levels = gwas_vec_r()),
        Gene  = slice$entity_id,
        stringsAsFactors = FALSE
      )
      if (!is_magma) out$Panel <- ifelse(is.na(slice$panel), "—", slice$panel)
      if (is_fusion) {
        out$Z <- .fmt_sig(slice$statistic, 3)
      } else if (is_smr) {
        out$BETA <- .fmt_sig(slice$statistic, 3)
        out$SE   <- .fmt_sig(slice$se, 3)
      }
      out$P       <- .fmt_sig(slice$p, 3)
      out$`P.FDR` <- .fmt_sig(slice$fdr, 3)
      if (is_fusion) {
        out$Colocalised <- ifelse(is.na(slice$evidence), "—",
                                    ifelse(slice$evidence, "Yes", "No"))
      } else if (is_smr) {
        out$`HEIDI supported` <- ifelse(is.na(slice$evidence), "—",
                                          ifelse(slice$evidence, "Yes", "No"))
      }
      # Sort by P.FDR by default. Column index depends on how many
      # extra columns are present.
      pfdr_col <- which(names(out) == "P.FDR") - 1L  # 0-based for DT
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(pfdr_col, "asc"))))
    }, server = TRUE)

    observeEvent(input$gene_compare_tbl_rows_selected, {
      sel <- input$gene_compare_tbl_rows_selected
      if (length(sel) != 1) return()
      slice <- gene_tbl_slice()
      if (nrow(slice) < sel) return()
      .show_entity_detail(as.character(slice$entity_id[sel]),
                            "gene", comparison_long(), gwas_vec_r())
    }, ignoreInit = TRUE)

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
                           "fdr" %||% "fdr",
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
              numericInput(ns("magma_sig_threshold"), "Significance threshold:",
                            value = 0.05, min = 1e-12, max = 1, step = 0.01)
            ),
            column(3,
              checkboxInput(ns("magma_only_recurrent"),
                            "Only show classes significant in ≥ k GWAS",
                            value = TRUE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("magma_only_recurrent")),
                uiOutput(ns("magma_k_slider_ui"))
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
      gd_legend(list(
        "Ring around a circle"         = "Nominal-significant (P < 0.05).",
        "Black square around a circle" = "FDR-significant (P.FDR < 0.05)."
      ), heading = "How to read this heatmap"),
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
              numericInput(ns("gsea_sig_threshold"), "Significance threshold:",
                            value = 0.05, min = 1e-12, max = 1, step = 0.01)
            ),
            column(3,
              checkboxInput(ns("gsea_only_recurrent"),
                            "Only show classes significant in ≥ k GWAS",
                            value = TRUE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("gsea_only_recurrent")),
                uiOutput(ns("gsea_k_slider_ui"))
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
        "Ring around a circle"         = "Nominal-significant (P < 0.05).",
        "Black square around a circle" = "FDR-significant (P.FDR < 0.05).",
        "Blank cell"                   = "ATC class was not tested in that GWAS.",
        "Blue"                         = "Direction 'Matches disease' — drug class shares the trait's TWAS signature.",
        "Red"                          = "Direction 'Opposes disease' — drug class counteracts the trait's TWAS signature."
      ), heading = "How to read this heatmap"),
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
  slice[, minus_log10_fdr := -log10(fdr)]
  slice[, nom_sig := !is.na(p)   & p   < 0.05]
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]

  # Drop untested cells so no point is drawn for (ATC, GWAS) pairs the
  # method never scored — previously they rendered as grey via
  # na.value, indistinguishable from a near-zero score.
  tested <- slice[!is.na(minus_log10_fdr)]
  if (nrow(tested) == 0) return(NULL)

  # Data-driven fill range so the ramp matches the observed values.
  max_val <- suppressWarnings(max(tested$minus_log10_fdr, na.rm = TRUE))
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1
  fill_lim <- c(0, max_val)

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = minus_log10_fdr), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(tested, ggplot2::aes(x = gwas, y = entity_label)) +
    base_layer +
    ggplot2::scale_fill_gradient(low = "#e6f4f1", high = .gd_teal,
                                   name = "-log10(FDR)",
                                   na.value = "transparent",
                                   limits = fill_lim)

  gg <- .compare_add_sig_overlays(gg, tested, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
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

  # Data-driven symmetric limits so the colour ramp matches the observed
  # range (fixed c(-12, 12) previously washed out small signed scores).
  max_abs <- suppressWarnings(max(abs(tested_dat$signed_score), na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  fill_lim <- c(-max_abs, max_abs)

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = signed_score), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(tested_dat, ggplot2::aes(x = gwas, y = entity_label)) +
    base_layer +
    ggplot2::scale_fill_gradient2(
      low = .gd_red, mid = "white", high = .gd_blue, midpoint = 0,
      limits = fill_lim, na.value = "transparent",
      name = "signed -log10(FDR)  (+ matches / − opposes)"
    )

  gg <- .compare_add_sig_overlays(gg, tested_dat, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
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

    # Server-side k sliders so max always tracks the current selection.
    output$magma_k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$magma_k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("magma_k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
    })
    output$gsea_k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$gsea_k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("gsea_k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
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
        sig_basis      = "fdr",
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

    .sync_dl_dims(session, magma_dims, "magma")

    output$atc_magma_plot_ui <- renderUI({
      if (is.null(magma_plot())) return(.compare_empty_state("ATC classes"))
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

    atc_magma_tbl_slice <- reactive({
      req(comparison_long())
      comparison_long()[method == "MAGMA-ATC" & gwas %in% magma_gwas_vec()]
    })

    output$atc_magma_tbl <- DT::renderDT({
      slice <- atc_magma_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS   = factor(slice$gwas, levels = magma_gwas_vec()),
        `ATC Code`  = slice$entity_id,
        Class       = slice$entity_label,
        `N Drugs`  = .fmt_int(slice$n_units),
        P           = .fmt_sig(slice$p, 3),
        `P.FDR`     = .fmt_sig(slice$fdr, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
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
                           "fdr" %||% "fdr",
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
        sig_basis      = "fdr",
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

    .sync_dl_dims(session, gsea_dims, "gsea")

    output$atc_gsea_plot_ui <- renderUI({
      if (is.null(gsea_plot())) return(.compare_empty_state("ATC classes"))
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

    atc_gsea_tbl_slice <- reactive({
      req(comparison_long())
      comparison_long()[method == "TWAS-GSEA-ATC" & gwas %in% gsea_gwas_vec()]
    })

    output$atc_gsea_tbl <- DT::renderDT({
      slice <- atc_gsea_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS       = factor(slice$gwas, levels = gsea_gwas_vec()),
        `ATC Code` = slice$entity_id,
        Class      = slice$entity_label,
        Panel      = slice$panel,
        Estimate   = .fmt_sig(slice$statistic, 3),
        Direction  = slice$direction,
        Reversal_Z = .fmt_sig(slice$reversal_z, 3),
        P          = .fmt_sig(slice$p, 3),
        `P.FDR`    = .fmt_sig(slice$fdr, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
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
                           "fdr" %||% "fdr",
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

########################################
# DRUG COMPARE
########################################
#
# Two sub-tabs (MAGMA and TWAS-GSEA) mirroring the ATC compare view but keyed
# on drug Name. Because the neuropsych bundle carries ~1500 drugs per GWAS,
# the default is to show all drugs with the row cap (50) trimming the
# heatmap; users toggle "only recurrent" to filter to drugs sig in >= k GWAS.

drug_compare_ui <- function(ns) {
  tabsetPanel(
    tabPanel("MAGMA", br(),
      p("Cross-GWAS MAGMA drug-level enrichment. Each row is one drug; ",
        "cell colour is -log10(FDR); cells with FDR < 0.05 are outlined ",
        "with a black square."),
      hr(),
      tags$details(class = "gd-details",
        tags$summary("Filter data"),
        tags$div(class = "gd-details-body",
          fluidRow(
            column(3,
              numericInput(ns("magma_sig_threshold"), "Significance threshold:",
                            value = 0.05, min = 1e-12, max = 1, step = 0.01)
            ),
            column(3,
              checkboxInput(ns("magma_only_recurrent"),
                            "Only show drugs significant in ≥ k GWAS",
                            value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("magma_only_recurrent")),
                uiOutput(ns("magma_k_slider_ui"))
              ),
              numericInput(ns("magma_row_cap"), "Rows shown in heatmap:",
                            value = 50, min = 5, max = 500, step = 5)
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
        uiOutput(ns("drug_magma_plot_ui"))
      ),
      gd_legend(list(
        "Ring around a circle"         = "Nominal-significant (P < 0.05).",
        "Black square around a circle" = "FDR-significant (P.FDR < 0.05)."
      ), heading = "How to read this heatmap"),
      br(),
      tags$div(style = "max-width: 1100px;",
        h4("Underlying data"),
        DT::DTOutput(ns("drug_magma_tbl"))
      )
    ),
    tabPanel("TWAS-GSEA", br(),
      p("Cross-GWAS TWAS-GSEA drug-level enrichment. Cell colour is signed ",
        "by direction of effect (blue = matches disease signature, red = ",
        "opposes disease signature). Intensity is -log10(FDR). Blank cell = ",
        "drug was not tested in that GWAS."),
      hr(),
      tags$details(class = "gd-details",
        tags$summary("Filter data"),
        tags$div(class = "gd-details-body",
          fluidRow(
            column(3,
              numericInput(ns("gsea_sig_threshold"), "Significance threshold:",
                            value = 0.05, min = 1e-12, max = 1, step = 0.01)
            ),
            column(3,
              checkboxInput(ns("gsea_only_recurrent"),
                            "Only show drugs significant in ≥ k GWAS",
                            value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("gsea_only_recurrent")),
                uiOutput(ns("gsea_k_slider_ui"))
              ),
              numericInput(ns("gsea_row_cap"), "Rows shown in heatmap:",
                            value = 50, min = 5, max = 500, step = 5)
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
        uiOutput(ns("drug_gsea_plot_ui"))
      ),
      gd_legend(list(
        "Ring around a circle"         = "Nominal-significant (P < 0.05).",
        "Black square around a circle" = "FDR-significant (P.FDR < 0.05).",
        "Blank cell"                   = "Drug was not tested in that GWAS.",
        "Blue"                         = "Direction 'Matches disease' — drug's TWAS signature matches the trait's.",
        "Red"                          = "Direction 'Opposes disease' — drug's TWAS signature counteracts the trait's."
      ), heading = "How to read this heatmap"),
      br(),
      tags$div(style = "max-width: 1100px;",
        h4("Underlying data"),
        DT::DTOutput(ns("drug_gsea_tbl"))
      )
    )
  )
}

# Reuses the ATC row-order helper (.atc_row_order) since the row-sorting
# logic is method-agnostic.

.drug_magma_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                  sig_basis = "fdr", sig_threshold = 0.05,
                                  row_cap = 50, font_size = 11,
                                  point_size = 4) {
  slice <- long[method == "MAGMA-drug" & gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)
  rec <- .atc_row_order(slice, sig_basis, sig_threshold)
  if (isTRUE(only_recurrent)) rec <- rec[n_sig >= k_min]
  if (nrow(rec) == 0) return(NULL)
  cap <- max(1L, min(as.integer(row_cap), nrow(rec)))
  rec <- rec[seq_len(cap)]
  slice <- slice[entity_id %in% rec$entity_id]

  slice[, entity_id := factor(entity_id, levels = rev(rec$entity_id))]
  slice[, gwas := factor(gwas, levels = gwas_vec)]
  slice[, minus_log10_fdr := -log10(fdr)]
  slice[, nom_sig := !is.na(p)   & p   < 0.05]
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]

  # Drop untested cells so no point is drawn.
  tested <- slice[!is.na(minus_log10_fdr)]
  if (nrow(tested) == 0) return(NULL)

  max_val <- suppressWarnings(max(tested$minus_log10_fdr, na.rm = TRUE))
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1
  fill_lim <- c(0, max_val)

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = minus_log10_fdr), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(tested, ggplot2::aes(x = gwas, y = entity_id)) +
    base_layer +
    ggplot2::scale_fill_gradient(low = "#e6f4f1", high = .gd_teal,
                                   name = "-log10(FDR)",
                                   na.value = "transparent",
                                   limits = fill_lim)

  gg <- .compare_add_sig_overlays(gg, tested, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

.drug_gsea_frame <- function(long, gwas_vec, only_recurrent, k_min,
                                panel_pick = "__best__",
                                sig_basis = "fdr", sig_threshold = 0.05,
                                row_cap = 50) {
  slice <- long[method == "TWAS-GSEA-drug" & gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)

  if (identical(panel_pick, "__best__") || is.null(panel_pick) ||
      !nzchar(panel_pick)) {
    slice <- pick_best_per_cell(slice, c("gwas", "entity_id"))
  } else {
    slice <- slice[panel == panel_pick]
  }
  if (nrow(slice) == 0) return(NULL)

  rec <- .atc_row_order(slice, sig_basis, sig_threshold)
  if (isTRUE(only_recurrent)) rec <- rec[n_sig >= k_min]
  if (nrow(rec) == 0) return(NULL)
  cap <- max(1L, min(as.integer(row_cap), nrow(rec)))
  rec <- rec[seq_len(cap)]
  slice <- slice[entity_id %in% rec$entity_id]

  # Blank cells for "not tested" — we simply leave them out; the heatmap
  # renders no point where there's no data.
  slice[, gwas := factor(gwas, levels = gwas_vec)]
  slice[, entity_id := factor(entity_id, levels = rev(rec$entity_id))]
  slice[, nom_sig := !is.na(p)   & p   < 0.05]
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]
  neg_log_fdr <- pmin(-log10(slice$fdr), 12)
  slice[, signed_score := ifelse(
    is.na(direction), NA_real_,
    ifelse(direction == "Matches disease",  neg_log_fdr,
    ifelse(direction == "Opposes disease", -neg_log_fdr, NA_real_)))]
  slice
}

.drug_gsea_ggplot <- function(long, gwas_vec, only_recurrent, k_min,
                                 panel_pick, sig_basis, sig_threshold,
                                 row_cap = 50, font_size = 11,
                                 point_size = 4) {
  slice <- .drug_gsea_frame(long, gwas_vec, only_recurrent, k_min,
                              panel_pick, sig_basis, sig_threshold, row_cap)
  if (is.null(slice)) return(NULL)

  # Untested (drug, GWAS) cells → NA signed_score → drop entirely so no
  # point is drawn. Previously we set na.value on the gradient scale,
  # which rendered grey circles for untested cells and made "not tested"
  # indistinguishable from "near-zero score".
  tested <- slice[!is.na(signed_score)]
  if (nrow(tested) == 0) return(NULL)

  # Data-driven symmetric limits so the colour ramp matches the observed
  # range. Fixed c(-12, 12) previously left almost every cell near the
  # scale's midpoint when actual scores were small.
  max_abs <- suppressWarnings(max(abs(tested$signed_score), na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  fill_lim <- c(-max_abs, max_abs)

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = signed_score), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(tested, ggplot2::aes(x = gwas, y = entity_id)) +
    base_layer +
    ggplot2::scale_fill_gradient2(
      low = .gd_red, mid = "white", high = .gd_blue, midpoint = 0,
      limits = fill_lim, na.value = "transparent",
      name = "signed -log10(FDR)  (+ matches / − opposes)"
    )

  gg <- .compare_add_sig_overlays(gg, tested, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

drug_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                   comparison_long) {
  moduleServer(id, function(input, output, session) {

    # Refresh the TWAS-GSEA panel dropdown with panels actually present.
    observeEvent(comparison_long(), {
      long <- comparison_long()
      panels <- unique(long[method == "TWAS-GSEA-drug", panel])
      panels <- panels[!is.na(panels)]
      choices <- c("Best-per-cell (min P)" = "__best__",
                    setNames(panels, panels))
      updateSelectInput(session, "gsea_panel",
                         choices = choices, selected = "__best__")
    })

    output$magma_k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$magma_k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("magma_k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
    })
    output$gsea_k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$gsea_k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("gsea_k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
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
      .drug_magma_ggplot(
        long           = comparison_long(),
        gwas_vec       = magma_gwas_vec(),
        only_recurrent = isTRUE(input$magma_only_recurrent),
        k_min          = if (is.null(input$magma_k_min)) 2L else as.integer(input$magma_k_min),
        sig_basis      = "fdr",
        sig_threshold  = if (is.null(input$magma_sig_threshold)) 0.05 else as.numeric(input$magma_sig_threshold),
        row_cap        = if (is.null(input$magma_row_cap)) 50L else as.integer(input$magma_row_cap),
        font_size      = if (is.null(input$magma_plot_font_size)) 11 else as.numeric(input$magma_plot_font_size),
        point_size     = if (is.null(input$magma_plot_point_size)) 4 else as.numeric(input$magma_plot_point_size)
      )
    })

    magma_dims <- reactive({
      p <- magma_plot()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_id))
      .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = length(magma_gwas_vec()),
        font_size = if (is.null(input$magma_plot_font_size)) 11 else as.numeric(input$magma_plot_font_size),
        y_labels  = if (is.null(p)) NULL else as.character(unique(p$data$entity_id))
      )
    })

    .sync_dl_dims(session, magma_dims, "magma")

    output$drug_magma_plot_ui <- renderUI({
      if (is.null(magma_plot())) return(.compare_empty_state("drugs"))
      dims <- magma_dims()
      plotOutput(session$ns("drug_magma_plot"),
                  height = dims$height, width = dims$width)
    })

    output$drug_magma_plot <- renderPlot({
      p <- magma_plot()
      if (is.null(p)) {
        plot.new(); title("No drugs meet the current filter")
        return(invisible())
      }
      print(p)
    })

    drug_magma_tbl_slice <- reactive({
      req(comparison_long())
      comparison_long()[method == "MAGMA-drug" & gwas %in% magma_gwas_vec()]
    })

    output$drug_magma_tbl <- DT::renderDT({
      slice <- drug_magma_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS   = factor(slice$gwas, levels = magma_gwas_vec()),
        Drug   = slice$entity_id,
        `N Genes` = .fmt_int(slice$n_units),
        BETA   = .fmt_sig(slice$statistic, 3),
        SE     = .fmt_sig(slice$se, 3),
        P      = .fmt_sig(slice$p, 3),
        `P.FDR` = .fmt_sig(slice$fdr, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(6, "asc"))))
    }, server = TRUE)

    output$magma_download_plot <- downloadHandler(
      filename = function() sprintf("drug_magma_compare_%s.%s",
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
      filename = function() sprintf("drug_magma_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- magma_gwas_vec()
        long <- comparison_long()[method == "MAGMA-drug" & gwas %in% gwas_vec]
        wide <- pivot_matrix(long, "fdr", gwas_vec)
        header <- sprintf("# GenoDisc drug MAGMA compare CSV | %s | sig_basis=%s threshold=%g k_min=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           "fdr" %||% "fdr",
                           input$magma_sig_threshold %||% 0.05,
                           input$magma_k_min %||% 2L)
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: P.FDR", con)
        out <- data.frame(Drug = rownames(wide),
                          format(wide, scientific = TRUE, digits = 3),
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )

    # ---------- TWAS-GSEA plot / table ----------

    gsea_plot <- reactive({
      req(comparison_long())
      .drug_gsea_ggplot(
        long           = comparison_long(),
        gwas_vec       = gsea_gwas_vec(),
        only_recurrent = isTRUE(input$gsea_only_recurrent),
        k_min          = if (is.null(input$gsea_k_min)) 2L else as.integer(input$gsea_k_min),
        panel_pick     = input$gsea_panel,
        sig_basis      = "fdr",
        sig_threshold  = if (is.null(input$gsea_sig_threshold)) 0.05 else as.numeric(input$gsea_sig_threshold),
        row_cap        = if (is.null(input$gsea_row_cap)) 50L else as.integer(input$gsea_row_cap),
        font_size      = if (is.null(input$gsea_plot_font_size)) 11 else as.numeric(input$gsea_plot_font_size),
        point_size     = if (is.null(input$gsea_plot_point_size)) 4 else as.numeric(input$gsea_plot_point_size)
      )
    })

    gsea_dims <- reactive({
      p <- gsea_plot()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_id))
      .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = length(gsea_gwas_vec()),
        font_size = if (is.null(input$gsea_plot_font_size)) 11 else as.numeric(input$gsea_plot_font_size),
        y_labels  = if (is.null(p)) NULL else as.character(unique(p$data$entity_id))
      )
    })

    .sync_dl_dims(session, gsea_dims, "gsea")

    output$drug_gsea_plot_ui <- renderUI({
      if (is.null(gsea_plot())) return(.compare_empty_state("drugs"))
      dims <- gsea_dims()
      plotOutput(session$ns("drug_gsea_plot"),
                  height = dims$height, width = dims$width)
    })

    output$drug_gsea_plot <- renderPlot({
      p <- gsea_plot()
      if (is.null(p)) {
        plot.new(); title("No drugs meet the current filter")
        return(invisible())
      }
      print(p)
    })

    drug_gsea_tbl_slice <- reactive({
      req(comparison_long())
      comparison_long()[method == "TWAS-GSEA-drug" & gwas %in% gsea_gwas_vec()]
    })

    output$drug_gsea_tbl <- DT::renderDT({
      slice <- drug_gsea_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS       = factor(slice$gwas, levels = gsea_gwas_vec()),
        Drug       = slice$entity_id,
        Panel      = slice$panel,
        `N Genes`  = .fmt_int(slice$n_units),
        Estimate   = .fmt_sig(slice$statistic, 3),
        Direction  = slice$direction,
        Reversal_Z = .fmt_sig(slice$reversal_z, 3),
        P          = .fmt_sig(slice$p, 3),
        `P.FDR`    = .fmt_sig(slice$fdr, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(8, "asc"))))
    }, server = TRUE)

    output$gsea_download_plot <- downloadHandler(
      filename = function() sprintf("drug_gsea_compare_%s.%s",
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
      filename = function() sprintf("drug_gsea_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- gsea_gwas_vec()
        long <- comparison_long()[method == "TWAS-GSEA-drug" & gwas %in% gwas_vec]
        best <- pick_best_per_cell(long, c("gwas", "entity_id"))
        best[, signed := ifelse(!is.na(direction) & direction == "Opposes disease",
                                 -(-log10(fdr)), -log10(fdr))]
        wide <- pivot_matrix(best, "signed", gwas_vec)
        display <- ifelse(is.na(wide), "NT", sprintf("%+.2f", wide))
        display <- matrix(display, nrow = nrow(wide), dimnames = dimnames(wide))
        header <- sprintf("# GenoDisc drug TWAS-GSEA compare CSV | %s | sig_basis=%s threshold=%g k_min=%s panel=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           "fdr" %||% "fdr",
                           input$gsea_sig_threshold %||% 0.05,
                           input$gsea_k_min %||% 2L,
                           input$gsea_panel %||% "__best__")
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: signed -log10(FDR): + = matches, - = opposes. NT = not tested.", con)
        out <- data.frame(Drug = rownames(display), display,
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )
  })
}

########################################
# GENCOR COMPARE (Bivariate LDSC)
########################################
#
# Cross-trait genetic-correlation heatmap. Reads from each selected GWAS's
# gwas_qc$ldsc_gencor_dat block (built by build_gencor_long). Bundles that
# didn't run bivariate LDSC get a "no data" placeholder instead of an empty
# plot.

gencor_compare_ui <- function(ns) {
  tagList(
    br(),
    p(
      "Cross-GWAS view of bivariate-LDSC genetic correlations (rG) with ",
      "external reference traits. Default plot lays out each bundle GWAS in ",
      "its own column facet as a forest plot (rG ± 95% CI, one row per ",
      "reference trait, coloured by nominal / FDR significance). Switch to ",
      "the heatmap when comparing many GWAS."
    ),
    hr(),
    tags$details(class = "gd-details",
      tags$summary("Filter data"),
      tags$div(class = "gd-details-body",
        fluidRow(
          column(3,
            radioButtons(ns("plot_type"), "Plot type:",
                          choices = c("Facet by GWAS" = "facet",
                                       "Heatmap" = "heatmap"),
                          selected = "facet", inline = TRUE),
            # `gwas_pick`, `ref_pick`, `row_facet` are rendered by
            # `renderUI` in the server so their choices are baked in at
            # widget-creation time. `updateSelectInput` sent before the
            # widget exists on the client is silently dropped (this UI
            # lives inside a renderUI swap in mod_h2_gencor.R), so
            # option-population was breaking before this pattern.
            uiOutput(ns("gwas_pick_ui")),
            selectInput(ns("gwas_sort"), "Order GWAS by:",
                         choices = .compare_gwas_sort_choices,
                         selected = "as_selected")
          ),
          column(3,
            uiOutput(ns("ref_pick_ui")),
            selectInput(ns("trait_order"), "Order traits by:",
                         choices = c("Category (then label)" = "category",
                                      "Alphabetical" = "alphabetical",
                                      "|rG| max across GWAS" = "rg_max",
                                      "Significance (-log10 P) max across GWAS" = "sig_max"),
                         selected = "category"),
            uiOutput(ns("row_facet_ui"))
          ),
          column(3,
            sliderInput(ns("plot_font_size"), "Font size (pt):",
                         min = 8, max = 20, value = 11, step = 1),
            sliderInput(ns("plot_point_size"), "Point size:",
                         min = 2, max = 14, value = 3, step = 1)
          ),
          .dl_and_download_column(ns, "gencor", default_w = 12, default_h = 10)
        )
      )
    ),
    br(),
    tags$div(style = "max-width: 1100px; overflow-x: auto;",
      uiOutput(ns("gencor_plot_ui"))
    ),
    gd_legend(list(
      "Facet mode"                   = "One panel per GWAS; each row is a reference trait; horizontal bar shows rG ± 95% CI.",
      "Red square (facet)"           = "FDR-significant (p.FDR < 0.05).",
      "Orange triangle (facet)"      = "Nominal-significant only (p < 0.05).",
      "Grey circle (facet)"          = "Not significant.",
      "Ring around a circle"         = "Heatmap: nominal-significant (p < 0.05).",
      "Black square around a circle" = "Heatmap: FDR-significant (p.FDR < 0.05).",
      "Blank cell"                   = "Heatmap: rG could not be estimated for this pair."
    ), heading = "How to read this plot"),
    br(),
    tags$div(style = "max-width: 1100px;",
      h4("Underlying data"),
      DT::DTOutput(ns("gencor_tbl"))
    )
  )
}

#' Faceted per-GWAS forest plot for the gencor compare view.
#'
#' Mirrors the single-GWAS build_gencor_plot() but places each selected
#' GWAS in its own column facet. Reference-trait order is shared across
#' facets; when a row-facet column is given the traits also stack into
#' semantic bands (e.g. by trait_category).
#'
#' @param long comparison gencor long tibble.
#' @param gwas_vec Which bundle GWAS to include (facet columns, in order).
#' @param trait_order One of "category", "alphabetical", "rg_max", "sig_max".
#'   `category` sorts traits by trait_category (if present) then label —
#'   matches the pre-existing "group by trait category" behaviour.
#'   `rg_max` orders by max |rg| across the selected GWAS.
#'   `sig_max` orders by max -log10(p) across the selected GWAS.
#' @param row_facet_col Optional column name from `long` to use as a
#'   horizontal row facet (e.g. "trait_category"). NULL to skip row facets.
#' @param show_facet_col_strip TRUE to keep the column-facet strip labels.
#'   FALSE hides them (used when only one GWAS is being plotted).
.gencor_facet_ggplot <- function(long, gwas_vec, font_size = 11, point_size = 3,
                                   trait_order = "category",
                                   row_facet_col = NULL,
                                   show_facet_col_strip = TRUE) {
  if (nrow(long) == 0) return(NULL)
  slice <- long[gwas %in% gwas_vec & !is.na(rg) & !is.na(rg_se)]
  if (nrow(slice) == 0) return(NULL)

  slice[, ci_lo := rg - 1.96 * rg_se]
  slice[, ci_hi := rg + 1.96 * rg_se]
  slice[, tier := data.table::fifelse(
        !is.na(rg_p_fdr) & rg_p_fdr < 0.05, "FDR significant",
        data.table::fifelse(!is.na(rg_p) & rg_p < 0.05, "Nominal (p < 0.05)",
                                                        "Non-significant"))]
  tier_levels <- c("FDR significant", "Nominal (p < 0.05)", "Non-significant")
  slice[, tier := factor(tier, levels = tier_levels)]

  # Row order — shared across all facet columns so recurrent traits line up.
  ref_meta <- unique(slice[, .(ref_label, trait_category)])
  ord <- switch(trait_order %||% "category",
    alphabetical = ref_meta[order(ref_label)],
    rg_max = {
      agg <- slice[, .(v = suppressWarnings(max(abs(rg), na.rm = TRUE))),
                    by = ref_label]
      agg <- agg[order(-v)]
      merge(agg[, .(ref_label)], ref_meta, by = "ref_label",
             sort = FALSE)
    },
    sig_max = {
      agg <- slice[, .(v = suppressWarnings(max(-log10(pmax(rg_p, 1e-300)),
                                                   na.rm = TRUE))),
                    by = ref_label]
      agg <- agg[order(-v)]
      merge(agg[, .(ref_label)], ref_meta, by = "ref_label",
             sort = FALSE)
    },
    # default = "category"
    ref_meta[order(is.na(trait_category), trait_category, ref_label)]
  )
  # First factor level renders at the BOTTOM of the plot, so reverse the
  # order for a top-first read.
  slice[, ref_label := factor(ref_label, levels = rev(ord$ref_label))]
  slice[, gwas      := factor(gwas, levels = gwas_vec)]

  tier_cols <- c("FDR significant"    = "#d62728",
                 "Nominal (p < 0.05)" = "#ff7f0e",
                 "Non-significant"    = "#7f7f7f")
  tier_shapes <- c("FDR significant"    = 15,
                   "Nominal (p < 0.05)" = 17,
                   "Non-significant"    = 16)

  # Ghost rows keep the legend showing every tier even if some tiers are
  # absent from the current slice. Must carry every faceting variable.
  ghost_cols <- list(
    rg    = rep(slice$rg[1],        length(tier_levels)),
    ci_lo = rep(slice$rg[1],        length(tier_levels)),
    ci_hi = rep(slice$rg[1],        length(tier_levels)),
    ref_label = rep(slice$ref_label[1], length(tier_levels)),
    gwas  = rep(slice$gwas[1],      length(tier_levels)),
    tier  = factor(tier_levels, levels = tier_levels)
  )
  if (!is.null(row_facet_col) && row_facet_col %in% names(slice)) {
    ghost_cols[[row_facet_col]] <- rep(slice[[row_facet_col]][1],
                                          length(tier_levels))
  }
  ghost <- data.table::as.data.table(ghost_cols)

  gg <- ggplot2::ggplot(slice,
                          ggplot2::aes(x = rg, y = ref_label,
                                        colour = tier, shape = tier)) +
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
    ggplot2::labs(x = expression("Genetic correlation ("*r[g]*")"), y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey93", colour = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )

  use_row_facet <- !is.null(row_facet_col) &&
                    row_facet_col %in% names(slice) &&
                    length(unique(slice[[row_facet_col]])) > 1
  if (use_row_facet) {
    # Bucket NA into "Unknown" so the row-facet strip is never blank.
    vals <- as.character(slice[[row_facet_col]])
    vals[is.na(vals) | vals == ""] <- "Unknown"
    slice[[row_facet_col]] <- vals
    gg$data <- slice
    # Both rows and cols must be `vars()` lists — facet_grid rejects a
    # mix of formula/vars.
    gg <- gg + ggplot2::facet_grid(
      rows = ggplot2::vars(.data[[row_facet_col]]),
      cols = ggplot2::vars(gwas),
      scales = "free_y", space = "free_y"
    ) + ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold")
    )
  } else {
    gg <- gg + ggplot2::facet_grid(cols = ggplot2::vars(gwas),
                                     scales = "free_x")
  }
  if (!isTRUE(show_facet_col_strip)) {
    gg <- gg + ggplot2::theme(strip.text.x = ggplot2::element_blank())
  }
  gg
}

.gencor_ggplot <- function(long, gwas_vec, font_size = 11, point_size = 6,
                             group_by_category = TRUE) {
  if (nrow(long) == 0) return(NULL)
  slice <- long[gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)

  # Row order: group by trait_category (alphabetical) then trait label
  # alphabetical within category, so semantically related traits sit
  # together. If the user turns grouping off, sort by trait label only.
  ref_order <- unique(slice[, .(ref_label, trait_category)])
  if (isTRUE(group_by_category)) {
    ref_order <- ref_order[order(is.na(trait_category), trait_category, ref_label)]
  } else {
    ref_order <- ref_order[order(ref_label)]
  }

  slice[, gwas := factor(gwas, levels = gwas_vec)]
  slice[, ref_label := factor(ref_label, levels = rev(ref_order$ref_label))]
  slice[, nom_sig := !is.na(rg_p)     & rg_p     < 0.05]
  slice[, fdr_sig := !is.na(rg_p_fdr) & rg_p_fdr < 0.05]
  # Clamp rg to [-1.2, 1.2] so occasional out-of-bounds estimates don't
  # push the fill scale off the diverging endpoints.
  slice[, rg_plot := pmax(-1.2, pmin(rg, 1.2))]

  # Drop cells with no rG estimate — they render as truly blank (no point
  # drawn), matching the caption's "Blank: rG not estimated."
  slice <- slice[!is.na(rg)]
  if (nrow(slice) == 0) return(NULL)

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = rg_plot), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(slice, ggplot2::aes(x = gwas, y = ref_label)) +
    base_layer +
    ggplot2::scale_fill_gradient2(
      low = .gd_red, mid = "white", high = .gd_blue, midpoint = 0,
      limits = c(-1.2, 1.2), na.value = .gd_grey,
      name = expression("Genetic correlation ("*r[g]*")")
    )

  gg <- .compare_add_sig_overlays(gg, slice, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(drop = FALSE,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

gencor_compare_server <- function(id, gwas_data, selected_gwas_multi) {
  moduleServer(id, function(input, output, session) {

    gencor_long <- reactive({
      req(gwas_data())
      req(length(selected_gwas_multi()) >= 1)
      build_gencor_long(gwas_data(), selected_gwas_multi())
    })

    # Render the three picker widgets with their choices baked in at
    # creation time. Using renderUI (rather than updateSelectInput on a
    # pre-declared static widget) is essential here because
    # gencor_compare_ui itself lives inside a renderUI swap in
    # mod_h2_gencor.R — a `updateSelectInput` message sent before the
    # target widget appears on the client is silently dropped, which
    # is what caused the pickers to render as empty selects.
    output$gwas_pick_ui <- renderUI({
      choices <- selected_gwas_multi()
      req(length(choices) >= 1)
      cur <- isolate(input$gwas_pick)
      keep <- if (is.null(cur)) choices else intersect(cur, choices)
      if (length(keep) == 0) keep <- choices
      selectInput(session$ns("gwas_pick"), "Include GWAS:",
                   choices = choices, selected = keep, multiple = TRUE)
    })

    output$ref_pick_ui <- renderUI({
      long <- gencor_long()
      if (nrow(long) == 0) return(NULL)
      ref_pairs <- unique(long[, .(ref_name, ref_label)])
      ref_pairs <- ref_pairs[order(ref_label)]
      # selectInput takes a named character vector: names are shown to
      # the user, values are what the input returns.
      opts <- stats::setNames(as.character(ref_pairs$ref_name),
                                as.character(ref_pairs$ref_label))
      cur <- isolate(input$ref_pick)
      keep <- if (is.null(cur)) ref_pairs$ref_name else intersect(cur, ref_pairs$ref_name)
      if (length(keep) == 0) keep <- ref_pairs$ref_name
      selectInput(session$ns("ref_pick"), "Include reference traits:",
                   choices = opts, selected = keep, multiple = TRUE)
    })

    output$row_facet_ui <- renderUI({
      long <- gencor_long()
      if (nrow(long) == 0) return(NULL)
      core <- c("gwas", "ref_name", "ref_label", "rg", "rg_se", "rg_p",
                 "rg_p_fdr", "n_snps", "gcov_int")
      extras <- setdiff(names(long), core)
      # Drop columns that are entirely NA — no sense faceting on them.
      extras <- extras[vapply(extras, function(c) any(!is.na(long[[c]])),
                                logical(1))]
      opts <- c("None" = "__none__", stats::setNames(extras, extras))
      cur <- isolate(input$row_facet)
      default <- if ("trait_category" %in% extras) "trait_category" else "__none__"
      sel <- if (is.null(cur) || !(cur %in% opts)) default else cur
      selectInput(session$ns("row_facet"), "Row facet (metadata):",
                   choices = opts, selected = sel)
    })

    # Effective GWAS vector for plots — intersect user's per-view pick with
    # the sidebar's selected_gwas_multi() to keep both filters authoritative.
    gwas_vec_r <- reactive({
      base <- selected_gwas_multi()
      pick <- input$gwas_pick
      if (is.null(pick) || length(pick) == 0) pick <- base
      chosen <- intersect(base, pick)
      if (length(chosen) == 0) chosen <- base
      order_gwas(chosen,
                  if (is.null(input$gwas_sort)) "as_selected" else input$gwas_sort,
                  NULL)
    })

    # Long slice with reference-trait filter applied.
    gencor_long_filt <- reactive({
      long <- gencor_long()
      req(nrow(long) > 0)
      picked <- input$ref_pick
      if (is.null(picked) || length(picked) == 0) return(long)
      long[ref_name %in% picked]
    })

    plot_obj <- reactive({
      long <- gencor_long_filt()
      req(nrow(long) > 0)
      if (identical(input$plot_type %||% "facet", "facet")) {
        row_facet <- input$row_facet
        if (isTRUE(row_facet == "__none__")) row_facet <- NULL
        .gencor_facet_ggplot(
          long                 = long,
          gwas_vec             = gwas_vec_r(),
          font_size            = if (is.null(input$plot_font_size)) 11 else as.numeric(input$plot_font_size),
          point_size           = if (is.null(input$plot_point_size)) 3 else as.numeric(input$plot_point_size),
          trait_order          = input$trait_order %||% "category",
          row_facet_col        = row_facet,
          # Hide the column strip when only a single GWAS is being plotted
          # (a bare strip labelled with one trait is just noise).
          show_facet_col_strip = length(gwas_vec_r()) > 1L
        )
      } else {
        .gencor_ggplot(
          long              = long,
          gwas_vec          = gwas_vec_r(),
          font_size         = if (is.null(input$plot_font_size)) 11 else as.numeric(input$plot_font_size),
          point_size        = if (is.null(input$plot_point_size)) 6 else as.numeric(input$plot_point_size),
          group_by_category = identical(input$trait_order %||% "category", "category")
        )
      }
    })

    plot_dims <- reactive({
      long <- gencor_long_filt()
      ref_labels <- unique(long[gwas %in% gwas_vec_r(), ref_label])
      is_facet <- identical(input$plot_type %||% "facet", "facet")
      n_cols_dim <- if (is_facet) max(2L, length(gwas_vec_r())) else length(gwas_vec_r())
      dims <- .compare_plot_dims(
        n_rows    = length(ref_labels),
        n_cols    = n_cols_dim,
        font_size = if (is.null(input$plot_font_size)) 11 else as.numeric(input$plot_font_size),
        y_labels  = ref_labels,
        min_width = 400, min_height = 300
      )
      # Facet mode needs extra width per facet (each panel is a mini forest
      # plot with its own x-axis, not a single-cell heatmap column).
      if (is_facet) dims$width <- dims$width * 1.6
      dims
    })

    .sync_dl_dims(session, plot_dims, "gencor")

    output$gencor_plot_ui <- renderUI({
      long <- gencor_long()
      if (nrow(long) == 0) {
        return(tags$div(
          class = "well",
          style = "max-width: 720px; margin-top: 12px;",
          tags$p(
            tags$b("No genetic-correlation data in this bundle."),
            " Bivariate LDSC was not run for the selected GWAS. Enable ",
            tags$code("gencor_gwas_list"),
            " in the pipeline config to populate this view."
          )
        ))
      }
      dims <- plot_dims()
      plotOutput(session$ns("gencor_plot"),
                  height = dims$height, width = dims$width)
    })

    output$gencor_plot <- renderPlot({
      p <- plot_obj()
      if (is.null(p)) { plot.new(); title("No gencor data"); return(invisible()) }
      print(p)
    })

    gencor_tbl_slice <- reactive({
      long <- gencor_long()
      if (nrow(long) == 0) return(long)
      long[gwas %in% gwas_vec_r()]
    })

    output$gencor_tbl <- DT::renderDT({
      slice <- gencor_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS          = factor(slice$gwas, levels = gwas_vec_r()),
        `Ref trait`   = slice$ref_label,
        Category      = ifelse(is.na(slice$trait_category), "—",
                                slice$trait_category),
        rG            = .fmt_sig(slice$rg, 3),
        SE            = .fmt_sig(slice$rg_se, 3),
        P             = .fmt_sig(slice$rg_p, 3),
        `P.FDR`       = .fmt_sig(slice$rg_p_fdr, 3),
        `N SNPs`      = .fmt_int(slice$n_snps),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(5, "asc"))))
    }, server = TRUE)

    observeEvent(input$gencor_tbl_rows_selected, {
      sel <- input$gencor_tbl_rows_selected
      if (length(sel) != 1) return()
      slice <- gencor_tbl_slice()
      if (nrow(slice) < sel) return()
      # Detail modal: show all bundle GWAS' rG values for the picked ref
      # trait (rows = GWAS, so users can compare a single trait across
      # selected GWAS at a glance).
      ref_lab <- as.character(slice$ref_label[sel])
      long <- gencor_long()
      d <- long[ref_label == ref_lab & gwas %in% gwas_vec_r()]
      if (nrow(d) == 0) return()
      d <- d[order(match(gwas, gwas_vec_r()))]
      out <- data.frame(
        GWAS   = factor(d$gwas, levels = gwas_vec_r()),
        rG     = signif(d$rg, 3),
        SE     = signif(d$rg_se, 3),
        P      = signif(d$rg_p, 3),
        `P.FDR` = signif(d$rg_p_fdr, 3),
        `N SNPs` = .fmt_int(d$n_snps),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      showModal(modalDialog(
        title = sprintf("%s (%s)", ref_lab, unique(d$ref_name)[1L]),
        size  = "l", easyClose = TRUE, footer = modalButton("Close"),
        DT::datatable(out, rownames = FALSE,
                      options = list(pageLength = 25, dom = "tip",
                                      order = list(list(0, "asc"))))
      ))
    }, ignoreInit = TRUE)

    output$gencor_download_plot <- downloadHandler(
      filename = function() sprintf("gencor_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$gencor_dl_format),
      content = function(file) {
        p <- plot_obj()
        if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$gencor_dl_format
        w <- input$gencor_dl_width; h <- input$gencor_dl_height
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$gencor_download_csv <- downloadHandler(
      filename = function() sprintf("gencor_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- gwas_vec_r()
        long <- gencor_long()[gwas %in% gwas_vec]
        # Reshape to reference-trait x GWAS wide matrix; retain trait label +
        # category as leading columns so the CSV is self-describing.
        wide <- data.table::dcast(long, ref_name + ref_label + trait_category ~ gwas,
                                    value.var = "rg")
        header <- sprintf("# GenoDisc gencor compare CSV | %s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Rows: external reference traits.  Columns: bundle GWAS.  Cells: rG (LDSC).", con)
        utils::write.csv(wide, con, row.names = FALSE, quote = FALSE, na = "")
      }
    )
  })
}

########################################
# CMAP COMPARE
########################################
#
# Two sub-tabs (Perturbation and MOA) mirroring the ATC / drug compare
# pattern. Both use signed colouring (Direction from CMap output).
# Perturbation rows are keyed on cmap_name; MOA rows on the MOA string.
# For both, best-per-cell reduces over the compound panel key that
# encodes Panel / cell line / dose / time.

cmap_compare_ui <- function(ns) {
  tabsetPanel(
    tabPanel("Perturbation", br(),
      p("Cross-GWAS CMap perturbation compare. Rows are compound / shRNA ",
        "perturbations (cmap_name). Best-per-cell reduction picks the ",
        "smallest p across (panel, cell line, treatment time, dose). Cell ",
        "colour is signed by Direction (blue = matches disease, red = ",
        "opposes disease). Ring = nominal-sig (P<0.05); black square = ",
        "FDR-sig; blank = not tested in that GWAS."),
      hr(),
      tags$details(class = "gd-details",
        tags$summary("Filter data"),
        tags$div(class = "gd-details-body",
          fluidRow(
            column(3,
              numericInput(ns("pert_sig_threshold"), "Significance threshold:",
                            value = 0.05, min = 1e-12, max = 1, step = 0.01)
            ),
            column(3,
              checkboxInput(ns("pert_only_recurrent"),
                            "Only show perturbations significant in ≥ k GWAS",
                            value = TRUE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("pert_only_recurrent")),
                uiOutput(ns("pert_k_slider_ui"))
              ),
              numericInput(ns("pert_row_cap"), "Rows shown in heatmap:",
                            value = 50, min = 5, max = 500, step = 5)
            ),
            column(3,
              selectInput(ns("pert_gwas_sort"), "Order GWAS by:",
                           choices = .compare_gwas_sort_choices,
                           selected = "as_selected"),
              sliderInput(ns("pert_plot_font_size"), "Font size (pt):",
                           min = 8, max = 20, value = 11, step = 1),
              sliderInput(ns("pert_plot_point_size"), "Point size:",
                           min = 2, max = 10, value = 4, step = 1)
            ),
            .dl_and_download_column(ns, "pert", default_h = 12)
          )
        )
      ),
      br(),
      tags$div(style = "max-width: 1100px; overflow-x: auto;",
        uiOutput(ns("cmap_pert_plot_ui"))
      ),
      gd_legend(list(
        "Ring around a circle"         = "Nominal-significant (P < 0.05).",
        "Black square around a circle" = "FDR-significant (P.FDR < 0.05).",
        "Blank cell"                   = "Perturbation was not tested in that GWAS.",
        "Blue"                         = "Direction 'Matches disease' — perturbation shares the trait's TWAS signature.",
        "Red"                          = "Direction 'Opposes disease' — perturbation counteracts the trait's TWAS signature."
      ), heading = "How to read this heatmap"),
      br(),
      tags$div(style = "max-width: 1100px;",
        h4("Underlying data"),
        DT::DTOutput(ns("cmap_pert_tbl"))
      )
    ),
    tabPanel("MOA", br(),
      p("Cross-GWAS CMap mechanism-of-action compare. Rows are MOA classes. ",
        "Best-per-cell reduction picks the smallest p across (panel, cell ",
        "line). Same colouring / overlay conventions as the Perturbation ",
        "view."),
      hr(),
      tags$details(class = "gd-details",
        tags$summary("Filter data"),
        tags$div(class = "gd-details-body",
          fluidRow(
            column(3,
              numericInput(ns("moa_sig_threshold"), "Significance threshold:",
                            value = 0.05, min = 1e-12, max = 1, step = 0.01)
            ),
            column(3,
              checkboxInput(ns("moa_only_recurrent"),
                            "Only show MOAs significant in ≥ k GWAS",
                            value = TRUE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("moa_only_recurrent")),
                uiOutput(ns("moa_k_slider_ui"))
              ),
              numericInput(ns("moa_row_cap"), "Rows shown in heatmap:",
                            value = 50, min = 5, max = 500, step = 5)
            ),
            column(3,
              selectInput(ns("moa_gwas_sort"), "Order GWAS by:",
                           choices = .compare_gwas_sort_choices,
                           selected = "as_selected"),
              sliderInput(ns("moa_plot_font_size"), "Font size (pt):",
                           min = 8, max = 20, value = 11, step = 1),
              sliderInput(ns("moa_plot_point_size"), "Point size:",
                           min = 2, max = 10, value = 4, step = 1)
            ),
            .dl_and_download_column(ns, "moa", default_h = 12)
          )
        )
      ),
      br(),
      tags$div(style = "max-width: 1100px; overflow-x: auto;",
        uiOutput(ns("cmap_moa_plot_ui"))
      ),
      gd_legend(list(
        "Ring around a circle"         = "Nominal-significant (P < 0.05).",
        "Black square around a circle" = "FDR-significant (P.FDR < 0.05).",
        "Blank cell"                   = "MOA was not tested in that GWAS.",
        "Blue"                         = "Direction 'Matches disease' — MOA shares the trait's TWAS signature.",
        "Red"                          = "Direction 'Opposes disease' — MOA counteracts the trait's TWAS signature."
      ), heading = "How to read this heatmap"),
      br(),
      tags$div(style = "max-width: 1100px;",
        h4("Underlying data"),
        DT::DTOutput(ns("cmap_moa_tbl"))
      )
    )
  )
}

# Shared point-style signed heatmap for CMap. `method_pick` selects
# "CMAP-perturbation" or "CMAP-MOA".
.cmap_signed_ggplot <- function(long, gwas_vec, method_pick,
                                    only_recurrent, k_min,
                                    sig_basis = "fdr", sig_threshold = 0.05,
                                    row_cap = 50, font_size = 11,
                                    point_size = 4) {
  slice <- long[method == method_pick & entity_type == "cmap" &
                    gwas %in% gwas_vec]
  if (nrow(slice) == 0) return(NULL)
  slice <- pick_best_per_cell(slice, c("gwas", "entity_id"))

  rec <- .atc_row_order(slice, sig_basis, sig_threshold)
  if (isTRUE(only_recurrent)) rec <- rec[n_sig >= k_min]
  if (nrow(rec) == 0) return(NULL)
  cap <- max(1L, min(as.integer(row_cap), nrow(rec)))
  rec <- rec[seq_len(cap)]
  slice <- slice[entity_id %in% rec$entity_id]

  slice[, entity_id := factor(entity_id, levels = rev(rec$entity_id))]
  slice[, gwas := factor(gwas, levels = gwas_vec)]
  slice[, nom_sig := !is.na(p)   & p   < 0.05]
  slice[, fdr_sig := !is.na(fdr) & fdr < 0.05]
  neg_log_fdr <- -log10(slice$fdr)
  slice[, signed_score := ifelse(
    is.na(direction), NA_real_,
    ifelse(direction == "Matches disease",  neg_log_fdr,
    ifelse(direction == "Opposes disease", -neg_log_fdr, NA_real_)))]
  # Drop rows with no signed score so blanks are truly blank (matches the
  # caption "Blank: not tested in that GWAS.").
  slice <- slice[!is.na(signed_score)]
  if (nrow(slice) == 0) return(NULL)

  # Data-driven symmetric limits so the ramp matches the observed range.
  max_abs <- suppressWarnings(max(abs(slice$signed_score), na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  fill_lim <- c(-max_abs, max_abs)

  base_layer <- ggplot2::geom_point(
    ggplot2::aes(fill = signed_score), shape = 21, stroke = 0,
    size = point_size)

  gg <- ggplot2::ggplot(slice, ggplot2::aes(x = gwas, y = entity_id)) +
    base_layer +
    ggplot2::scale_fill_gradient2(
      low = .gd_red, mid = "white", high = .gd_blue, midpoint = 0,
      limits = fill_lim, na.value = "transparent",
      name = "signed -log10(FDR)  (+ matches / − opposes)"
    )

  gg <- .compare_add_sig_overlays(gg, slice, base_layer,
                                    nominal_flag = "nom_sig",
                                    fdr_flag = "fdr_sig",
                                    point_size = point_size)

  gg +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid.major = ggplot2::element_line(colour = "#eef1f6"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(t = 60, r = 20, b = 10, l = 10, unit = "pt")
    )
}

.cmap_pert_ggplot <- function(...) .cmap_signed_ggplot(..., method_pick = "CMAP-perturbation")
.cmap_moa_ggplot  <- function(...) .cmap_signed_ggplot(..., method_pick = "CMAP-MOA")

cmap_compare_server <- function(id, gwas_data, selected_gwas_multi,
                                    comparison_long) {
  moduleServer(id, function(input, output, session) {

    output$pert_k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$pert_k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("pert_k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
    })
    output$moa_k_slider_ui <- renderUI({
      n_sel <- length(selected_gwas_multi())
      cur <- isolate(input$moa_k_min)
      cur <- if (is.null(cur)) 2L else as.integer(cur)
      sliderInput(session$ns("moa_k_min"), "k:",
                   min = 1, max = max(n_sel, 1L),
                   value = min(cur, max(n_sel, 1L)), step = 1)
    })

    pert_gwas_vec <- reactive({
      order_gwas(selected_gwas_multi(),
                  if (is.null(input$pert_gwas_sort)) "as_selected" else input$pert_gwas_sort,
                  NULL)
    })
    moa_gwas_vec <- reactive({
      order_gwas(selected_gwas_multi(),
                  if (is.null(input$moa_gwas_sort)) "as_selected" else input$moa_gwas_sort,
                  NULL)
    })

    # ---------- Perturbation plot / table ----------

    pert_plot <- reactive({
      req(comparison_long())
      .cmap_pert_ggplot(
        long           = comparison_long(),
        gwas_vec       = pert_gwas_vec(),
        only_recurrent = isTRUE(input$pert_only_recurrent),
        k_min          = if (is.null(input$pert_k_min)) 2L else as.integer(input$pert_k_min),
        sig_basis      = "fdr",
        sig_threshold  = if (is.null(input$pert_sig_threshold)) 0.05 else as.numeric(input$pert_sig_threshold),
        row_cap        = if (is.null(input$pert_row_cap)) 50L else as.integer(input$pert_row_cap),
        font_size      = if (is.null(input$pert_plot_font_size)) 11 else as.numeric(input$pert_plot_font_size),
        point_size     = if (is.null(input$pert_plot_point_size)) 4 else as.numeric(input$pert_plot_point_size)
      )
    })

    pert_dims <- reactive({
      p <- pert_plot()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_id))
      .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = length(pert_gwas_vec()),
        font_size = if (is.null(input$pert_plot_font_size)) 11 else as.numeric(input$pert_plot_font_size),
        y_labels  = if (is.null(p)) NULL else as.character(unique(p$data$entity_id))
      )
    })

    .sync_dl_dims(session, pert_dims, "pert")

    output$cmap_pert_plot_ui <- renderUI({
      if (is.null(pert_plot())) return(.compare_empty_state("perturbations"))
      dims <- pert_dims()
      plotOutput(session$ns("cmap_pert_plot"),
                  height = dims$height, width = dims$width)
    })

    output$cmap_pert_plot <- renderPlot({
      p <- pert_plot()
      if (is.null(p)) { plot.new(); title("No perturbations meet the current filter"); return(invisible()) }
      print(p)
    })

    pert_tbl_slice <- reactive({
      req(comparison_long())
      comparison_long()[method == "CMAP-perturbation" & entity_type == "cmap" &
                          gwas %in% pert_gwas_vec()]
    })

    output$cmap_pert_tbl <- DT::renderDT({
      slice <- pert_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS       = factor(slice$gwas, levels = pert_gwas_vec()),
        Perturbation = slice$entity_id,
        Panel      = slice$panel,
        Estimate   = .fmt_sig(slice$statistic, 3),
        SE         = .fmt_sig(slice$se, 3),
        Direction  = slice$direction,
        Reversal_Z = .fmt_sig(slice$reversal_z, 3),
        P          = .fmt_sig(slice$p, 3),
        `P.FDR`    = .fmt_sig(slice$fdr, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(8, "asc"))))
    }, server = TRUE)

    observeEvent(input$cmap_pert_tbl_rows_selected, {
      sel <- input$cmap_pert_tbl_rows_selected
      if (length(sel) != 1) return()
      slice <- pert_tbl_slice()
      if (nrow(slice) < sel) return()
      .show_entity_detail(as.character(slice$entity_id[sel]),
                            "cmap", comparison_long(), pert_gwas_vec())
    }, ignoreInit = TRUE)

    output$pert_download_plot <- downloadHandler(
      filename = function() sprintf("cmap_pert_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$pert_dl_format),
      content = function(file) {
        p <- pert_plot(); if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$pert_dl_format
        w <- input$pert_dl_width; h <- input$pert_dl_height
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$pert_download_csv <- downloadHandler(
      filename = function() sprintf("cmap_pert_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- pert_gwas_vec()
        long <- comparison_long()[method == "CMAP-perturbation" &
                                    entity_type == "cmap" &
                                    gwas %in% gwas_vec]
        best <- pick_best_per_cell(long, c("gwas", "entity_id"))
        best[, signed := ifelse(!is.na(direction) & direction == "Opposes disease",
                                 -(-log10(fdr)), -log10(fdr))]
        wide <- pivot_matrix(best, "signed", gwas_vec)
        display <- ifelse(is.na(wide), "NT", sprintf("%+.2f", wide))
        display <- matrix(display, nrow = nrow(wide), dimnames = dimnames(wide))
        header <- sprintf("# GenoDisc CMap perturbation compare CSV | %s | sig_basis=%s threshold=%g k_min=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           "fdr" %||% "fdr",
                           input$pert_sig_threshold %||% 0.05,
                           input$pert_k_min %||% 2L)
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: signed -log10(FDR): + = matches, - = opposes. NT = not tested.", con)
        out <- data.frame(Perturbation = rownames(display), display,
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )

    # ---------- MOA plot / table ----------

    moa_plot <- reactive({
      req(comparison_long())
      .cmap_moa_ggplot(
        long           = comparison_long(),
        gwas_vec       = moa_gwas_vec(),
        only_recurrent = isTRUE(input$moa_only_recurrent),
        k_min          = if (is.null(input$moa_k_min)) 2L else as.integer(input$moa_k_min),
        sig_basis      = "fdr",
        sig_threshold  = if (is.null(input$moa_sig_threshold)) 0.05 else as.numeric(input$moa_sig_threshold),
        row_cap        = if (is.null(input$moa_row_cap)) 50L else as.integer(input$moa_row_cap),
        font_size      = if (is.null(input$moa_plot_font_size)) 11 else as.numeric(input$moa_plot_font_size),
        point_size     = if (is.null(input$moa_plot_point_size)) 4 else as.numeric(input$moa_plot_point_size)
      )
    })

    moa_dims <- reactive({
      p <- moa_plot()
      n_rows <- if (is.null(p)) 1L else length(unique(p$data$entity_id))
      .compare_plot_dims(
        n_rows    = n_rows,
        n_cols    = length(moa_gwas_vec()),
        font_size = if (is.null(input$moa_plot_font_size)) 11 else as.numeric(input$moa_plot_font_size),
        y_labels  = if (is.null(p)) NULL else as.character(unique(p$data$entity_id))
      )
    })

    .sync_dl_dims(session, moa_dims, "moa")

    output$cmap_moa_plot_ui <- renderUI({
      if (is.null(moa_plot())) return(.compare_empty_state("MOAs"))
      dims <- moa_dims()
      plotOutput(session$ns("cmap_moa_plot"),
                  height = dims$height, width = dims$width)
    })

    output$cmap_moa_plot <- renderPlot({
      p <- moa_plot()
      if (is.null(p)) { plot.new(); title("No MOAs meet the current filter"); return(invisible()) }
      print(p)
    })

    moa_tbl_slice <- reactive({
      req(comparison_long())
      comparison_long()[method == "CMAP-MOA" & entity_type == "cmap" &
                          gwas %in% moa_gwas_vec()]
    })

    output$cmap_moa_tbl <- DT::renderDT({
      slice <- moa_tbl_slice()
      if (nrow(slice) == 0) return(NULL)
      out <- data.frame(
        GWAS       = factor(slice$gwas, levels = moa_gwas_vec()),
        MOA        = slice$entity_id,
        Panel      = slice$panel,
        `N Drugs`  = .fmt_int(slice$n_units),
        Estimate   = .fmt_sig(slice$statistic, 3),
        Direction  = slice$direction,
        Reversal_Z = .fmt_sig(slice$reversal_z, 3),
        P          = .fmt_sig(slice$p, 3),
        `P.FDR`    = .fmt_sig(slice$fdr, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     selection = "none",
                     options = list(pageLength = 20, server = TRUE,
                                     order = list(list(8, "asc"))))
    }, server = TRUE)

    observeEvent(input$cmap_moa_tbl_rows_selected, {
      sel <- input$cmap_moa_tbl_rows_selected
      if (length(sel) != 1) return()
      slice <- moa_tbl_slice()
      if (nrow(slice) < sel) return()
      .show_entity_detail(as.character(slice$entity_id[sel]),
                            "cmap", comparison_long(), moa_gwas_vec())
    }, ignoreInit = TRUE)

    output$moa_download_plot <- downloadHandler(
      filename = function() sprintf("cmap_moa_compare_%s.%s",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     input$moa_dl_format),
      content = function(file) {
        p <- moa_plot(); if (is.null(p)) { grDevices::png(file); dev.off(); return() }
        fmt <- input$moa_dl_format
        w <- input$moa_dl_width; h <- input$moa_dl_height
        if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300)
        else if (fmt == "pdf") grDevices::pdf(file, width = w, height = h)
        else grDevices::svg(file, width = w, height = h)
        print(p); grDevices::dev.off()
      }
    )

    output$moa_download_csv <- downloadHandler(
      filename = function() sprintf("cmap_moa_compare_matrix_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        gwas_vec <- moa_gwas_vec()
        long <- comparison_long()[method == "CMAP-MOA" &
                                    entity_type == "cmap" &
                                    gwas %in% gwas_vec]
        best <- pick_best_per_cell(long, c("gwas", "entity_id"))
        best[, signed := ifelse(!is.na(direction) & direction == "Opposes disease",
                                 -(-log10(fdr)), -log10(fdr))]
        wide <- pivot_matrix(best, "signed", gwas_vec)
        display <- ifelse(is.na(wide), "NT", sprintf("%+.2f", wide))
        display <- matrix(display, nrow = nrow(wide), dimnames = dimnames(wide))
        header <- sprintf("# GenoDisc CMap MOA compare CSV | %s | sig_basis=%s threshold=%g k_min=%s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           "fdr" %||% "fdr",
                           input$moa_sig_threshold %||% 0.05,
                           input$moa_k_min %||% 2L)
        con <- file(file, "w"); on.exit(close(con))
        writeLines(header, con)
        writeLines("# Cells: signed -log10(FDR): + = matches, - = opposes. NT = not tested.", con)
        out <- data.frame(MOA = rownames(display), display,
                          check.names = FALSE, stringsAsFactors = FALSE)
        utils::write.csv(out, con, row.names = FALSE, quote = FALSE)
      }
    )
  })
}
