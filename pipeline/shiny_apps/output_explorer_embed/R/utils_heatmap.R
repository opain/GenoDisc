#' Create group labels from method and type
#'
#' Pastes method and type with newline, then cleans up special cases:
#' 'SNP\nFine-mapping\n' -> 'SuSiE', 'MAGMA\n' -> 'MAGMA', 'Nearest\nGene\n' -> 'Nearest\nGene'
#'
#' @param method Character vector of method names
#' @param type Character vector of type names
#' @return Character vector of cleaned group labels
make_group_labels <- function(method, type) {
  res_group <- paste0(method, '\n', type)
  res_group[res_group == 'SNP\nFine-mapping\n'] <- 'SuSiE'
  res_group[res_group == 'MAGMA\n'] <- 'MAGMA'
  res_group[res_group == 'Nearest\nGene\n'] <- 'Nearest\nGene'
  res_group
}

#' Measure rendered width of strings at the given point size, returned in pt.
#' Uses a null pdf device so no side effects on the caller.
strwidth_pt <- function(s, ps = 14) {
  if (length(s) == 0) return(numeric(0))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  graphics::par(ps = ps)
  graphics::strwidth(as.character(s), units = "inches") * 72
}

#' Target width in pt for a single facet group
#'
#' Shared by `calc_plot_dims` (width budget) and the render code (absolute
#' panel widths) so the two stay in lockstep. Floor is the strip title
#' width (widest line, measured with `strwidth` at the plot's font size, so
#' longer labels at larger fonts get proportionally wider facets). Panel
#' names are only used for the arity `n_x`; overflow from long rotated tick
#' labels is absorbed by extra `plot.margin` on the left of the plot, not by
#' widening panels.
facet_target_width_pt <- function(group_label, panel_names, font_size = 14) {
  n_x <- length(panel_names)
  lines <- unlist(strsplit(group_label, "\n"))
  title_pt <- max(strwidth_pt(lines, ps = font_size)) + 10
  # Per-panel budget also scales with font so points don't cramp at big fonts.
  per_panel <- 25 * (font_size / 14)
  max(n_x * per_panel, title_pt)
}

#' Calculate plot dimensions from data
#'
#' Computes height and width for faceted heatmap plots based on:
#' - Height: number of y-axis items + x-axis label rotation space
#' - Width: y-axis label space + per-facet width + right-side legend allowance
#' - Left pad: horizontal projection of the leftmost facet's widest 45°
#'   tick label so the label sits fully in the plot's left margin (no
#'   overlap with the y-axis label area)
#'
#' @param df Data frame containing the plot data
#' @param y_col Column used for y-axis labels (default "ID")
#' @param x_col Column used for x-axis categories (default "Panel")
#' @param facet_col Column used for facet_wrap (default "Method")
#' @param facet_order Optional character vector giving the left-to-right facet
#'   order in the rendered plot. Used to identify the leftmost facet for the
#'   overflow calculation; if `NULL`, uses `unique(df[[facet_col]])`.
#' @param font_size Base text size in pt. Facet titles, y-axis text and
#'   rotated x-tick labels all scale with this, so per-facet widths and the
#'   overall plot width grow proportionally when the user picks a larger font.
#' @return List with height (pt), width (pt), and left_pad_pt (pt) for the
#'   `plot.margin(l = ...)` needed to absorb tick-label overflow.
calc_plot_dims <- function(df, y_col = "ID", x_col = "Panel", facet_col = "Method",
                           facet_order = NULL, font_size = 14, min_height = 0,
                           min_panel_h_pt = 0) {
  if (nrow(df) == 0) {
    return(list(height = 100, width = 100, left_pad_pt = 0, panel_h_pt = 0))
  }

  # Height: x-label rotation space + rows + padding, scaled by font size.
  fs_ratio <- font_size / 14
  num_row <- length(unique(df[[y_col]]))
  # Natural panel-only height (num_row * 20 * fs_ratio). We floor it at
  # `min_panel_h_pt` so short plots have a panel row tall enough to fit
  # the right-hand legend without clipping — that means row spacing
  # spreads a bit for very-few-row plots, in exchange for the legend
  # staying visible and vertically centred on the panel.
  panel_h_pt <- max(num_row * 20 * fs_ratio, min_panel_h_pt)
  # +30pt covers the fixed 15pt top + 15pt bottom padding rows that
  # `reserve_legend_space` inserts, so the device has room for them
  # without squeezing the panel.
  plot_height <- (max(nchar(df[[x_col]]), na.rm = TRUE) * 3 * fs_ratio) +
                 panel_h_pt + 100 + 30
  plot_height <- max(plot_height, min_height)

  # Facet order (left→right)
  present <- unique(df[[facet_col]])
  facets <- if (!is.null(facet_order)) {
    c(intersect(facet_order, present), setdiff(present, facet_order))
  } else present

  # Width: y-label margin + sum of per-facet widths
  y_label_width <- max(strwidth_pt(unique(df[[y_col]]), ps = font_size)) + 10
  total_facet_width <- 0
  for (f in facets) {
    panel_names <- unique(df[[x_col]][df[[facet_col]] == f])
    total_facet_width <- total_facet_width + facet_target_width_pt(f, panel_names, font_size)
  }

  # Left pad = horizontal projection of the leftmost facet's widest 45°
  # rotated tick label, minus a small tick-padding buffer (the label may
  # extend into the ~15pt gap between panel and y-axis without issue).
  # Reserving the (near-)full label slot means the label sits mostly in
  # the margin without overlapping the y-axis label area.
  left_pad_pt <- 0
  if (length(facets) > 0) {
    leftmost_panels <- unique(df[[x_col]][df[[facet_col]] == facets[1]])
    if (length(leftmost_panels) > 0) {
      label_slot <- max(strwidth_pt(leftmost_panels, ps = font_size)) * cos(pi / 4)
      left_pad_pt <- max(0, label_slot - 15)
    }
  }

  # Right-side allowance for Z-score legend + panel spacing + margins.
  plot_width <- y_label_width + total_facet_width + 150 + left_pad_pt

  list(height = plot_height, width = plot_width, left_pad_pt = left_pad_pt,
       panel_h_pt = panel_h_pt)
}

#' Reserve legend space by pinning the panel row height and letting the
#' legend column stretch below it.
#'
#' When a heatmap has only a handful of rows the panel is very short, so the
#' right-hand colourbar legend gets clipped at the top. This helper:
#'   1. fixes the panel row to `panel_h_pt` (its natural pt height so row
#'      spacing doesn't change), then
#'   2. adds a `1null` filler row at the bottom to absorb any extra device
#'      height, then
#'   3. extends the guide-box's bottom row index to that new filler so the
#'      legend has vertical room to draw at full size.
#' No-op if there is no guide-box, if `panel_h_pt` is `NA`/`NULL`, or if
#' the gtable has more than one panel row (`facet_grid` layouts). For
#' multi-row facets `panel_h_pt` is the total across all rows, not per
#' row, so pinning every row to it would blow the layout up — the legend
#' also isn't at risk of clipping there because vertical space is already
#' abundant.
reserve_legend_space <- function(gt, panel_h_pt) {
  if (is.null(panel_h_pt) || is.na(panel_h_pt) || panel_h_pt <= 0) return(gt)
  panel_t <- unique(gt$layout$t[grepl("^panel", gt$layout$name)])
  if (length(panel_t) != 1) return(gt)
  gt$heights[panel_t] <- grid::unit(panel_h_pt, "pt")
  # Guaranteed 15pt padding at absolute top and bottom (breathing room
  # above the facet strip and below the axis text), plus a `1null` row
  # at each end so any additional device slack from `min_height` splits
  # evenly top/bottom instead of piling under the plot.
  gt <- gtable::gtable_add_rows(gt, grid::unit(15, "pt"), pos = 0)
  gt <- gtable::gtable_add_rows(gt, grid::unit(15, "pt"), pos = -1)
  gt <- gtable::gtable_add_rows(gt, grid::unit(1, "null"), pos = 0)
  gt <- gtable::gtable_add_rows(gt, grid::unit(1, "null"), pos = -1)
  gt
}

#' Get the standard group ordering vector for molecular association heatmaps
#'
#' @return Character vector of group names in standard order
get_group_order <- function() {
  c('SuSiE', 'FUSION\nExpr.', 'FUSION\nSplice', 'SMR\nExpr.',
    'FUSION\nProtein', 'SMR\nProtein', 'MAGMA', 'Nearest\nGene')
}

#' Group features into genomic loci by a fixed base-pair window
#'
#' Sorts features by chromosome then position and starts a new locus at each
#' chromosome change or whenever the gap from the previous feature exceeds
#' `window`. Each locus gets a coordinate-range label (e.g. "chr3:12.1-12.4Mb")
#' and an integer order that sorts loci by genomic position.
#'
#' @param pos_df data.frame with columns ID, CHR, BP (one row per feature)
#' @param window Base-pair gap that starts a new locus (default 5e5 = 500 kb)
#' @return data.frame with columns ID, CHR, BP, Locus, locus_order (empty if no
#'   positioned features are supplied)
assign_loci <- function(pos_df, window = 5e5) {
  empty <- data.frame(ID = character(0), CHR = numeric(0), BP = numeric(0),
                      Locus = character(0), locus_order = integer(0),
                      stringsAsFactors = FALSE)
  if (is.null(pos_df) || nrow(pos_df) == 0) return(empty)

  d <- pos_df[!is.na(pos_df$CHR) & !is.na(pos_df$BP), c("ID", "CHR", "BP")]
  if (nrow(d) == 0) return(empty)
  d <- d[order(d$CHR, d$BP), ]

  n <- nrow(d)
  brk <- rep(TRUE, n)
  if (n > 1) {
    brk[2:n] <- (d$CHR[2:n] != d$CHR[1:(n - 1)]) |
                ((d$BP[2:n] - d$BP[1:(n - 1)]) > window)
  }
  d$locus_order <- cumsum(brk)

  loc_min <- tapply(d$BP, d$locus_order, min)
  loc_max <- tapply(d$BP, d$locus_order, max)
  loc_chr <- tapply(d$CHR, d$locus_order, function(x) x[1])
  lab <- sprintf("chr%s:%.1f–%.1fMb", loc_chr, loc_min / 1e6, loc_max / 1e6)
  names(lab) <- names(loc_min)
  d$Locus <- lab[as.character(d$locus_order)]
  d
}
