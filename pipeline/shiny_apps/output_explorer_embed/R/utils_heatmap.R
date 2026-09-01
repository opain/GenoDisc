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

#' Fill missing Panel+Method combinations with NA rows
#'
#' For each Panel+Method combo, find IDs from other combos that are missing
#' and insert placeholder rows with NA values.
#'
#' @param df Data frame with columns Panel, Method, and the id column, plus Z, Sig, Coloc, Type, Group
#' @param id_col Character name of the ID column (default "ID")
#' @return Data frame with missing rows inserted
fill_missing_panel_method <- function(df, id_col = "ID") {
  result <- NULL
  for (i in unique(df$Panel)) {
    for (j in unique(df$Method[df$Panel == i])) {

      df_panel <- df[df$Panel == i & df$Method == j, ]
      df_other <- df[!(df$Panel %in% df_panel$Panel) & !(df$Method %in% df_panel$Method), ]
      df_other <- df_other[!(df_other[[id_col]] %in% df_panel[[id_col]]), ]
      df_other <- unique(df_other[[id_col]])

      if (length(df_other) > 0) {
        df_panel_missing <- data.frame(x = df_other)
        names(df_panel_missing) <- id_col

        df_panel_missing$Panel <- i
        df_panel_missing[[id_col]] <- df_other
        df_panel_missing$Z <- NA
        df_panel_missing$Sig <- 0
        df_panel_missing$Coloc <- 0
        df_panel_missing$Method <- j
        df_panel_missing$Type <- df_panel$Type[1]
        df_panel_missing$Group <- df_panel$Group[1]

        df_panel_missing <- df_panel_missing[, names(df_panel)]

        result <- rbind(result, df_panel_missing)
      }

      result <- rbind(result, df_panel)
    }
  }
  result
}

#' Insert placeholder rows for each unique Panel+Method combo
#'
#' Creates one placeholder row per unique Panel+Method combination with
#' ID='Placeholder', Z=NA, Sig=NA, Coloc=NA.
#'
#' @param df Data frame with columns Panel, Method, and id_col, Z, Sig, Coloc
#' @param id_col Character name of the ID column (default "ID")
#' @return Data frame with placeholder rows prepended
insert_placeholders <- function(df, id_col = "ID") {
  na_rows <- df[!(duplicated(paste0(df$Panel, df$Method))), ]
  na_rows[[id_col]] <- 'Placeholder'
  na_rows$Z <- NA
  na_rows$Sig <- NA
  na_rows$Coloc <- NA
  rbind(na_rows, df)
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
#' - Left pad: only added if the leftmost facet's leftmost 45° tick label
#'   would extend past the y-axis area
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
                           facet_order = NULL, font_size = 14, min_height = 0) {
  if (nrow(df) == 0) {
    return(list(height = 100, width = 100, left_pad_pt = 0, panel_h_pt = 0))
  }

  # Height: x-label rotation space + rows + padding, scaled by font size.
  fs_ratio <- font_size / 14
  num_row <- length(unique(df[[y_col]]))
  # Natural panel-only height (what ggplot would give the panel row with
  # no min_height floor). Used by `reserve_legend_space` to pin the panel
  # row when we bump plot_height up to fit the legend.
  panel_h_pt <- num_row * 20 * fs_ratio
  plot_height <- (max(nchar(df[[x_col]]), na.rm = TRUE) * 3 * fs_ratio) +
                 panel_h_pt + 100
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

  # Left pad: only needed if the leftmost facet's leftmost 45° tick label
  # extends past the y-axis area. Otherwise the y-axis absorbs the overflow.
  left_pad_pt <- 0
  if (length(facets) > 0) {
    leftmost <- facets[1]
    leftmost_panels <- unique(df[[x_col]][df[[facet_col]] == leftmost])
    n_x_left <- length(leftmost_panels)
    if (n_x_left > 0) {
      label_slot <- max(strwidth_pt(leftmost_panels, ps = font_size)) * cos(pi / 4)
      leftmost_panel_width <- facet_target_width_pt(leftmost, leftmost_panels, font_size)
      leftmost_tick_from_edge <- leftmost_panel_width / (2 * n_x_left)
      overflow_past_panel <- max(0, label_slot - leftmost_tick_from_edge)
      # Approximate rendered y-axis area: strwidth of longest y-label at the
      # plot font size + tick length + text-to-tick padding (~15pt combined).
      y_axis_actual_pt <- max(strwidth_pt(unique(df[[y_col]]), ps = font_size)) + 15
      left_pad_pt <- max(0, overflow_past_panel - y_axis_actual_pt)
    }
  }

  # Right-side allowance for Z-score legend + panel spacing + margins.
  plot_width <- y_label_width + total_facet_width + 150 + left_pad_pt

  list(height = plot_height, width = plot_width, left_pad_pt = left_pad_pt,
       panel_h_pt = panel_h_pt)
}

#' Calculate group sizes and proportional widths
#'
#' @param df Data frame with Panel column and the group column
#' @param group_col Character name of the group/method column (default "Group")
#' @param groups Character vector of group levels in order
#' @return Data frame with columns Group, Size, Prop, Width
calc_group_widths <- function(df, group_col = "Group", groups = NULL) {
  if (is.null(groups)) {
    groups <- unique(df[[group_col]])
  }
  group_siz <- NULL
  for (i in groups) {
    group_siz <- rbind(group_siz, data.frame(Group = i,
                                              Size = length(unique(df$Panel[df[[group_col]] == i]))))
  }
  # Set minimum size to 2 to allow space for labels
  group_siz$Size[group_siz$Size < 2] <- 2
  group_siz$Prop <- group_siz$Size / sum(group_siz$Size)
  group_siz$Width <- 4 * group_siz$Prop
  group_siz
}

#' Create the Z-score colour scale
#'
#' @param z_values Numeric vector of Z-score values
#' @return A ggplot2 scale_colour_gradientn object
z_score_colour_scale <- function(z_values) {
  x <- c(-max(abs(z_values), na.rm = TRUE), 0, max(abs(z_values), na.rm = TRUE))
  x <- (x - min(x)) / (max(x) - min(x))
  ggplot2::scale_colour_gradientn(
    colours = c("#0066FF", "#0099FF", "#FFFFFF", "#FF6666", "#FF0000"),
    na.value = NA,
    name = "Z-score",
    limits = c(-max(abs(z_values), na.rm = TRUE), max(abs(z_values), na.rm = TRUE)),
    values = x
  )
}

#' Reserve legend space by pinning the panel row height and letting the
#' legend column stretch below it.
#'
#' When a heatmap has only a handful of rows the panel is very short, so the
#' right-hand colourbar legend gets clipped at the top. This helper:
#'   1. fixes the panel row(s) to `panel_h_pt` (their natural pt height so
#'      row spacing doesn't change), then
#'   2. adds a `1null` filler row at the bottom to absorb any extra device
#'      height, then
#'   3. extends the guide-box's bottom row index to that new filler so the
#'      legend has vertical room to draw at full size.
#' No-op if there is no guide-box (`grep("guide", gt$layout$name)` is empty)
#' or if `panel_h_pt` is `NA`/`NULL`.
reserve_legend_space <- function(gt, panel_h_pt) {
  if (is.null(panel_h_pt) || is.na(panel_h_pt) || panel_h_pt <= 0) return(gt)
  panel_t <- unique(gt$layout$t[grepl("^panel", gt$layout$name)])
  if (length(panel_t) == 0) return(gt)
  gt$heights[panel_t] <- grid::unit(panel_h_pt, "pt")
  gt <- gtable::gtable_add_rows(gt, grid::unit(1, "null"), pos = -1)
  gb_idx <- grep("guide", gt$layout$name)
  if (length(gb_idx) > 0) gt$layout$b[gb_idx] <- nrow(gt)
  gt
}

#' Manipulate ggplot gtable to adjust facet panel widths and draw
#'
#' @param plot A ggplot object
#' @param group_widths Data frame with Width column, one row per facet panel
draw_faceted_heatmap <- function(plot, group_widths) {
  gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
  for (i in 1:nrow(group_widths)) {
    gt$widths[gt$layout$l[grep(paste0('panel-', i, '-1'), gt$layout$name)]] <-
      group_widths$Width[i] * gt$widths[gt$layout$l[grep(paste0('panel-', i, '-1'), gt$layout$name)]]
  }
  grid::grid.draw(gt)
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
