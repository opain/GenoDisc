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

#' Calculate plot dimensions from data
#'
#' @param df Data frame with Panel, Method, Type columns and the id column
#' @param id_col Character name of the ID column (default "ID")
#' @return List with height and width
calc_plot_dims <- function(df, id_col = "ID") {
  if (nrow(df) > 0) {
    num_row <- length(unique(df[[id_col]]))
    plot_height <- (max(nchar(df$Panel)) * 3) + (num_row * 20) + 100
    num_col <- length(unique(paste0(df$Panel, '_', df$Method, '_', df$Type)))
    num_pan <- length(unique(df$Method))
    plot_width <- 120 + (max(nchar(df[[id_col]]), na.rm = TRUE) * 4) + (num_col * 27) + (num_pan * 15)
    plot_width <- max(plot_width, (length(unique(df$Method)) * 140))
  } else {
    plot_height <- 100
    plot_width <- 100
  }
  list(height = plot_height, width = plot_width)
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
