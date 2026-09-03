#' Resolve the "Plot options" theme dropdown value to a ggplot theme function.
enr_theme_fn <- function(name) {
  switch(name %||% "bw",
    bw = ggplot2::theme_bw,
    minimal = ggplot2::theme_minimal,
    classic = ggplot2::theme_classic,
    light = ggplot2::theme_light,
    ggplot2::theme_bw
  )
}

#' Build the tissue-enrichment lollipop plot as a ggplot object.
#'
#' Shared by the on-screen `renderPlot` and the download handler so saved
#' files match what's on screen. Returns `NULL` if there's no data.
#'
#' @param d Filtered tissue data (Tissue, negLog10P, FDR_Sig, Retained).
#' @param n_total Total number of tissues before filtering — used to derive
#'   the Bonferroni threshold line.
#' @param font_size Base font size in pt.
#' @param point_size Point size for the dots.
#' @param theme_fn A ggplot theme function (e.g. `ggplot2::theme_bw`).
#' @param title Optional plot title; empty string draws no title.
build_tissue_plot <- function(d, n_total, sort_choice = "significance",
                               font_size = 13, point_size = 3,
                               theme_fn = ggplot2::theme_bw, title = "") {
  if (is.null(d) || nrow(d) == 0) return(NULL)

  # For y-axis, first factor level renders at the BOTTOM of the plot.
  # Significance: least-significant at bottom, most at top.
  # Alphabetical: A at top, Z at bottom → reverse-alphabetical factor order.
  if (identical(sort_choice, "alphabetical")) {
    d <- d[order(as.character(d$Tissue), decreasing = TRUE), ]
  } else {
    d <- d[order(d$negLog10P), ]
  }
  d$Label <- ifelse(d$Retained, paste0(d$Tissue, " *"), d$Tissue)
  d$Label <- factor(d$Label, levels = d$Label)
  face_vec <- ifelse(d$Retained, "bold", "plain")

  nom_line  <- -log10(0.05)
  bonf_line <- -log10(0.05 / max(n_total, 1))

  ggplot2::ggplot(d, ggplot2::aes(x = negLog10P, y = Label)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = negLog10P, yend = Label),
                          colour = "grey78", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = nom_line,  linetype = "dashed", colour = "grey55") +
    ggplot2::geom_vline(xintercept = bonf_line, linetype = "dotted", colour = "grey35") +
    ggplot2::geom_point(ggplot2::aes(fill = FDR_Sig),
                        shape = 21, size = point_size, colour = "black", stroke = 0.4) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "white", `TRUE` = "#0f766e"),
      labels = c(`FALSE` = "Not FDR-significant", `TRUE` = "FDR-significant"),
      name = NULL) +
    ggplot2::labs(x = expression(-log[10](italic(P))), y = NULL,
                  title = if (nzchar(title)) title else NULL) +
    theme_fn(base_size = font_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = face_vec),
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

#' Build the drug-enrichment summary heatmap as a gtable
#'
#' Shared by the on-screen `renderPlot` and the download handler so downloaded
#' files match what's on screen. Returns `NULL` if there's no data.
build_tx_drug_gtable <- function(all_gs, sort_choice = "Alphabetical",
                                  font_size = 14, point_size = 5,
                                  theme_fn = ggplot2::theme_bw, title = "",
                                  panel_h_pt = NULL, left_pad_pt = 0) {
  if (is.null(all_gs) || nrow(all_gs) == 0) return(NULL)
  all_gs$`ATC Code` <- NULL

  # Insert missing Panel × Method combinations so the heatmap grid is complete.
  all_gs_all <- NULL
  for (i in unique(all_gs$Panel)) {
    for (j in unique(all_gs$Method[all_gs$Panel == i])) {
      all_gs_panel <- all_gs[all_gs$Panel == i & all_gs$Method == j, ]
      all_gs_other <- all_gs[!(all_gs$Panel %in% all_gs_panel$Panel) & !(all_gs$Method %in% all_gs_panel$Method), ]
      all_gs_other <- all_gs_other[!(all_gs_other$Name %in% all_gs_panel$Name), ]
      all_gs_other <- unique(all_gs_other$Name)

      if (length(all_gs_other) > 0) {
        all_gs_panel_missing <- data.frame(Name = all_gs_other)
        all_gs_panel_missing$Panel  <- i
        all_gs_panel_missing$Name   <- all_gs_other
        all_gs_panel_missing$Z      <- NA
        all_gs_panel_missing$P      <- NA
        all_gs_panel_missing$P.FDR  <- NA
        all_gs_panel_missing$Method <- j
        all_gs_panel_missing <- all_gs_panel_missing[, names(all_gs_panel)]
        all_gs_all <- rbind(all_gs_all, all_gs_panel_missing)
      }
      all_gs_all <- rbind(all_gs_all, all_gs_panel)
    }
  }
  all_gs_all <- all_gs_all[all_gs_all$Name != 'Placeholder', ]
  # Filters may leave only Placeholder rows (e.g. no FDR-significant drugs in
  # the selected methods). Bail before ggplot's facet_wrap errors on an empty
  # factor; the "no drugs" UI message explains the situation to the user.
  if (nrow(all_gs_all) == 0) return(NULL)

  methods <- c('MAGMA','GCSC','TWAS-GSEA','TWAS-GSEA (non-dir)')
  methods <- methods[methods %in% all_gs_all$Method]
  all_gs_all$Method <- factor(all_gs_all$Method, levels = methods)

  # Sort levels by user choice
  if (sort_choice == 'Alphabetical') {
    all_gs_all$Name <- factor(all_gs_all$Name, levels = unique(all_gs_all$Name[rev(order(all_gs_all$Name))]))
  } else if (sort_choice == 'All - Z') {
    all_gs_all$Name <- factor(all_gs_all$Name, levels = rev(unique(rev(all_gs_all$Name[order(all_gs_all$Z, na.last = FALSE)]))))
  } else if (sort_choice == 'TWAS-GSEA - Z') {
    all_gs_all$Name <- factor(all_gs_all$Name, levels = rev(unique(rev(all_gs_all$Name[all_gs_all$Method == 'TWAS-GSEA'][order(all_gs_all$Z[all_gs_all$Method == 'TWAS-GSEA'], na.last = FALSE)]))))
  } else if (sort_choice == 'TWAS-GSEA (non-dir) - Z') {
    all_gs_all$Name <- factor(all_gs_all$Name, levels = rev(unique(rev(all_gs_all$Name[all_gs_all$Method == 'TWAS-GSEA (non-dir)'][order(all_gs_all$Z[all_gs_all$Method == 'TWAS-GSEA (non-dir)'], na.last = FALSE)]))))
  } else if (sort_choice == 'MAGMA - Z') {
    all_gs_all$Name <- factor(all_gs_all$Name, levels = unique(all_gs_all$Name[all_gs_all$Method == 'MAGMA'][order(all_gs_all$Z[all_gs_all$Method == 'MAGMA'], na.last = FALSE)]))
  } else if (sort_choice == 'GCSC - Z') {
    all_gs_all$Name <- factor(all_gs_all$Name, levels = unique(all_gs_all$Name[all_gs_all$Method == 'GCSC'][order(all_gs_all$Z[all_gs_all$Method == 'GCSC'], na.last = FALSE)]))
  }

  # Per-facet width scaling (min-2 floor, proportional to panel count).
  group_siz <- do.call(rbind, lapply(methods, function(m)
    data.frame(Group = m, Size = length(unique(all_gs_all$Panel[all_gs_all$Method == m])))))
  group_siz$Size[group_siz$Size < 2] <- 2
  group_siz$Prop <- group_siz$Size / sum(group_siz$Size)
  group_siz$Width <- 4 * group_siz$Prop

  all_gs_all <- data.table::data.table(all_gs_all)
  dir_data <- all_gs_all[Method == 'TWAS-GSEA']
  pos_data <- all_gs_all[Method != 'TWAS-GSEA']
  dir_max <- if (nrow(dir_data) > 0) suppressWarnings(max(abs(dir_data$Z), na.rm = TRUE)) else NA
  pos_max <- if (nrow(pos_data) > 0) suppressWarnings(max(pos_data$Z,      na.rm = TRUE)) else NA
  if (!is.finite(dir_max)) dir_max <- NA
  if (!is.finite(pos_max)) pos_max <- NA

  heatmap <- ggplot2::ggplot(data = all_gs_all, ggplot2::aes(x = Panel, y = Name)) +
    theme_fn()
  if (nrow(dir_data) > 0 && is.finite(dir_max)) {
    heatmap <- heatmap +
      ggplot2::geom_point(data = dir_data, ggplot2::aes(fill = Z), shape = 21, stroke = 0, size = point_size) +
      ggplot2::scale_fill_gradientn(colours = c("#0066FF","#0099FF","#FFFFFF","#FF6666","#FF0000"),
                                    na.value = NA, name = "TWAS-GSEA\nZ-score",
                                    limits = c(-dir_max, dir_max))
  }
  if (nrow(pos_data) > 0 && is.finite(pos_max)) {
    heatmap <- heatmap +
      ggplot2::geom_point(data = pos_data, ggplot2::aes(colour = Z), shape = 16, size = point_size) +
      ggplot2::scale_colour_gradientn(colours = c("#FFFFFF","#00CC66"),
                                      na.value = NA, name = "Enrichment\nZ-score",
                                      limits = c(0, pos_max))
  }
  heatmap <- heatmap +
    ggplot2::geom_point(data = all_gs_all[which(all_gs_all$P < 0.05), ], ggplot2::aes(x = Panel, y = Name), colour = 'black', fill = NA, size = point_size + 1) +
    ggplot2::geom_point(data = all_gs_all[which(all_gs_all$P.FDR < 0.05), ], ggplot2::aes(x = Panel, y = Name), colour = 'black', fill = NA, size = point_size + 2, shape = 15)
  if (nrow(dir_data) > 0) {
    heatmap <- heatmap + ggplot2::geom_point(data = dir_data, ggplot2::aes(fill = Z), shape = 21, stroke = 0, size = point_size)
  }
  if (nrow(pos_data) > 0) {
    heatmap <- heatmap + ggplot2::geom_point(data = pos_data, ggplot2::aes(colour = Z), shape = 16, size = point_size)
  }
  heatmap <- heatmap +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(x = '', y = '', title = if (nzchar(title)) title else NULL) +
    ggplot2::facet_wrap(~ Method, nrow = 1, scales = "free_x") +
    ggplot2::theme(text = ggplot2::element_text(size = font_size),
                   plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5,
                                                 l = left_pad_pt, unit = "pt"))

  gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(heatmap))
  for (i in seq_len(nrow(group_siz))) {
    panel_l <- gt$layout$l[grep(paste0('panel-', i, '-1'), gt$layout$name)]
    gt$widths[panel_l] <- group_siz$Width[i] * gt$widths[panel_l]
  }
  reserve_legend_space(gt, panel_h_pt)
}

#' Build the ATC-class enrichment summary heatmap as a gtable
build_tx_atc_gtable <- function(all_gs_atc, sort_choice = "Alphabetical",
                                 font_size = 14, point_size = 5,
                                 theme_fn = ggplot2::theme_bw, title = "",
                                 panel_h_pt = NULL, left_pad_pt = 0) {
  if (is.null(all_gs_atc) || nrow(all_gs_atc) == 0) return(NULL)

  all_gs_atc_all <- NULL
  for (i in unique(all_gs_atc$Panel)) {
    for (j in unique(all_gs_atc$Method[all_gs_atc$Panel == i])) {
      all_gs_atc_panel <- all_gs_atc[all_gs_atc$Panel == i & all_gs_atc$Method == j, ]
      all_gs_atc_other <- all_gs_atc[!(all_gs_atc$Panel %in% all_gs_atc_panel$Panel) & !(all_gs_atc$Method %in% all_gs_atc_panel$Method), ]
      all_gs_atc_other <- all_gs_atc_other[!(all_gs_atc_other$Name %in% all_gs_atc_panel$Name), ]
      all_gs_atc_other <- unique(all_gs_atc_other$Name)

      if (length(all_gs_atc_other) > 0) {
        all_gs_atc_panel_missing <- data.frame(Name = all_gs_atc_other)
        all_gs_atc_panel_missing$Panel   <- i
        all_gs_atc_panel_missing$Name    <- all_gs_atc_other
        all_gs_atc_panel_missing$Z       <- NA
        all_gs_atc_panel_missing$FDR_Sig <- NA
        all_gs_atc_panel_missing$Nom_Sig <- NA
        all_gs_atc_panel_missing$Method  <- j
        all_gs_atc_panel_missing <- all_gs_atc_panel_missing[, names(all_gs_atc_panel)]
        all_gs_atc_all <- rbind(all_gs_atc_all, all_gs_atc_panel_missing)
      }
      all_gs_atc_all <- rbind(all_gs_atc_all, all_gs_atc_panel)
    }
  }
  all_gs_atc_all <- all_gs_atc_all[all_gs_atc_all$Name != 'Placeholder', ]
  # See build_tx_drug_gtable: bail before facet_wrap errors on an empty factor.
  if (nrow(all_gs_atc_all) == 0) return(NULL)

  methods <- c('MAGMA','GCSC','TWAS-GSEA','TWAS-GSEA (non-dir)')
  methods <- methods[methods %in% all_gs_atc_all$Method]
  all_gs_atc_all$Method <- factor(all_gs_atc_all$Method, levels = methods)

  # Shorten long ATC descriptions
  for (i in unique(all_gs_atc_all$Name)) {
    if (nchar(i) > 30) {
      i_new <- paste0(substr(i, 1, 27), '...')
      i_new <- gsub(' \\.\\.\\.', '...', i_new)
      all_gs_atc_all$Name[all_gs_atc_all$Name == i] <- i_new
    }
  }

  if (sort_choice == 'Alphabetical') {
    all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = unique(all_gs_atc_all$Name[rev(order(all_gs_atc_all$Name))]))
  } else if (sort_choice == 'All - Z') {
    all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = rev(unique(rev(all_gs_atc_all$Name[order(all_gs_atc_all$Z, na.last = FALSE)]))))
  } else if (sort_choice == 'TWAS-GSEA - Z') {
    all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = rev(unique(rev(all_gs_atc_all$Name[all_gs_atc_all$Method == 'TWAS-GSEA'][order(all_gs_atc_all$Z[all_gs_atc_all$Method == 'TWAS-GSEA'], na.last = FALSE)]))))
  } else if (sort_choice == 'TWAS-GSEA (non-dir) - Z') {
    all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = rev(unique(rev(all_gs_atc_all$Name[all_gs_atc_all$Method == 'TWAS-GSEA (non-dir)'][order(all_gs_atc_all$Z[all_gs_atc_all$Method == 'TWAS-GSEA (non-dir)'], na.last = FALSE)]))))
  } else if (sort_choice == 'MAGMA - Z') {
    all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = unique(all_gs_atc_all$Name[all_gs_atc_all$Method == 'MAGMA'][order(all_gs_atc_all$Z[all_gs_atc_all$Method == 'MAGMA'], na.last = FALSE)]))
  } else if (sort_choice == 'GCSC - Z') {
    all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = unique(all_gs_atc_all$Name[all_gs_atc_all$Method == 'GCSC'][order(all_gs_atc_all$Z[all_gs_atc_all$Method == 'GCSC'], na.last = FALSE)]))
  }

  group_siz <- do.call(rbind, lapply(methods, function(m)
    data.frame(Group = m, Size = length(unique(all_gs_atc_all$Panel[all_gs_atc_all$Method == m])))))
  group_siz$Size[group_siz$Size < 2] <- 2
  group_siz$Prop <- group_siz$Size / sum(group_siz$Size)
  group_siz$Width <- 4 * group_siz$Prop

  all_gs_atc_all <- data.table::data.table(all_gs_atc_all)
  dir_data_atc <- all_gs_atc_all[Method == 'TWAS-GSEA']
  pos_data_atc <- all_gs_atc_all[Method != 'TWAS-GSEA']
  dir_max_atc <- if (nrow(dir_data_atc) > 0) suppressWarnings(max(abs(dir_data_atc$Z), na.rm = TRUE)) else NA
  pos_max_atc <- if (nrow(pos_data_atc) > 0) suppressWarnings(max(pos_data_atc$Z,      na.rm = TRUE)) else NA
  if (!is.finite(dir_max_atc)) dir_max_atc <- NA
  if (!is.finite(pos_max_atc)) pos_max_atc <- NA

  heatmap <- ggplot2::ggplot(data = all_gs_atc_all, ggplot2::aes(x = Panel, y = Name)) +
    theme_fn()
  if (nrow(dir_data_atc) > 0 && is.finite(dir_max_atc)) {
    heatmap <- heatmap +
      ggplot2::geom_point(data = dir_data_atc, ggplot2::aes(fill = Z), shape = 21, stroke = 0, size = point_size) +
      ggplot2::scale_fill_gradientn(colours = c("#0066FF","#0099FF","#FFFFFF","#FF6666","#FF0000"),
                                    na.value = NA, name = "TWAS-GSEA\nZ-score",
                                    limits = c(-dir_max_atc, dir_max_atc))
  }
  if (nrow(pos_data_atc) > 0 && is.finite(pos_max_atc)) {
    heatmap <- heatmap +
      ggplot2::geom_point(data = pos_data_atc, ggplot2::aes(colour = Z), shape = 16, size = point_size) +
      ggplot2::scale_colour_gradientn(colours = c("#FFFFFF","#00CC66"),
                                      na.value = NA, name = "Enrichment\nZ-score",
                                      limits = c(0, pos_max_atc))
  }
  heatmap <- heatmap +
    ggplot2::geom_point(data = all_gs_atc_all[which(all_gs_atc_all$Nom_Sig == TRUE), ], ggplot2::aes(x = Panel, y = Name), colour = 'black', fill = NA, size = point_size + 1) +
    ggplot2::geom_point(data = all_gs_atc_all[which(all_gs_atc_all$FDR_Sig == TRUE), ], ggplot2::aes(x = Panel, y = Name), colour = 'black', fill = NA, size = point_size + 2, shape = 15)
  if (nrow(dir_data_atc) > 0) {
    heatmap <- heatmap + ggplot2::geom_point(data = dir_data_atc, ggplot2::aes(fill = Z), shape = 21, stroke = 0, size = point_size)
  }
  if (nrow(pos_data_atc) > 0) {
    heatmap <- heatmap + ggplot2::geom_point(data = pos_data_atc, ggplot2::aes(colour = Z), shape = 16, size = point_size)
  }
  heatmap <- heatmap +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(x = '', y = '', title = if (nzchar(title)) title else NULL) +
    ggplot2::facet_wrap(~ Method, nrow = 1, scales = "free_x") +
    ggplot2::theme(text = ggplot2::element_text(size = font_size),
                   plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5,
                                                 l = left_pad_pt, unit = "pt"))

  gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(heatmap))
  for (i in seq_len(nrow(group_siz))) {
    panel_l <- gt$layout$l[grep(paste0('panel-', i, '-1'), gt$layout$name)]
    gt$widths[panel_l] <- group_siz$Width[i] * gt$widths[panel_l]
  }
  reserve_legend_space(gt, panel_h_pt)
}

enrichmentUI <- function(id) {
  ns <- NS(id)

  tabPanel(
    title="Enrichment Analysis",
    br(),
    p("This tab shows enrichment analysis results. Select the tabs below to see results for your desired gene annotations and methods"),
    hr(),
    uiOutput(ns("enrichment_tabs"))
  )
}

enrichmentServer <- function(id, gwas_data, selected_gwas, config_flags,
                              selected_gwas_multi = NULL,
                              comparison_mode = NULL,
                              shared_filters = NULL,
                              comparison_long = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Cross-GWAS compare-mode sub-modules. Registered unconditionally so their
    # outputs exist for the DOM; they self-guard on comparison_mode() and
    # only render content when >=2 GWAS are selected.
    if (!is.null(selected_gwas_multi) && !is.null(comparison_long)) {
      tissue_compare_server("tissue_compare",
                             gwas_data, selected_gwas_multi,
                             shared_filters, comparison_long)
      atc_compare_server("atc_compare",
                          gwas_data, selected_gwas_multi,
                          shared_filters, comparison_long)
    }

    # Column-guide legend for an enrichment results table. Text is centralised
    # here so each table sub-tab just inserts enr_legend("<key>").
    enr_legend <- function(which) {
      lg_p    <- "Enrichment p-value (smaller = stronger enrichment)."
      lg_pfdr <- "Benjamini-Hochberg FDR-adjusted p-value; < 0.05 is the usual significance threshold."
      lg_panel<- "Gene-expression reference panel (tissue or dataset) used for the TWAS."
      lg_est  <- "Enrichment effect size from the gene-set analysis."
      lg_dir  <- paste0(
        "Direction of effect (relative to the trait being analysed): ",
        "'Opposes disease' = the drug's expression signature counteracts the trait's ",
        "predicted signature; 'Matches disease' = it mimics the trait's signature. ",
        "(Column labels use 'disease' generically; interpretation is the same for any trait.)")
      items <- switch(which,
        drug_magma = list(
          "Name" = "Drug name.",
          "ATC Code" = "Anatomical Therapeutic Chemical (ATC) classification code.",
          "BETA" = "Enrichment effect size (MAGMA competitive gene-set test).",
          "SE" = "Standard error of BETA.",
          "P" = lg_p, "P.FDR" = lg_pfdr),
        drug_gcsc = list(
          "Name" = "Drug name.",
          "ATC Code" = "ATC classification code.",
          "Enrichment" = "Gene co-regulation score (GCSC) enrichment effect size.",
          "SE" = "Standard error of the enrichment.",
          "Z" = "Enrichment Z-score (higher = stronger).",
          "P" = lg_p, "P.FDR" = lg_pfdr),
        drug_twas = list(
          "Name" = "Drug name.",
          "Panel" = lg_panel,
          "N Genes" = "Number of target genes tested for this drug.",
          "Estimate" = lg_est,
          "SE" = "Standard error of Estimate.",
          "Direction" = lg_dir,
          "P" = lg_p, "P.FDR" = lg_pfdr,
          "ATC Code" = "ATC classification code.",
          "ATC Description" = "Text label for the ATC class.",
          "ChEMBL" = "ChEMBL database identifier for the drug."),
        atc_pval = list(
          "Name" = "Drug class (ATC code: description).",
          "N Drugs" = "Number of drugs in the class that were tested.",
          "P" = lg_p, "P.FDR" = lg_pfdr),
        atc_twas = list(
          "Name" = "Drug class (ATC code: description).",
          "Panel" = lg_panel,
          "N Drugs" = "Number of drugs in the class.",
          "Estimate" = lg_est,
          "Direction" = lg_dir,
          "P" = lg_p, "P.FDR" = lg_pfdr),
        cmap_drug = list(
          "cmap_name" = "Compound (CMAP drug signature).",
          "cell_iname" = "Cell line the signature was measured in.",
          "pert_itime" = "Treatment duration.",
          "pert_idose" = "Treatment dose.",
          "moa" = "Mechanism of action.",
          "Panel" = lg_panel,
          "Estimate" = lg_est,
          "Direction" = lg_dir,
          "P" = lg_p, "P.FDR" = lg_pfdr),
        cmap_moa = list(
          "MOA" = "Mechanism of action (drugs grouped by target/mechanism).",
          "Cell_Line" = "Cell line the signatures were measured in.",
          "Panel" = lg_panel,
          "Estimate" = lg_est,
          "Direction" = lg_dir,
          "P" = lg_p, "P.FDR" = lg_pfdr),
        tissue = list(
          "Tissue" = "GTEx v8 tissue.",
          "N Gene" = "Number of genes used in the test for this tissue.",
          "BETA" = "Tissue-specific expression enrichment effect size.",
          "SE" = "Standard error of BETA.",
          "P" = lg_p, "P.FDR" = lg_pfdr,
          "Retained" = "TRUE = the tissue stays FDR-significant after conditioning on the other significant tissues.")
      )
      gd_legend(items, heading = "Column guide")
    }

    ########
    # Dynamic tab rendering
    ########

    output$enrichment_tabs <- renderUI({
      req(config_flags())
      cf <- config_flags()

      drug_targetor_available <- any(cf$magma_drugtargetor, cf$gcsc,
                                     cf$twas_gsea_drugtargetor,
                                     cf$twas_gsea_drugtargetor_nondirectional)

      # Pre-compute Drug choices (may be NULL if no drug enrichment was run)
      drug_data <- if (drug_targetor_available) tx_drug_summary_data() else NULL
      drug_methods <- unique(drug_data$Method)
      drug_expr_panels <- unique(drug_data$Panel[grepl('^TWAS-GSEA', drug_data$Method)])

      # Pre-compute ATC choices
      atc_data <- if (drug_targetor_available) tx_atc_summary_data() else NULL
      atc_methods <- unique(atc_data$Method)
      atc_expr_panels <- unique(atc_data$Panel[grepl('^TWAS-GSEA', atc_data$Method)])

      # Build Drug sub-tabs
      drug_tabs <- list(
        tabPanel(title="Summary", br(),
          p("This tab shows a heatmap summarising, for each candidate drug, whether the genes it targets are enriched for association with the trait across every method and reference panel included in the analysis. Use ", tags$b("Filter data"), " to control which results appear and ", tags$b("Plot options"), " to customise or download the figure."),
          hr(),
          tags$details(class = "gd-details",
            tags$summary("Filter data"),
            tags$div(class = "gd-details-body",
              tags$p(class = "gd-details-intro",
                "Choose which drug-enrichment results appear in the heatmap, ",
                "restrict the view to specific drugs or ATC classes, and choose ",
                "how to sort the rows."),
              fluidRow(
                column(4,
                  selectInput(ns("selected_methods_drug"), "Include results from these methods:", choices=drug_methods, selected=drug_methods, multiple=T),
                  selectInput(ns("selected_expr_panels_drug"), "Include expression / splicing panels:", choices=drug_expr_panels, selected=drug_expr_panels, multiple=T)
                ),
                column(4,
                  textInput(ns("drugInput_drug"), "Show only these drugs (comma- or space-separated):"),
                  textInput(ns("atcInput_drug"), "Show only these ATC codes (comma- or space-separated):")
                ),
                column(4,
                  selectInput(ns("selected_sort_drug"), "Sort rows by:", '', multiple = F),
                  radioButtons(ns("conf_only_drug"), "Show FDR-significant drugs only :",
                               choices = c("True" = T, "False" = F), selected = T)
                )
              )
            )
          ),
          tags$details(class = "gd-details",
            tags$summary("Plot options"),
            tags$div(class = "gd-details-body",
              tags$p(class = "gd-details-intro",
                "Customise how the heatmap looks (title, theme, font size, point size) ",
                "and download it as a PNG, PDF, or SVG at the size and resolution you choose."),
              fluidRow(
                column(4,
                  textInput(ns("plot_title_drug"), "Plot title (optional):", value = ""),
                  selectInput(ns("plot_theme_drug"), "Theme:",
                              choices = c("Black & white" = "bw", "Minimal" = "minimal",
                                          "Classic" = "classic", "Light" = "light"),
                              selected = "bw")
                ),
                column(4,
                  sliderInput(ns("plot_font_size_drug"), "Font size (pt):",
                              min = 10, max = 20, value = 14, step = 1),
                  sliderInput(ns("plot_point_size_drug"), "Point size:",
                              min = 2, max = 8, value = 5, step = 1)
                ),
                column(4,
                  selectInput(ns("dl_format_drug"), "Download format:",
                              choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                              selected = "png"),
                  numericInput(ns("dl_width_drug"), "Width (inches):",
                               value = 12, min = 2, max = 40, step = 0.5),
                  numericInput(ns("dl_height_drug"), "Height (inches):",
                               value = 8, min = 2, max = 40, step = 0.5),
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'png'", ns("dl_format_drug")),
                    numericInput(ns("dl_dpi_drug"), "Resolution (DPI, PNG only):",
                                 value = 300, min = 72, max = 600, step = 25)
                  ),
                  downloadButton(ns("download_plot_drug"), "Download plot")
                )
              )
            )
          ),
          br(),
          uiOutput(ns("message_too_large_drug")),
          uiOutput(ns("message_no_drugs_drug")),
          uiOutput(ns("tx_drug_plot.ui"))
        )
      )
      if (cf$magma_drugtargetor) {
        drug_tabs <- c(drug_tabs, list(tabPanel(title="MAGMA", br(),
          p("This tab shows MAGMA drug enrichment results."), hr(), br(),
          fluidRow(column(width=9, dataTableOutput(ns("tx_drug_magma_table")))), enr_legend("drug_magma"), br()
        )))
      }
      if (cf$gcsc) {
        drug_tabs <- c(drug_tabs, list(tabPanel(title="GCSC", br(),
          p("This tab shows GCSC drug enrichment results."), hr(), br(),
          fluidRow(column(width=9, dataTableOutput(ns("tx_drug_gcsc_table")))), enr_legend("drug_gcsc"), br()
        )))
      }
      if (cf$twas_gsea_drugtargetor) {
        drug_tabs <- c(drug_tabs, list(tabPanel(title="TWAS-GSEA", br(),
          p("This tab shows TWAS-GSEA drug enrichment results."), hr(), br(),
          fluidRow(column(width=9, dataTableOutput(ns("tx_drug_twas_gsea_table")))), enr_legend("drug_twas"), br()
        )))
      }
      if (cf$twas_gsea_drugtargetor_nondirectional) {
        drug_tabs <- c(drug_tabs, list(tabPanel(title="TWAS-GSEA (non-directional)", br(),
          p("This tab shows TWAS-GSEA drug enrichment results using the full DrugTargetor gene-set file (no direction of effect; comparable to MAGMA)."), hr(), br(),
          fluidRow(column(width=9, dataTableOutput(ns("tx_drug_twas_gsea_nondir_table")))), enr_legend("drug_twas"), br()
        )))
      }

      # Build ATC sub-tabs
      atc_tabs <- list(
        tabPanel(title="Summary", br(),
          p("This tab shows a heatmap summarising, for each ATC drug class, whether drugs in that class collectively target genes enriched for trait association, across every method and reference panel included in the analysis. Use ", tags$b("Filter data"), " to control which results appear and ", tags$b("Plot options"), " to customise or download the figure."),
          hr(),
          tags$details(class = "gd-details",
            tags$summary("Filter data"),
            tags$div(class = "gd-details-body",
              tags$p(class = "gd-details-intro",
                "Choose which ATC-class enrichment results appear in the heatmap, ",
                "restrict the view to specific ATC codes, and choose how to sort ",
                "the rows."),
              fluidRow(
                column(4,
                  selectInput(ns("selected_methods_atc"), "Include results from these methods:", choices=atc_methods, selected=atc_methods, multiple=T),
                  selectInput(ns("selected_expr_panels_atc"), "Include expression / splicing panels:", choices=atc_expr_panels, selected=atc_expr_panels, multiple=T)
                ),
                column(4,
                  textInput(ns("atcInput_atc"), "Show only these ATC codes (comma- or space-separated):")
                ),
                column(4,
                  selectInput(ns("selected_sort_atc"), "Sort rows by:", '', multiple = F),
                  radioButtons(ns("conf_only_atc"), "Show FDR-significant classes only :",
                               choices = c("True" = T, "False" = F), selected = T)
                )
              )
            )
          ),
          tags$details(class = "gd-details",
            tags$summary("Plot options"),
            tags$div(class = "gd-details-body",
              tags$p(class = "gd-details-intro",
                "Customise how the heatmap looks (title, theme, font size, point size) ",
                "and download it as a PNG, PDF, or SVG at the size and resolution you choose."),
              fluidRow(
                column(4,
                  textInput(ns("plot_title_atc"), "Plot title (optional):", value = ""),
                  selectInput(ns("plot_theme_atc"), "Theme:",
                              choices = c("Black & white" = "bw", "Minimal" = "minimal",
                                          "Classic" = "classic", "Light" = "light"),
                              selected = "bw")
                ),
                column(4,
                  sliderInput(ns("plot_font_size_atc"), "Font size (pt):",
                              min = 10, max = 20, value = 14, step = 1),
                  sliderInput(ns("plot_point_size_atc"), "Point size:",
                              min = 2, max = 8, value = 5, step = 1)
                ),
                column(4,
                  selectInput(ns("dl_format_atc"), "Download format:",
                              choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                              selected = "png"),
                  numericInput(ns("dl_width_atc"), "Width (inches):",
                               value = 12, min = 2, max = 40, step = 0.5),
                  numericInput(ns("dl_height_atc"), "Height (inches):",
                               value = 8, min = 2, max = 40, step = 0.5),
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'png'", ns("dl_format_atc")),
                    numericInput(ns("dl_dpi_atc"), "Resolution (DPI, PNG only):",
                                 value = 300, min = 72, max = 600, step = 25)
                  ),
                  downloadButton(ns("download_plot_atc"), "Download plot")
                )
              )
            )
          ),
          br(),
          uiOutput(ns("message_too_large_atc")),
          uiOutput(ns("message_no_atcs_atc")),
          uiOutput(ns("tx_atc_plot.ui"))
        )
      )
      if (cf$magma_drugtargetor) {
        atc_tabs <- c(atc_tabs, list(tabPanel(title="MAGMA", br(),
          p("This tab shows MAGMA ATC enrichment results."), hr(), br(),
          fluidRow(column(width=6, dataTableOutput(ns("tx_atc_magma_table")))), enr_legend("atc_pval"), br()
        )))
      }
      if (cf$gcsc) {
        atc_tabs <- c(atc_tabs, list(tabPanel(title="GCSC", br(),
          p("This tab shows GCSC ATC enrichment results."), hr(), br(),
          fluidRow(column(width=6, dataTableOutput(ns("tx_atc_gcsc_table")))), enr_legend("atc_pval"), br()
        )))
      }
      if (cf$twas_gsea_drugtargetor) {
        atc_tabs <- c(atc_tabs, list(tabPanel(title="TWAS-GSEA", br(),
          p("This tab shows TWAS-GSEA ATC enrichment results."), hr(), br(),
          fluidRow(column(width=8, dataTableOutput(ns("tx_atc_twas_gsea_table")))), enr_legend("atc_twas"), br()
        )))
      }
      if (cf$twas_gsea_drugtargetor_nondirectional) {
        atc_tabs <- c(atc_tabs, list(tabPanel(title="TWAS-GSEA (non-directional)", br(),
          p("This tab shows TWAS-GSEA ATC enrichment results using the full DrugTargetor gene-set file (no direction of effect; comparable to MAGMA)."), hr(), br(),
          fluidRow(column(width=8, dataTableOutput(ns("tx_atc_twas_gsea_nondir_table")))), enr_legend("atc_twas"), br()
        )))
      }

      # Build CMAP sub-tabs (per-signature drug heatmap + per-MOA heatmap +
      # raw tables). Only assembled when the run produced CMAP results.
      cmap_drug_data <- if (cf$twas_gsea_cmap) build_cmap_drug_summary_data(gwas_data(), selected_gwas()) else NULL
      cmap_moa_data  <- if (cf$twas_gsea_cmap) build_cmap_moa_summary_data(gwas_data(), selected_gwas())  else NULL
      cmap_panels    <- unique(c(cmap_drug_data$Panel, cmap_moa_data$Panel))

      cmap_core_cells <- c("A375", "HA1E", "HCC515", "HT29", "MCF7", "PC3", "VCAP", "HEPG2", "A549")

      cmap_tab <- NULL
      if (cf$twas_gsea_cmap && (!is.null(cmap_drug_data) || !is.null(cmap_moa_data))) {
        cmap_inner <- list()

        cmap_inner <- c(cmap_inner, list(tabPanel(
          title="Drug",
          br(),
          fluidPage(
            sidebarPanel(
              selectInput(ns("selected_expr_panels_cmap_drug"), "Select expression panels", choices=cmap_panels, selected=cmap_panels, multiple=T),
              selectInput(ns("selected_cell_lines_cmap_drug"), "Select cell lines",
                          choices=unique(cmap_drug_data$cell_iname),
                          selected=intersect(cmap_core_cells, unique(cmap_drug_data$cell_iname)),
                          multiple=T),
              radioButtons(ns("conf_only_cmap_drug"), "Show FDR significant only :",
                           choices = c("True" = T, "False" = F), selected = T),
              textInput(ns("drugInput_cmap_drug"), "Search drug (whitespace- or comma-separated):"),
              textInput(ns("moaInput_cmap_drug"),  "Search MOA (whitespace- or comma-separated):"),
              selectInput(ns("selected_sort_cmap_drug"), "Sort by:", c('Z (best panel)', 'Alphabetical'), selected='Z (best panel)')
            ),
            mainPanel(
              uiOutput(ns("message_too_large_cmap_drug")),
              uiOutput(ns("message_no_cmap_drug")),
              uiOutput(ns("tx_cmap_drug_plot.ui"))
            )
          )
        )))

        cmap_inner <- c(cmap_inner, list(tabPanel(
          title="MOA",
          br(),
          fluidPage(
            sidebarPanel(
              selectInput(ns("selected_expr_panels_cmap_moa"), "Select expression panels", choices=cmap_panels, selected=cmap_panels, multiple=T),
              selectInput(ns("selected_cell_lines_cmap_moa"), "Select cell lines",
                          choices=unique(cmap_moa_data$Cell_Line),
                          selected=intersect(cmap_core_cells, unique(cmap_moa_data$Cell_Line)),
                          multiple=T),
              radioButtons(ns("conf_only_cmap_moa"), "Show FDR significant only :",
                           choices = c("True" = T, "False" = F), selected = T),
              textInput(ns("moaInput_cmap_moa"), "Search MOA (whitespace- or comma-separated):"),
              selectInput(ns("selected_sort_cmap_moa"), "Sort by:", c('Z (best panel)', 'Alphabetical'), selected='Z (best panel)')
            ),
            mainPanel(
              uiOutput(ns("message_too_large_cmap_moa")),
              uiOutput(ns("message_no_cmap_moa")),
              uiOutput(ns("tx_cmap_moa_plot.ui"))
            )
          )
        )))

        cmap_inner <- c(cmap_inner, list(tabPanel(
          title="Drug table", br(),
          p("All per-signature CMAP TWAS-GSEA results."), hr(), br(),
          fluidRow(column(width=12, dataTableOutput(ns("tx_cmap_drug_table")))), enr_legend("cmap_drug"), br()
        )))

        cmap_inner <- c(cmap_inner, list(tabPanel(
          title="MOA table", br(),
          p("All per-MOA CMAP TWAS-GSEA enrichment results."), hr(), br(),
          fluidRow(column(width=10, dataTableOutput(ns("tx_cmap_moa_table")))), enr_legend("cmap_moa"), br()
        )))

        cmap_tab <- do.call(tabPanel,
                            c(list(title="CMAP", br(),
                                   p("Drug repurposing using TWAS-GSEA against reprocessed CMAP level5 drug signatures. Each compound was assayed in multiple cell lines, durations and doses, so per-signature results live under the 'Drug' subtab; per-mechanism aggregation (computed separately per cell line) lives under 'MOA'.")),
                              list(do.call(tabsetPanel, cmap_inner))))
      }

      # Build Tissue tab (MAGMA tissue-specific enrichment)
      # Compare-mode swap: when >=2 GWAS are selected, replace the tissue tab
      # body with the cross-GWAS heatmap. The single-GWAS UI is unchanged.
      in_compare <- !is.null(comparison_mode) && isTRUE(comparison_mode())

      tissue_tab <- NULL
      if (cf$tissue_magma && in_compare) {
        tissue_tab <- tabPanel(
          title = "Tissue", br(),
          tissue_compare_ui(NS(ns("tissue_compare")))
        )
      } else if (cf$tissue_magma) {
        tissue_data <- build_tissue_data(gwas_data(), selected_gwas())
        if (!is.null(tissue_data) && nrow(tissue_data) > 0) {
          tissue_tab <- tabPanel(
            title="Tissue", br(),
            p("MAGMA tissue-specific enrichment across GTEx v8 tissues. Each dot is a tissue; further right means stronger enrichment. Filled teal dots are FDR-significant; a bold label with * means the tissue is retained in the conditional analysis. Use ", tags$b("Filter data"), " to restrict which tissues appear and ", tags$b("Plot options"), " to customise or download the figure."),
            hr(),
            tags$details(class = "gd-details",
              tags$summary("Filter data"),
              tags$div(class = "gd-details-body",
                tags$p(class = "gd-details-intro",
                  "Restrict which tissues appear in the plot and how they are ordered."),
                fluidRow(
                  column(4,
                    selectInput(ns("selected_tissues_tissue"),
                                "Include these tissues:",
                                choices = sort(unique(tissue_data$Tissue)),
                                selected = sort(unique(tissue_data$Tissue)),
                                multiple = TRUE)
                  ),
                  column(4,
                    selectInput(ns("sort_tissue"), "Sort tissues by:",
                                choices = c("Significance" = "significance",
                                            "Alphabetical" = "alphabetical"),
                                selected = "significance")
                  ),
                  column(4,
                    radioButtons(ns("conf_only_tissue"),
                                 "Show FDR-significant tissues only :",
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
                    textInput(ns("plot_title_tissue"), "Plot title (optional):", value = ""),
                    selectInput(ns("plot_theme_tissue"), "Theme:",
                                choices = c("Black & white" = "bw", "Minimal" = "minimal",
                                            "Classic" = "classic", "Light" = "light"),
                                selected = "bw")
                  ),
                  column(4,
                    sliderInput(ns("plot_font_size_tissue"), "Font size (pt):",
                                min = 8, max = 20, value = 13, step = 1),
                    sliderInput(ns("plot_point_size_tissue"), "Point size:",
                                min = 1, max = 8, value = 3, step = 1)
                  ),
                  column(4,
                    selectInput(ns("dl_format_tissue"), "Download format:",
                                choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                                selected = "png"),
                    numericInput(ns("dl_width_tissue"), "Width (inches):",
                                 value = 8, min = 2, max = 40, step = 0.5),
                    numericInput(ns("dl_height_tissue"), "Height (inches):",
                                 value = 9, min = 2, max = 40, step = 0.5),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'png'", ns("dl_format_tissue")),
                      numericInput(ns("dl_dpi_tissue"), "Resolution (DPI, PNG only):",
                                   value = 300, min = 72, max = 600, step = 25)
                    ),
                    downloadButton(ns("download_plot_tissue"), "Download plot")
                  )
                )
              )
            ),
            br(),
            tags$div(style = "max-width: 900px;",
              plotOutput(ns("tx_tissue_plot"), height = "700px")
            ),
            gd_legend(list(
              "X-axis" = "-log10(p-value) for tissue-specific expression enrichment; further right = stronger.",
              "Y-axis" = "GTEx v8 tissue, ordered by significance or alphabetically (chosen in Filter data).",
              "Filled teal dot" = "FDR-significant (P.FDR < 0.05).",
              "Bold label with *" = "Retained in the conditional analysis (still significant after conditioning on the other significant tissues).",
              "Dashed / dotted vertical lines" = "Nominal significance (p = 0.05) and the Bonferroni threshold."
            ), heading = "How to read this plot"),
            br(),
            fluidRow(column(width = 8, dataTableOutput(ns("tx_tissue_table")))),
            enr_legend("tissue"),
            br()
          )
        }
      }

      outer_tabs <- list()
      if (!is.null(tissue_tab)) outer_tabs <- c(outer_tabs, list(tissue_tab))
      if (drug_targetor_available) {
        # Compare-mode swap for the ATC inner tab. Drug-level compare view is
        # deferred to a later phase, so the Drug inner tab keeps single-GWAS
        # behaviour with a small "Viewing:" caption above the existing content.
        atc_body <- if (in_compare) {
          tagList(
            tags$p(style = "color: var(--gd-text-mute); margin-bottom: 6px;",
                    "Cross-GWAS comparison view. Adjust k, threshold and basis at the top of the page."),
            atc_compare_ui(NS(ns("atc_compare")))
          )
        } else {
          do.call(tabsetPanel, atc_tabs)
        }
        drug_body <- if (in_compare) {
          tagList(
            tags$p(style = "color: var(--gd-text-mute);",
                    "Showing single-GWAS view for ", tags$b(selected_gwas()),
                    ". Drug-level cross-GWAS comparison arrives in a later phase."),
            do.call(tabsetPanel, drug_tabs)
          )
        } else {
          do.call(tabsetPanel, drug_tabs)
        }
        drug_targetor_inner <- list(
          tabPanel(title = "Drug", br(), drug_body),
          tabPanel(title = "ATC",  br(), atc_body)
        )
        outer_tabs <- c(outer_tabs, list(
          do.call(tabPanel, c(list(title="Drug Targetor", br()), list(do.call(tabsetPanel, drug_targetor_inner))))
        ))
      }
      if (!is.null(cmap_tab)) outer_tabs <- c(outer_tabs, list(cmap_tab))

      do.call(tabsetPanel, outer_tabs)
    })

    #######
    # Prepare data for drug-specific association tables
    #######

    output$tx_drug_magma_table<-renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )

      tmp<-gd_read(gwas_data(), selected_gwas(), "tx/drug")$magma
      if(is.null(tmp)) return(NULL)
      tmp$BETA<-round(tmp$BETA,3)
      tmp$SE<-round(tmp$SE,3)

      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 4:5)
          )),
        escape = FALSE
      )
    })

    output$tx_drug_gcsc_table<-renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )

      tmp<-gd_read(gwas_data(), selected_gwas(), "tx/drug")$gcsc
      if(is.null(tmp)) return(NULL)
      tmp$Enrichment<-round(tmp$Enrichment, 3)
      tmp$SE<-round(tmp$SE, 3)
      tmp$Z<-round(tmp$Z, 3)

      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 4:5)
          )),
        escape = FALSE
      )
    })

    render_twas_gsea_drug_table <- function(slot){
      req(gwas_data(), selected_gwas())
      # Column order (0-indexed for the JS callback): Name, Panel, N Genes,
      # Estimate, SE, Direction, P, P.FDR, ATC Code, ATC Description, ChEMBL.
      # P (6) and P.FDR (7) are rendered in scientific notation.
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[6];",
        "  $('td:eq(6)', row).html(x.toExponential(2));",
        "  var y = data[7];",
        "  $('td:eq(7)', row).html(y.toExponential(2));",
        "}"
      )

      tmp<-gd_read(gwas_data(), selected_gwas(), "tx/drug")[[slot]]
      if(is.null(tmp)) return(NULL)
      tmp$Estimate<-round(tmp$Estimate, 3)
      tmp$SE<-round(tmp$SE, 3)

      # Explicit column selection: show Direction alongside Estimate, hide
      # Reversal_Z (used for plotting only; Direction conveys the same info
      # in a more human-readable form).
      keep <- intersect(c('Name','Panel','N Genes','Estimate','SE','Direction',
                          'P','P.FDR','ATC Code','ATC Description','ChEMBL'),
                        names(tmp))
      tmp <- tmp[, keep, with = FALSE]

      datatable(
        tmp,
        rownames = F,
        options = list(
          rowCallback = JS(js),
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 6:7)
          )),
        escape = FALSE
      )
    }

    output$tx_drug_twas_gsea_table<-renderDataTable({ render_twas_gsea_drug_table('twas_gsea') })
    output$tx_drug_twas_gsea_nondir_table<-renderDataTable({ render_twas_gsea_drug_table('twas_gsea_nondir') })

    #######
    # Prepare data for drug enrichment plot
    #######

    tx_drug_summary_data <- reactive({
      req(gwas_data(), selected_gwas())
      build_drug_summary_data(gwas_data(), selected_gwas())
    })

    tx_drug_summary_data_filtered<-reactive({
      req(input$conf_only_drug, input$selected_methods_drug)
      all_gs<-tx_drug_summary_data()

      # Filter results table by user specified methods
      all_gs<-all_gs[all_gs$Method %in% input$selected_methods_drug,]

      # Filter results table by user specified expression (applies to both directional and non-directional TWAS-GSEA)
      gsea_rows <- grepl('^TWAS-GSEA', all_gs$Method)
      if(any(gsea_rows)){
        all_gs<-all_gs[!(gsea_rows & !(all_gs$Panel %in% input$selected_expr_panels_drug)),]
      }

      # Insert NA rows for all panels and methods so when filtering by drug, all selected panels and methods remain
      na_rows<-all_gs[!(duplicated(paste0(all_gs$Panel, all_gs$Method))),]
      na_rows$Name<-'Placeholder'
      na_rows$`ATC Code`<-'Placeholder'
      na_rows$Z<-NA
      na_rows$P<-NA
      na_rows$P.FDR<-NA

      all_gs<-rbind(na_rows, all_gs)

      # Filter results table if user specifies high confidence genes only
      if(input$conf_only_drug){
        hc_drugs<-all_gs$Name[
          which(
            (all_gs$Method == 'TWAS-GSEA' & all_gs$P.FDR < 0.05) |
              (all_gs$Method == 'TWAS-GSEA (non-dir)' & all_gs$P.FDR < 0.05 & all_gs$Z > 0) |
              (all_gs$Method == 'MAGMA' & all_gs$P.FDR < 0.05 & all_gs$Z > 0) |
              (all_gs$Method == 'GCSC' & all_gs$P.FDR < 0.05 & all_gs$Z > 0)
          )]

        all_gs<-all_gs[all_gs$Name %in% c(hc_drugs,'Placeholder'),]
      }

      input_drugs <- unlist(strsplit(input$drugInput_drug, "[, ]"))
      selected_drugs <- input_drugs[input_drugs != ""]

      input_atc <- unlist(strsplit(input$atcInput_drug, "[, ]"))
      selected_atc <- input_atc[input_atc != ""]

      if(length(selected_drugs) > 0){
        if(sum(grepl(paste(selected_drugs, collapse='|'), all_gs$Name, ignore.case = T)) > 0){
          selected_drugs<-c(selected_drugs, 'Placeholder')
          all_gs<-all_gs[grepl(paste(selected_drugs, collapse='|'), all_gs$Name, ignore.case = T) & !is.na(all_gs$Name),]
        } else {
          all_gs<-data.frame(matrix(nrow=0, ncol=5))
        }
      }

      if(length(selected_atc) > 0){
        if(sum(grepl(paste(selected_atc, collapse='|'), all_gs$`ATC Code`, ignore.case = T)) > 0){
          selected_atc<-c(selected_atc, 'Placeholder')
          all_gs<-all_gs[grepl(paste(selected_atc, collapse='|'), all_gs$`ATC Code`, ignore.case = T) & !is.na(all_gs$`ATC Code`),]
        } else {
          all_gs<-data.frame(matrix(nrow=0, ncol=5))
        }
      }

      return(all_gs)
    }) %>% debounce(1000)

    # Identify number of drugs
    plot_dim_drug<-reactive({
      all_gs<-tx_drug_summary_data_filtered()
      # `min_height = 350` reserves overall device height. `min_panel_h_pt
      # = 200` floors the panel row so the two stacked colourbar legends
      # (centred on the panel) don't clip even when only a handful of
      # drugs are visible — row spacing spreads slightly in that
      # edge case, in exchange for a fully-visible legend.
      calc_plot_dims(all_gs, y_col = "Name", x_col = "Panel", facet_col = "Method",
                     font_size = input$plot_font_size_drug %||% 14,
                     min_height = 350, min_panel_h_pt = 200)
    })

    observeEvent(tx_drug_summary_data_filtered(), {
      tmp<-tx_drug_summary_data_filtered()
      choices<-'All - Z'
      if(any(tmp$Method == 'MAGMA')){
        choices<-c(choices, 'MAGMA - Z')
      }
      if(any(tmp$Method == 'GCSC')){
        choices<-c(choices, 'GCSC - Z')
      }
      if(any(tmp$Method == 'TWAS-GSEA')){
        choices<-c(choices, 'TWAS-GSEA - Z')
      }
      if(any(tmp$Method == 'TWAS-GSEA (non-dir)')){
        choices<-c(choices, 'TWAS-GSEA (non-dir) - Z')
      }
      if(length(unique(tmp$Method)) == 1){
        choices<-choices[choices != 'All - Z']
      }
      choices<-c(choices, 'Alphabetical')

      updateSelectInput(session, "selected_sort_drug", choices = choices, selected=choices[1])
    })

    output$tx_drug_plot <- renderPlot({
      all_gs <- tx_drug_summary_data_filtered()
      if (plot_dim_drug()[['height']] >= 10000 || nrow(all_gs) == 0) return(NULL)
      gt <- build_tx_drug_gtable(
        all_gs = all_gs,
        sort_choice = input$selected_sort_drug %||% "Alphabetical",
        font_size = input$plot_font_size_drug %||% 14,
        point_size = input$plot_point_size_drug %||% 5,
        theme_fn = enr_theme_fn(input$plot_theme_drug),
        title = input$plot_title_drug %||% "",
        panel_h_pt = plot_dim_drug()[['panel_h_pt']],
        left_pad_pt = plot_dim_drug()[['left_pad_pt']]
      )
      if (!is.null(gt)) grid::grid.draw(gt)
    })

    # Keep download width/height in sync with the on-screen plot so the
    # default download matches the current window. User overrides get
    # replaced whenever the on-screen dims change.
    observeEvent(plot_dim_drug(), {
      dims <- plot_dim_drug()
      if (dims$width > 100 && dims$height < 10000) {
        updateNumericInput(session, "dl_width_drug",  value = round(dims$width  / 72, 1))
        updateNumericInput(session, "dl_height_drug", value = round(dims$height / 72, 1))
      }
    })

    output$download_plot_drug <- downloadHandler(
      filename = function() {
        sprintf("drug_enrichment_%s.%s",
                format(Sys.time(), "%Y%m%d_%H%M%S"),
                input$dl_format_drug %||% "png")
      },
      content = function(file) {
        gt <- build_tx_drug_gtable(
          all_gs = tx_drug_summary_data_filtered(),
          sort_choice = input$selected_sort_drug %||% "Alphabetical",
          font_size = input$plot_font_size_drug %||% 14,
          point_size = input$plot_point_size_drug %||% 5,
          theme_fn = enr_theme_fn(input$plot_theme_drug),
          title = input$plot_title_drug %||% "",
          panel_h_pt = plot_dim_drug()[['panel_h_pt']],
          left_pad_pt = plot_dim_drug()[['left_pad_pt']]
        )
        w <- input$dl_width_drug  %||% 12
        h <- input$dl_height_drug %||% 8
        fmt <- input$dl_format_drug %||% "png"
        switch(fmt,
          png = grDevices::png(file, width = w, height = h, units = "in",
                               res = input$dl_dpi_drug %||% 300),
          pdf = grDevices::pdf(file, width = w, height = h),
          svg = grDevices::svg(file, width = w, height = h)
        )
        if (!is.null(gt)) grid::grid.draw(gt)
        grDevices::dev.off()
      }
    )

    output$tx_drug_plot.ui <- renderUI({
      filtered <- tx_drug_summary_data_filtered()
      has_real <- any(filtered$Name != 'Placeholder')
      if(plot_dim_drug()[['height']] < 10000 && has_real){
        tagList(
          plotOutput(ns("tx_drug_plot"), height = plot_dim_drug()[['height']], width=plot_dim_drug()[['width']]),
          tags$div(style = "max-width: 800px; font-size: 0.9em; margin-top: 10px;",
            tags$b("Figure legend. "),
            "Each point shows a drug-enrichment Z-score for one method/expression-panel combination. ",
            tags$b("Directional TWAS-GSEA"), " uses a diverging palette: ",
            tags$span(style = "color:#FF0000;", tags$b("red")), " indicates the drug-target signature is ",
            tags$i("opposite"), " to the trait's TWAS signature, and ",
            tags$span(style = "color:#0066FF;", tags$b("blue")), " indicates the drug-target signature ",
            tags$i("matches"), " the trait's signature. ",
            tags$b("MAGMA, GCSC, and non-directional TWAS-GSEA"), " use a sequential white-to-",
            tags$span(style = "color:#00CC66;", tags$b("green")),
            " palette: green indicates that genes targeted by the drug are enriched for the trait's association signal (no direction-of-effect interpretation). ",
            "Hollow black circles mark nominally significant results (P < 0.05); black squares mark FDR-significant results (FDR < 0.05)."
          )
        )
      } else {
        NULL
      }
    })

    output$message_too_large_drug <- renderUI({
      if(plot_dim_drug()[['height']] > 10000){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "Plot is too large. Restrict to significant drugs or specify a list of drugs."
        ))
      }
    })

    output$message_no_drugs_drug <- renderUI({
      filtered <- tx_drug_summary_data_filtered()
      real_drugs <- filtered$Name[filtered$Name != 'Placeholder']
      if (length(real_drugs) == 0 && isTRUE(as.logical(input$conf_only_drug))) {
        return(HTML(paste0(
          "<div style='color: #666; font-style: italic; margin: 10px 0;'>",
          "No FDR-significant drugs match your current filters. ",
          "Set <b>Show FDR-significant drugs only</b> to <b>False</b> to see all results, ",
          "or broaden the method / panel / drug / ATC filters.",
          "</div>"
        )))
      }
      if (length(real_drugs) == 0) {
        return(HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "No drugs match your current filters. Try broadening the method / panel / drug / ATC filters."
        )))
      }
    })

    #######
    # Prepare data for atc-specific association tables
    #######

    output$tx_atc_magma_table<-renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[2];",
        "  $('td:eq(2)', row).html(x.toExponential(2));",
        "  var y = data[3];",
        "  $('td:eq(3)', row).html(y.toExponential(2));",
        "}"
      )

      tmp<-gd_read(gwas_data(), selected_gwas(), "tx/atc")$magma
      if(is.null(tmp)) return(NULL)
      tmp$Name<-paste0(tmp$`ATC Code`,': ',tmp$`ATC Description`)
      tmp<-tmp[,c('Name','N Drugs','P','P.FDR'), with=F]

      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = 0:3),
            list(width = '60px', targets = 2:3)
          )),
        escape = FALSE
      )
    })

    output$tx_atc_gcsc_table<-renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[2];",
        "  $('td:eq(2)', row).html(x.toExponential(2));",
        "  var y = data[3];",
        "  $('td:eq(3)', row).html(y.toExponential(2));",
        "}"
      )

      tmp<-gd_read(gwas_data(), selected_gwas(), "tx/atc")$gcsc
      if(is.null(tmp)) return(NULL)
      tmp$Name<-paste0(tmp$`ATC Code`,': ',tmp$`ATC Description`)
      tmp<-tmp[,c('Name','N Drugs','P','P.FDR'), with=F]

      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = 0:3),
            list(width = '60px', targets = 2:3)
          )),
        escape = FALSE
      )
    })

    render_twas_gsea_atc_table <- function(slot){
      req(gwas_data(), selected_gwas())
      tmp<-gd_read(gwas_data(), selected_gwas(), "tx/atc")[[slot]]
      if(is.null(tmp)) return(NULL)

      tmp$Name<-paste0(tmp$`ATC Code`,': ',tmp$`ATC Description`)
      tmp$Estimate<-round(tmp$Estimate,3)

      # Column order: Name, Panel, N Drugs, Estimate, Direction, P, P.FDR.
      # Direction is supplied by the read function (Direction = "Opposes
      # disease" / "Matches disease" / NA). Reversal_Z is used by the heatmap
      # only and is dropped from this table to keep it scannable.
      tmp<-tmp[,c("Name","Panel","N Drugs","Estimate","Direction","P","P.FDR"), with=F]

      # JS callback: P (5) and P.FDR (6) in scientific notation.
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[5];",
        "  $('td:eq(5)', row).html(x.toExponential(2));",
        "  var y = data[6];",
        "  $('td:eq(6)', row).html(y.toExponential(2));",
        "}"
      )

      datatable(
        tmp,
        rownames = F,
        options = list(
          rowCallback = JS(js),
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 5:6)
          )),
        escape = FALSE
      )
    }

    output$tx_atc_twas_gsea_table<-renderDataTable({ render_twas_gsea_atc_table('twas_gsea') })
    output$tx_atc_twas_gsea_nondir_table<-renderDataTable({ render_twas_gsea_atc_table('twas_gsea_nondir') })

    #######
    # Prepare data for atc enrichment plot
    #######

    tx_atc_summary_data <- reactive({
      req(gwas_data(), selected_gwas())
      build_atc_summary_data(gwas_data(), selected_gwas())
    })

    tx_atc_summary_data_filtered<-reactive({
      req(input$conf_only_atc, input$selected_methods_atc)
      all_gs_atc<-tx_atc_summary_data()

      # Filter results table by user specified methods
      all_gs_atc<-all_gs_atc[all_gs_atc$Method %in% input$selected_methods_atc,]

      # Filter results table by user specified expression
      gsea_rows_atc <- grepl('^TWAS-GSEA', all_gs_atc$Method)
      if(any(gsea_rows_atc)){
        all_gs_atc<-all_gs_atc[!(gsea_rows_atc & !(all_gs_atc$Panel %in% input$selected_expr_panels_atc)),]
      }

      # Insert NA rows for all panels and methods so when filtering by atc, all selected panels and methods remain
      na_rows<-all_gs_atc[!(duplicated(paste0(all_gs_atc$Panel, all_gs_atc$Method))),]
      na_rows$Name<-'Placeholder'
      na_rows$Z<-NA
      na_rows$FDR_Sig<-NA
      na_rows$Nom_Sig<-NA

      all_gs_atc<-rbind(na_rows, all_gs_atc)

      # Filter results table if user specifies high confidence genes only
      if(input$conf_only_atc){
        hc_atc<-all_gs_atc$Name[
          which(
            (all_gs_atc$Method == 'TWAS-GSEA' & all_gs_atc$FDR_Sig) |
              (all_gs_atc$Method == 'TWAS-GSEA (non-dir)' & all_gs_atc$FDR_Sig & all_gs_atc$Z > 0) |
              (all_gs_atc$Method == 'MAGMA' & all_gs_atc$FDR_Sig & all_gs_atc$Z > 0) |
              (all_gs_atc$Method == 'GCSC' & all_gs_atc$FDR_Sig & all_gs_atc$Z > 0)
          )]

        all_gs_atc<-all_gs_atc[all_gs_atc$Name %in% c(hc_atc,'Placeholder'),]
      }

      input_atcs <- unlist(strsplit(input$atcInput_atc, "[, ]"))
      selected_atcs <- input_atcs[input_atcs != ""]

      if(length(selected_atcs) > 0){
        if(sum(grepl(paste(selected_atcs, collapse='|'), all_gs_atc$Name, ignore.case = T)) > 0){
          selected_atcs<-c(selected_atcs, 'Placeholder')
          all_gs_atc<-all_gs_atc[grepl(paste(selected_atcs, collapse='|'), all_gs_atc$Name, ignore.case = T) & !is.na(all_gs_atc$Name),]
        } else {
          all_gs_atc<-data.frame(matrix(nrow=0, ncol=5))
        }
      }

      return(all_gs_atc)
    }) %>% debounce(1000)

    # Identify number of atcs
    plot_dim_atc<-reactive({
      all_gs_atc<-tx_atc_summary_data_filtered()
      # Truncate long ATC names to match what the plot displays
      long <- nchar(all_gs_atc$Name) > 30
      all_gs_atc$Name[long] <- paste0(substr(all_gs_atc$Name[long], 1, 27), '...')
      # Same rationale as plot_dim_drug — reserve device height for the
      # stacked colourbar legends and floor the panel row so they don't
      # clip when only a few ATC classes are visible.
      calc_plot_dims(all_gs_atc, y_col = "Name", x_col = "Panel", facet_col = "Method",
                     font_size = input$plot_font_size_atc %||% 14,
                     min_height = 350, min_panel_h_pt = 200)
    })

    observeEvent(tx_atc_summary_data_filtered(), {
      tmp<-tx_atc_summary_data_filtered()
      choices<-'All - Z'
      if(any(tmp$Method == 'MAGMA')){
        choices<-c(choices, 'MAGMA - Z')
      }
      if(any(tmp$Method == 'GCSC')){
        choices<-c(choices, 'GCSC - Z')
      }
      if(any(tmp$Method == 'TWAS-GSEA')){
        choices<-c(choices, 'TWAS-GSEA - Z')
      }
      if(any(tmp$Method == 'TWAS-GSEA (non-dir)')){
        choices<-c(choices, 'TWAS-GSEA (non-dir) - Z')
      }
      if(length(unique(tmp$Method)) == 1){
        choices<-choices[choices != 'All - Z']
      }
      choices<-c(choices, 'Alphabetical')

      updateSelectInput(session, "selected_sort_atc", choices = choices, selected=choices[1])
    })

    output$tx_atc_plot <- renderPlot({
      all_gs_atc <- tx_atc_summary_data_filtered()
      if (plot_dim_atc()[['height']] >= 10000 || nrow(all_gs_atc) == 0) return(NULL)
      gt <- build_tx_atc_gtable(
        all_gs_atc = all_gs_atc,
        sort_choice = input$selected_sort_atc %||% "Alphabetical",
        font_size = input$plot_font_size_atc %||% 14,
        point_size = input$plot_point_size_atc %||% 5,
        theme_fn = enr_theme_fn(input$plot_theme_atc),
        title = input$plot_title_atc %||% "",
        panel_h_pt = plot_dim_atc()[['panel_h_pt']],
        left_pad_pt = plot_dim_atc()[['left_pad_pt']]
      )
      if (!is.null(gt)) grid::grid.draw(gt)
    })

    observeEvent(plot_dim_atc(), {
      dims <- plot_dim_atc()
      if (dims$width > 100 && dims$height < 10000) {
        updateNumericInput(session, "dl_width_atc",  value = round(dims$width  / 72, 1))
        updateNumericInput(session, "dl_height_atc", value = round(dims$height / 72, 1))
      }
    })

    output$download_plot_atc <- downloadHandler(
      filename = function() {
        sprintf("atc_enrichment_%s.%s",
                format(Sys.time(), "%Y%m%d_%H%M%S"),
                input$dl_format_atc %||% "png")
      },
      content = function(file) {
        gt <- build_tx_atc_gtable(
          all_gs_atc = tx_atc_summary_data_filtered(),
          sort_choice = input$selected_sort_atc %||% "Alphabetical",
          font_size = input$plot_font_size_atc %||% 14,
          point_size = input$plot_point_size_atc %||% 5,
          theme_fn = enr_theme_fn(input$plot_theme_atc),
          title = input$plot_title_atc %||% "",
          panel_h_pt = plot_dim_atc()[['panel_h_pt']],
          left_pad_pt = plot_dim_atc()[['left_pad_pt']]
        )
        w <- input$dl_width_atc  %||% 12
        h <- input$dl_height_atc %||% 8
        fmt <- input$dl_format_atc %||% "png"
        switch(fmt,
          png = grDevices::png(file, width = w, height = h, units = "in",
                               res = input$dl_dpi_atc %||% 300),
          pdf = grDevices::pdf(file, width = w, height = h),
          svg = grDevices::svg(file, width = w, height = h)
        )
        if (!is.null(gt)) grid::grid.draw(gt)
        grDevices::dev.off()
      }
    )

    output$tx_atc_plot.ui <- renderUI({
      filtered <- tx_atc_summary_data_filtered()
      has_real <- any(filtered$Name != 'Placeholder')
      if(plot_dim_atc()[['height']] < 10000 && has_real){
        tagList(
          plotOutput(ns("tx_atc_plot"), height = plot_dim_atc()[['height']], width=plot_dim_atc()[['width']]),
          tags$div(style = "max-width: 800px; font-size: 0.9em; margin-top: 10px;",
            tags$b("Figure legend. "),
            "Each point shows an ATC-class enrichment Z-score for one method/expression-panel combination. ",
            tags$b("Directional TWAS-GSEA"), " uses a diverging palette: ",
            tags$span(style = "color:#FF0000;", tags$b("red")), " indicates drugs in the ATC class collectively ",
            tags$i("oppose"), " the trait's TWAS signature, and ",
            tags$span(style = "color:#0066FF;", tags$b("blue")), " indicates the class collectively ",
            tags$i("matches"), " the trait's signature. ",
            tags$b("MAGMA, GCSC, and non-directional TWAS-GSEA"), " use a sequential white-to-",
            tags$span(style = "color:#00CC66;", tags$b("green")),
            " palette: green indicates that drugs in the ATC class are enriched for the trait's association signal (no direction-of-effect interpretation). ",
            "Hollow black circles mark nominally significant results (P < 0.05); black squares mark FDR-significant results (FDR < 0.05)."
          )
        )
      } else {
        NULL
      }
    })

    output$message_too_large_atc <- renderUI({
      if(plot_dim_atc()[['height']] > 10000){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "Plot is too large. Restrict to significant atcs or specify a list of atcs."
        ))
      }
    })

    output$message_no_atcs_atc <- renderUI({
      filtered <- tx_atc_summary_data_filtered()
      real_atcs <- filtered$Name[filtered$Name != 'Placeholder']
      if (length(real_atcs) == 0 && isTRUE(as.logical(input$conf_only_atc))) {
        return(HTML(paste0(
          "<div style='color: #666; font-style: italic; margin: 10px 0;'>",
          "No FDR-significant ATC classes match your current filters. ",
          "Set <b>Show FDR-significant classes only</b> to <b>False</b> to see all results, ",
          "or broaden the method / panel / ATC-code filters.",
          "</div>"
        )))
      }
      if (length(real_atcs) == 0) {
        return(HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "No ATC classes match your current filters. Try broadening the method / panel / ATC-code filters."
        )))
      }
    })

    #######
    # CMAP TWAS-GSEA enrichment (per-signature drug + per-MOA)
    #######

    tx_cmap_drug_data <- reactive({
      req(gwas_data(), selected_gwas())
      build_cmap_drug_summary_data(gwas_data(), selected_gwas())
    })

    tx_cmap_moa_data <- reactive({
      req(gwas_data(), selected_gwas())
      build_cmap_moa_summary_data(gwas_data(), selected_gwas())
    })

    # ----- per-MOA -----
    tx_cmap_moa_filtered <- reactive({
      d <- tx_cmap_moa_data()
      if(is.null(d) || nrow(d) == 0) return(d)
      if(!is.null(input$selected_expr_panels_cmap_moa) && length(input$selected_expr_panels_cmap_moa) > 0)
        d <- d[d$Panel %in% input$selected_expr_panels_cmap_moa, ]
      if(!is.null(input$selected_cell_lines_cmap_moa) && length(input$selected_cell_lines_cmap_moa) > 0)
        d <- d[d$Cell_Line %in% input$selected_cell_lines_cmap_moa, ]
      if(isTRUE(as.logical(input$conf_only_cmap_moa))){
        keep_names <- unique(d$Name[d$FDR_Sig %in% TRUE])
        d <- d[d$Name %in% keep_names, ]
      }
      input_moas <- unlist(strsplit(input$moaInput_cmap_moa, "[, ]+"))
      selected_moas <- input_moas[input_moas != ""]
      if(length(selected_moas) > 0){
        d <- d[grepl(paste(selected_moas, collapse='|'), d$Name, ignore.case = TRUE) & !is.na(d$Name), ]
      }
      d
    }) %>% debounce(500)

    plot_dim_cmap_moa <- reactive({
      d <- tx_cmap_moa_filtered()
      calc_plot_dims(d, y_col = "Name", x_col = "Panel", facet_col = "Cell_Line")
    })

    # ----- per-signature drug -----
    tx_cmap_drug_filtered <- reactive({
      d <- tx_cmap_drug_data()
      if(is.null(d) || nrow(d) == 0) return(d)
      if(!is.null(input$selected_expr_panels_cmap_drug) && length(input$selected_expr_panels_cmap_drug) > 0)
        d <- d[d$Panel %in% input$selected_expr_panels_cmap_drug, ]
      if(!is.null(input$selected_cell_lines_cmap_drug) && length(input$selected_cell_lines_cmap_drug) > 0)
        d <- d[d$cell_iname %in% input$selected_cell_lines_cmap_drug, ]
      if(isTRUE(as.logical(input$conf_only_cmap_drug))){
        keep_names <- unique(d$Name[d$FDR_Sig %in% TRUE])
        d <- d[d$Name %in% keep_names, ]
      }
      input_drugs <- unlist(strsplit(input$drugInput_cmap_drug, "[, ]+"))
      selected_drugs <- input_drugs[input_drugs != ""]
      if(length(selected_drugs) > 0){
        d <- d[grepl(paste(selected_drugs, collapse='|'), d$Name, ignore.case = TRUE) & !is.na(d$Name), ]
      }
      input_moas <- unlist(strsplit(input$moaInput_cmap_drug, "[, ]+"))
      selected_moas <- input_moas[input_moas != ""]
      if(length(selected_moas) > 0){
        d <- d[grepl(paste(selected_moas, collapse='|'), d$moa, ignore.case = TRUE) & !is.na(d$moa), ]
      }
      d
    }) %>% debounce(500)

    plot_dim_cmap_drug <- reactive({
      d <- tx_cmap_drug_filtered()
      calc_plot_dims(d, y_col = "Name", x_col = "Panel", facet_col = "cell_iname")
    })

    # Shared CMAP heatmap helper. Faceted by cell line, reuses the red/blue
    # diverging palette + nominal/FDR overlays as the directional DrugTargetor
    # TWAS-GSEA plot. Z is already negated upstream so red = drug reverses
    # disease signature.
    cmap_heatmap <- function(d, sort_choice, facet_col = NULL, left_pad_pt = 0){
      if(is.null(d) || nrow(d) == 0) return(NULL)
      best_z <- tapply(d$Z, d$Name, function(x) max(x, na.rm = TRUE))
      best_z[!is.finite(best_z)] <- NA
      if(identical(sort_choice, 'Alphabetical')){
        lvl <- sort(names(best_z), decreasing = TRUE)
      } else {
        lvl <- names(sort(best_z, decreasing = FALSE, na.last = FALSE))
      }
      d$Name <- factor(d$Name, levels = lvl)

      z_max <- max(abs(d$Z), na.rm = TRUE)
      if(!is.finite(z_max) || z_max == 0) z_max <- 1

      p <- ggplot(d, aes(x = Panel, y = Name)) +
        theme_bw() +
        geom_point(aes(fill = Z), shape = 21, stroke = 0, size = 5) +
        scale_fill_gradientn(colours = c("#0066FF","#0099FF","#FFFFFF","#FF6666","#FF0000"),
                             na.value = NA, name = "CMAP\nZ-score",
                             limits = c(-z_max, z_max)) +
        geom_point(data = d[d$Nom_Sig %in% TRUE, ], colour = 'black', fill = NA, size = 6) +
        geom_point(data = d[d$FDR_Sig %in% TRUE, ], colour = 'black', fill = NA, size = 7, shape = 15) +
        geom_point(aes(fill = Z), shape = 21, stroke = 0, size = 5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title  = element_text(hjust = 0.5),
              text        = element_text(size = 14),
              plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = left_pad_pt, unit = "pt")) +
        labs(x = '', y = '')

      if (!is.null(facet_col)) {
        p <- p + facet_wrap(as.formula(paste("~", facet_col)), nrow = 1, scales = "free_x")
      }
      p
    }

    cmap_legend_div <- tags$div(style = "max-width: 800px; font-size: 0.9em; margin-top: 10px;",
      tags$b("Figure legend. "),
      "Each point shows a CMAP TWAS-GSEA Z-score for one expression panel. ",
      tags$span(style = "color:#FF0000;", tags$b("Red")), " = drug signature is ", tags$i("opposite"),
      " to the trait's TWAS signature. ",
      tags$span(style = "color:#0066FF;", tags$b("Blue")), " = drug signature ", tags$i("matches"),
      " the trait. ",
      "Hollow black circles mark nominally significant results (P < 0.05); black squares mark FDR-significant results (FDR < 0.05)."
    )

    output$tx_cmap_moa_plot <- renderPlot({
      d <- tx_cmap_moa_filtered()
      if(is.null(d) || nrow(d) == 0) return(NULL)
      if(plot_dim_cmap_moa()[['height']] >= 10000) return(NULL)
      cmap_heatmap(d, input$selected_sort_cmap_moa, facet_col = "Cell_Line",
                   left_pad_pt = plot_dim_cmap_moa()[['left_pad_pt']])
    })

    output$tx_cmap_moa_plot.ui <- renderUI({
      d <- tx_cmap_moa_filtered()
      if(is.null(d) || nrow(d) == 0) return(NULL)
      if(plot_dim_cmap_moa()[['height']] >= 10000) return(NULL)
      tagList(
        plotOutput(ns("tx_cmap_moa_plot"),
                   height = plot_dim_cmap_moa()[['height']],
                   width  = plot_dim_cmap_moa()[['width']]),
        cmap_legend_div
      )
    })

    output$message_too_large_cmap_moa <- renderUI({
      if(plot_dim_cmap_moa()[['height']] >= 10000)
        HTML("<div style='color: red;'>Plot is too large. Restrict to FDR-significant MOAs or pick specific MOAs in the search box.</div>")
    })

    output$message_no_cmap_moa <- renderUI({
      d <- tx_cmap_moa_filtered()
      if(is.null(d) || nrow(d) == 0)
        HTML("<div style='color: red;'>No MOAs to display. Disable the FDR filter or check that twas_gsea_cmap was run.</div>")
    })

    output$tx_cmap_drug_plot <- renderPlot({
      d <- tx_cmap_drug_filtered()
      if(is.null(d) || nrow(d) == 0) return(NULL)
      if(plot_dim_cmap_drug()[['height']] >= 10000) return(NULL)
      cmap_heatmap(d, input$selected_sort_cmap_drug, facet_col = "cell_iname",
                   left_pad_pt = plot_dim_cmap_drug()[['left_pad_pt']])
    })

    output$tx_cmap_drug_plot.ui <- renderUI({
      d <- tx_cmap_drug_filtered()
      if(is.null(d) || nrow(d) == 0) return(NULL)
      if(plot_dim_cmap_drug()[['height']] >= 10000) return(NULL)
      tagList(
        plotOutput(ns("tx_cmap_drug_plot"),
                   height = plot_dim_cmap_drug()[['height']],
                   width  = plot_dim_cmap_drug()[['width']]),
        cmap_legend_div
      )
    })

    output$message_too_large_cmap_drug <- renderUI({
      if(plot_dim_cmap_drug()[['height']] >= 10000)
        HTML("<div style='color: red;'>Too many signatures to plot. Restrict to FDR-significant rows, or use the search boxes to pick specific drugs / MOAs.</div>")
    })

    output$message_no_cmap_drug <- renderUI({
      d <- tx_cmap_drug_filtered()
      if(is.null(d) || nrow(d) == 0)
        HTML("<div style='color: red;'>No drug signatures to display. Disable the FDR filter or check that twas_gsea_cmap was run.</div>")
    })

    output$tx_cmap_drug_table <- renderDataTable({
      req(gwas_data(), selected_gwas())
      d <- gd_read(gwas_data(), selected_gwas(), "tx/cmap")$drug
      if(is.null(d)) return(NULL)
      # Hide Reversal_Z from the table (it is used by the heatmap; Direction
      # column conveys the same information in human-readable form).
      hide_idx <- which(names(d) == 'Reversal_Z') - 1L
      cdefs <- list(list(className = 'dt-center', targets = '_all'))
      if(length(hide_idx) == 1L) cdefs <- c(cdefs, list(list(visible = FALSE, targets = hide_idx)))
      datatable(d, rownames = FALSE,
                options = list(scrollX = TRUE, columnDefs = cdefs))
    })

    output$tx_cmap_moa_table <- renderDataTable({
      req(gwas_data(), selected_gwas())
      d <- gd_read(gwas_data(), selected_gwas(), "tx/cmap")$moa
      if(is.null(d)) return(NULL)
      hide_idx <- which(names(d) == 'Reversal_Z') - 1L
      cdefs <- list(list(className = 'dt-center', targets = '_all'))
      if(length(hide_idx) == 1L) cdefs <- c(cdefs, list(list(visible = FALSE, targets = hide_idx)))
      datatable(d, rownames = FALSE,
                options = list(scrollX = TRUE, columnDefs = cdefs))
    })

    #######
    # Tissue-specific enrichment (MAGMA tissue_spec + conditional retention)
    #######

    tx_tissue_data <- reactive({
      req(gwas_data(), selected_gwas())
      build_tissue_data(gwas_data(), selected_gwas())
    })

    tx_tissue_data_filtered <- reactive({
      d <- tx_tissue_data()
      if (is.null(d)) return(NULL)
      picked <- input$selected_tissues_tissue
      if (!is.null(picked) && length(picked) > 0) {
        d <- d[d$Tissue %in% picked, ]
      }
      if (isTRUE(as.logical(input$conf_only_tissue))) {
        d <- d[which(d$FDR_Sig), ]
      }
      d
    })

    output$tx_tissue_plot <- renderPlot({
      d <- tx_tissue_data_filtered()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      build_tissue_plot(
        d = d,
        n_total = nrow(tx_tissue_data()),
        sort_choice = input$sort_tissue %||% "significance",
        font_size = input$plot_font_size_tissue %||% 13,
        point_size = input$plot_point_size_tissue %||% 3,
        theme_fn = enr_theme_fn(input$plot_theme_tissue),
        title = input$plot_title_tissue %||% ""
      )
    })

    output$download_plot_tissue <- downloadHandler(
      filename = function() {
        sprintf("tissue_enrichment_%s.%s",
                format(Sys.time(), "%Y%m%d_%H%M%S"),
                input$dl_format_tissue %||% "png")
      },
      content = function(file) {
        p <- build_tissue_plot(
          d = tx_tissue_data_filtered(),
          n_total = nrow(tx_tissue_data() %||% data.frame()),
          sort_choice = input$sort_tissue %||% "significance",
          font_size = input$plot_font_size_tissue %||% 13,
          point_size = input$plot_point_size_tissue %||% 3,
          theme_fn = enr_theme_fn(input$plot_theme_tissue),
          title = input$plot_title_tissue %||% ""
        )
        if (is.null(p)) {
          grDevices::png(file, width = 4, height = 1, units = "in", res = 96)
          grid::grid.text("No tissues to plot.")
          grDevices::dev.off()
          return(invisible())
        }
        w <- input$dl_width_tissue  %||% 8
        h <- input$dl_height_tissue %||% 9
        fmt <- input$dl_format_tissue %||% "png"
        switch(fmt,
          png = grDevices::png(file, width = w, height = h, units = "in",
                               res = input$dl_dpi_tissue %||% 300),
          pdf = grDevices::pdf(file, width = w, height = h),
          svg = grDevices::svg(file, width = w, height = h)
        )
        print(p)
        grDevices::dev.off()
      }
    )

    output$tx_tissue_table <- renderDataTable({
      d <- tx_tissue_data_filtered()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )
      tmp <- d[, c("Tissue", "N Gene", "BETA", "SE", "P", "P.FDR", "Retained")]
      tmp$BETA <- round(tmp$BETA, 3)
      tmp$SE   <- round(tmp$SE, 3)
      datatable(
        tmp,
        rownames = FALSE,
        options = list(
          rowCallback = JS(js),
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 4:5)
          )
        ),
        escape = FALSE
      )
    })

  })
}
