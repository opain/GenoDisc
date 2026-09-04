#' Column-guide legend for a Molecular Associations results table.
#'
#' @param which One of "magma", "fusion", "smr_expr", "smr_protein", "panel".
#' @return A gd_legend() tag describing that table's columns.
mol_table_legend <- function(which) {
  lg_pfdr  <- "Benjamini-Hochberg FDR-adjusted p-value; < 0.05 is the usual significance threshold."
  lg_coloc <- paste0(
    "Colocalisation posterior probabilities: PP3 = expression and trait have distinct causal ",
    "variants; PP4 = they share one causal variant (high PP4 supports a genuine link).")
  items <- switch(which,
    magma = list(
      "CHR" = "Chromosome.",
      "START / STOP" = "Gene start and stop position (base pairs).",
      "ID" = "Gene symbol.",
      "P" = "MAGMA gene-based association p-value.",
      "P.FDR" = lg_pfdr),
    fusion = list(
      "PANEL" = "Expression/protein reference panel (tissue or dataset) used.",
      "CHR" = "Chromosome.",
      "P0 / P1" = "Gene start and end position (base pairs).",
      "Ensembl ID / Gene Symbol" = "Gene identifiers.",
      "Z" = "TWAS/PWAS Z-score: sign gives the direction and magnitude the strength of the molecular-trait association.",
      "P" = "TWAS/PWAS association p-value.",
      "P.FDR" = lg_pfdr,
      "COLOC.PP3 / COLOC.PP4" = lg_coloc,
      "High Confidence" = "TRUE = FDR-significant (P.FDR < 0.05) and colocalised (PP4-supported)."),
    smr_expr = list(
      "PANEL" = "Expression reference panel (eQTL dataset) used.",
      "CHR / BP" = "Chromosome and position of the probe/variant.",
      "Ensembl ID / Gene Symbol" = "Gene identifiers.",
      "BETA / SE" = "SMR effect estimate and its standard error.",
      "P" = "SMR association p-value.",
      "P.FDR" = lg_pfdr,
      "P (HEIDI)" = "HEIDI test p-value; > 0.05 favours a single shared causal variant over linkage.",
      "High Confidence" = "TRUE = FDR-significant SMR and HEIDI not rejected (P (HEIDI) > 0.05)."),
    smr_protein = list(
      "PANEL" = "Protein reference panel (pQTL dataset) used.",
      "CHR / BP" = "Chromosome and position of the probe/variant.",
      "Ensembl ID / Gene Symbol" = "Gene identifiers.",
      "BETA / SE" = "SMR effect estimate and its standard error.",
      "P" = "SMR association p-value.",
      "P.FDR" = lg_pfdr,
      "P (HEIDI)" = "HEIDI test p-value; > 0.05 favours a single shared causal variant over linkage."),
    panel = list(
      "Panel" = "Reference-panel name.",
      "Software" = "Method the panel is used with (FUSION or SMR).",
      "Type" = "Molecular feature type (expression, splicing or protein).",
      "N Individuals" = "Number of individuals in the panel's reference dataset.",
      "N Genes" = "Number of genes/features available in the panel.")
  )
  gd_legend(items, heading = "Column guide")
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

#' Resolve the "Plot options" theme dropdown value to a ggplot theme function.
mol_theme_fn <- function(name) {
  switch(name %||% "bw",
    bw = ggplot2::theme_bw,
    minimal = ggplot2::theme_minimal,
    classic = ggplot2::theme_classic,
    light = ggplot2::theme_light,
    ggplot2::theme_bw
  )
}

#' Build the Molecular Associations summary heatmap as a gtable
#'
#' Mirrors the drug/atc enrichment plot pattern (`build_tx_drug_gtable`):
#' returns a gtable that renderPlot draws via grid.draw, per-facet widths
#' scaled by a min-2 panel-count floor, and legend space reserved via
#' `reserve_legend_space` so the colourbar stays visible for short plots.
#'
#' @param all_func_res Filtered summary data (from `mol_assoc_summary_data_filtered()`).
#' @param locus_mode `TRUE` for the by-locus layout, `FALSE` for alphabetical.
#' @param locus_map Data frame with `ID`, `BP`, `Locus`, `locus_order` used
#'   in locus mode (ignored otherwise).
#' @param font_size Base font size in pt (default 14).
#' @param point_size Base point size for the main Z-score dots (default 4).
#' @param theme_fn A ggplot theme function (e.g. `ggplot2::theme_bw`).
#' @param title Optional plot title; empty string draws no title.
#' @param panel_h_pt Panel-row height (pt) from `plot_dim_mol()[["panel_h_pt"]]`;
#'   passed to `reserve_legend_space` so the colourbar stays visible on
#'   short plots.
#' @param left_pad_pt Left `plot.margin` in pt from
#'   `plot_dim_mol()[["left_pad_pt"]]` — absorbs rotated x-tick overflow of
#'   the leftmost facet so long panel names don't run off the screen.
build_mol_assoc_gtable <- function(all_func_res, locus_mode, locus_map,
                                    font_size = 14, point_size = 4,
                                    theme_fn = ggplot2::theme_bw, title = "",
                                    panel_h_pt = NULL, left_pad_pt = 0) {
  if (is.null(all_func_res) || nrow(all_func_res) == 0) return(NULL)

  # Fill missing Panel × Method × ID combos so every facet is present even
  # when the current gene subset has no data for it. Panel/Method slots
  # come from the ORIGINAL data (including Placeholder rows added upstream
  # in mol_assoc_summary_data_filtered), so a Panel/Method whose only
  # matching row is a Placeholder still shows up as an empty facet.
  dt <- data.table::as.data.table(all_func_res)
  if (nrow(dt) == 0) return(NULL)
  slots   <- unique(dt[, .(Panel, Method, Type)])
  all_ids <- unique(dt$ID[dt$ID != 'Placeholder'])
  if (length(all_ids) == 0) return(NULL)
  grid_dt <- slots[, .(ID = all_ids), by = .(Panel, Method, Type)]
  vals    <- dt[ID != 'Placeholder', .(Panel, Method, ID, Z, Sig, Coloc)]
  all_func_res_all <- vals[grid_dt, on = c("Panel", "Method", "ID")]
  all_func_res_all[is.na(Sig),   Sig   := 0]
  all_func_res_all[is.na(Coloc), Coloc := 0]

  all_func_res_all$Group <- make_group_labels(all_func_res_all$Method, all_func_res_all$Type)
  groups <- get_group_order()
  groups <- groups[groups %in% all_func_res_all$Group]
  all_func_res_all$Group <- factor(all_func_res_all$Group, levels = groups)

  # Locus mode requires positioned features; caller ensures locus_map is valid.
  if (locus_mode) {
    all_func_res_all <- merge(
      all_func_res_all, locus_map[, c("ID", "BP", "Locus", "locus_order")],
      by = "ID", all.x = TRUE)
    unplaced <- is.na(all_func_res_all$locus_order)
    max_order <- suppressWarnings(max(all_func_res_all$locus_order, na.rm = TRUE))
    if (!is.finite(max_order)) max_order <- 0
    all_func_res_all$Locus[unplaced] <- "Unplaced"
    all_func_res_all$locus_order[unplaced] <- max_order + 1
    all_func_res_all$BP[unplaced] <- Inf

    # NB `dt[i, "ID"]` on a data.table returns a data.table (not a vector),
    # so use `$ID[order]` to keep this compatible with the data.table pipeline.
    ord <- order(all_func_res_all$locus_order,
                 all_func_res_all$BP,
                 as.character(all_func_res_all$ID))
    gene_order <- unique(all_func_res_all$ID[ord])
    all_func_res_all$ID <- factor(all_func_res_all$ID, levels = rev(gene_order))
    locus_levels <- unique(all_func_res_all$Locus[order(all_func_res_all$locus_order)])
    all_func_res_all$Locus <- factor(all_func_res_all$Locus, levels = locus_levels)
  } else {
    all_func_res_all <- all_func_res_all[order(as.character(all_func_res_all$ID)), ]
    all_func_res_all$ID <- factor(all_func_res_all$ID, levels = rev(unique(as.character(all_func_res_all$ID))))
  }

  # Per-facet width scaling — mirrors build_tx_drug_gtable: min-2 floor on
  # panel count so 1-panel facets (MAGMA / SuSiE / Nearest Gene) keep enough
  # room for their multi-line strip titles.
  group_siz <- do.call(rbind, lapply(groups, function(g)
    data.frame(Group = g, Size = length(unique(all_func_res_all$Panel[all_func_res_all$Group == g])))))
  group_siz$Size[group_siz$Size < 2] <- 2
  group_siz$Prop  <- group_siz$Size / sum(group_siz$Size)
  group_siz$Width <- 4 * group_siz$Prop

  x <- c(-max(abs(all_func_res_all$Z), na.rm = TRUE), 0, max(abs(all_func_res_all$Z), na.rm = TRUE))
  x <- (x - min(x)) / (max(x) - min(x))

  all_func_res_all <- data.table::data.table(all_func_res_all)

  heatmap <- ggplot2::ggplot(data = all_func_res_all, ggplot2::aes(x = Panel, y = ID)) +
    theme_fn() +
    # geom_blank on the full data ensures every facet's x-scale is trained,
    # so empty facets (Panel/Method combos with no matching gene data) still
    # show their panel background, x-tick labels, and gridlines.
    ggplot2::geom_blank() +
    ggplot2::geom_point(data = all_func_res_all[all_func_res_all$Sig == TRUE, ], ggplot2::aes(x = Panel, y = ID), colour = 'black', size = point_size + 1) +
    ggplot2::geom_point(data = all_func_res_all[all_func_res_all$Coloc == TRUE & all_func_res_all$Sig == TRUE, ], ggplot2::aes(x = Panel, y = ID), colour = 'black', shape = 15, size = point_size + 2) +
    ggplot2::geom_point(data = all_func_res_all[(all_func_res_all$Method == 'SNP\nFine-mapping' | all_func_res_all$Method == 'Nearest\nGene') & !is.na(all_func_res_all$Z), ], ggplot2::aes(x = Panel, y = ID), colour = '#00FF00', size = point_size + 1) +
    ggplot2::geom_point(data = all_func_res_all[all_func_res_all$Method != 'SNP\nFine-mapping' & all_func_res_all$Method != 'Nearest\nGene', ], ggplot2::aes(colour = Z), size = point_size) +
    ggplot2::scale_colour_gradientn(colours = c("#0066FF", "#0099FF", "#FFFFFF", "#FF6666", "#FF0000"), na.value = NA, name = "Z-score", limits = c(-max(abs(all_func_res_all$Z), na.rm = TRUE), max(abs(all_func_res_all$Z), na.rm = TRUE)), values = x) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(x = '', y = '', title = if (nzchar(title)) title else NULL) +
    ggplot2::theme(text = ggplot2::element_text(size = font_size),
                   plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5,
                                                 l = left_pad_pt, unit = "pt"))

  if (locus_mode) {
    # `space = "free_y"` gives each locus row a height proportional to its
    # gene count, while leaving x-widths at 1null each — that lets the
    # `group_siz$Width` multiplier below apply proportional widths (with
    # the min-2 floor) exactly as it does in gene mode. Using `space =
    # "free"` here would pre-scale x-widths by arity and the multiplier
    # would double-scale.
    heatmap <- heatmap +
      ggplot2::facet_grid(Locus ~ Group, scales = "free", space = "free_y") +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0),
                     panel.spacing.y = grid::unit(4, "pt"))
  } else {
    heatmap <- heatmap +
      ggplot2::facet_wrap(~ Group, nrow = 1, scales = "free_x") +
      ggplot2::scale_y_discrete(limits = unique(rev(all_func_res_all$ID)))
  }

  gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(heatmap))
  # Multiply per-facet width by group_siz$Width so 1-panel facets get the
  # min-2 floor. Panel naming is `panel-<x_facet_index>-<y_facet_index>`
  # for both facet_wrap (y_index always 1) and facet_grid, so the same
  # regex targets column i in both modes.
  for (i in seq_len(nrow(group_siz))) {
    panel_l <- unique(gt$layout$l[grep(paste0('^panel-', i, '-\\d+$'), gt$layout$name)])
    if (length(panel_l) > 0) {
      gt$widths[panel_l] <- group_siz$Width[i] * gt$widths[panel_l]
    }
  }
  reserve_legend_space(gt, panel_h_pt)
}

#' Molecular Associations module UI
#'
#' @param id Module namespace id
#' @return UI elements for the Molecular Associations tab
molAssocUI <- function(id) {
  ns <- NS(id)

  tabPanel(
    title = "Molecular Associations",
    br(),
    p("This tab shows molecular association results. Select the tabs below to see a summary of results across methods, or method-specific results tables."),
    hr(),
    tabsetPanel(
      id = ns("mol_assoc_tabset"),
      tabPanel(
        title = "Summary",
        br(),
        # Body swaps between the single-GWAS UI (rendered server-side) and
        # the cross-GWAS gene compare UI when >= 2 GWAS are selected.
        uiOutput(ns("summary_body"))
      ),
      tabPanel(title = "MAGMA", br(),
        p("This tab shows MAGMA gene association results."), hr(), br(),
        fluidRow(column(width = 6, dataTableOutput(ns("mol_assoc_magma_table")))),
        mol_table_legend("magma"), br()
      ),
      tabPanel(title = "Expression - FUSION", br(),
        p("This tab shows differential expression association results from FUSION."), hr(), br(),
        fluidRow(column(width = 9, dataTableOutput(ns("mol_assoc_fusion_expr_table")))),
        mol_table_legend("fusion"), br()
      ),
      tabPanel(title = "Protein - FUSION", br(),
        p("This tab shows differential protein level association results from FUSION."), hr(), br(),
        fluidRow(column(width = 9, dataTableOutput(ns("mol_assoc_fusion_protein_table")))),
        mol_table_legend("fusion"), br()
      ),
      tabPanel(title = "Expression - SMR", br(),
        p("This tab shows differential expression association results from SMR"), hr(), br(),
        fluidRow(column(width = 9, dataTableOutput(ns("mol_assoc_smr_expr_table")))),
        mol_table_legend("smr_expr"), br()
      ),
      tabPanel(title = "Protein - SMR", br(),
        p("This tab shows differential protein level association results from SMR"), hr(), br(),
        fluidRow(column(width = 9, dataTableOutput(ns("mol_assoc_smr_protein_table")))),
        mol_table_legend("smr_protein"), br()
      ),
      tabPanel(title = "Panel Info.", br(),
        p("This tab shows the number of features and individuals for each panel."), hr(), br(),
        fluidRow(column(width = 7, dataTableOutput(ns("panel_info_table")))),
        mol_table_legend("panel"), br()
      )
    )
  )
}

#' Molecular Associations module server
#'
#' @param id Module namespace id
#' @param gwas_data Reactive returning the loaded GWAS data list
#' @param selected_gwas Reactive returning the currently selected GWAS name
#' @param config_flags Reactive returning a list from parse_config_flags()
molAssocServer <- function(id, gwas_data, selected_gwas, config_flags,
                            selected_gwas_multi = NULL,
                            comparison_mode = NULL,
                            comparison_long = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Cross-GWAS gene compare sub-module. Registered unconditionally so its
    # outputs exist even when we're in single-GWAS mode; it self-guards via
    # its own reactives.
    if (!is.null(selected_gwas_multi) && !is.null(comparison_long)) {
      gene_compare_server("gene_compare",
                            gwas_data, selected_gwas_multi, comparison_long)
    }

    # Summary tab body: single-GWAS UI unless >= 2 GWAS are selected, in
    # which case swap to the cross-GWAS gene compare UI. The rest of the
    # method-specific tabs (MAGMA, FUSION, SMR, etc.) continue to show
    # single-GWAS content for the first-selected GWAS.
    output$summary_body <- renderUI({
      in_compare <- !is.null(comparison_mode) && isTRUE(comparison_mode())
      if (in_compare) {
        return(gene_compare_ui(NS(ns("gene_compare"))))
      }
      tagList(
        p("This tab shows a heatmap summarising, for each gene, whether it is associated with the trait across every method and reference panel included in the analysis. Use ", tags$b("Filter data"), " to control which results appear and ", tags$b("Plot options"), " to customise or download the figure."),
        hr(),
        tags$details(class = "gd-details",
          tags$summary("Filter data"),
          tags$div(class = "gd-details-body",
            tags$p(class = "gd-details-intro",
              "Choose which molecular-association results appear in the heatmap, ",
              "restrict the view to specific genes, and decide how a gene is ",
              "labelled ", tags$em("high-confidence"), "."),
            fluidRow(
              column(4,
                selectInput(ns("selected_methods_mol"), "Include results from these methods:", "", multiple = T),
                selectInput(ns("selected_expr_panels_mol"), "Include expression / splicing panels:", "", multiple = T),
                selectInput(ns("selected_protein_panels_mol"), "Include protein panels:", "", multiple = T)
              ),
              column(4,
                textInput(ns("geneInput_mol"), "Show only these genes (comma- or space-separated):"),
                selectInput(ns("selected_group_hc_mol"),
                            "Define 'high-confidence' as significant in any of:",
                            "", multiple = T)
              ),
              column(4,
                radioButtons(ns("mol_layout"), "Feature ordering:",
                             choices = c("Alphabetical" = "alphabetical",
                                         "By locus (genomic position)" = "locus"),
                             selected = "alphabetical"),
                radioButtons(ns("conf_only_mol"), "Show high-confidence genes only :",
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
                textInput(ns("plot_title"), "Plot title (optional):", value = ""),
                selectInput(ns("plot_theme"), "Theme:",
                            choices = c("Black & white" = "bw",
                                        "Minimal" = "minimal",
                                        "Classic" = "classic",
                                        "Light" = "light"),
                            selected = "bw")
              ),
              column(4,
                sliderInput(ns("plot_font_size"), "Font size (pt):",
                            min = 10, max = 20, value = 14, step = 1),
                sliderInput(ns("plot_point_size"), "Point size:",
                            min = 2, max = 8, value = 4, step = 1)
              ),
              column(4,
                selectInput(ns("dl_format"), "Download format:",
                            choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                            selected = "png"),
                numericInput(ns("dl_width"), "Width (inches):",
                             value = 12, min = 2, max = 40, step = 0.5),
                numericInput(ns("dl_height"), "Height (inches):",
                             value = 8, min = 2, max = 40, step = 0.5),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'png'", ns("dl_format")),
                  numericInput(ns("dl_dpi"), "Resolution (DPI, PNG only):",
                               value = 300, min = 72, max = 600, step = 25)
                ),
                downloadButton(ns("download_plot"), "Download plot")
              )
            )
          )
        ),
        br(),
        uiOutput(ns("message_too_large_mol")),
        uiOutput(ns("message_no_genes_mol")),
        uiOutput(ns("message_no_positions_mol")),
        uiOutput(ns("mol_assoc_plot.ui"))
      )
    })

    ########
    # Show/hide tabs based on config
    ########

    # Per-method sub-tabs that duplicate what the compare Summary tab
    # already shows via its Method dropdown. Hidden in compare mode; shown
    # (subject to config-flag gates below) in single-GWAS mode.
    .redundant_method_tabs <- c("MAGMA", "Expression - FUSION",
                                  "Protein - FUSION", "Expression - SMR",
                                  "Protein - SMR")

    apply_mol_assoc_visibility <- function() {
      cf <- config_flags()
      if (is.null(cf)) return()
      tab_id <- "mol_assoc_tabset"
      in_compare <- !is.null(comparison_mode) && isTRUE(comparison_mode())

      if (in_compare) {
        for (t in .redundant_method_tabs) hideTab(tab_id, t)
        # If the user was viewing one of the hidden method sub-tabs,
        # jump to the compare Summary.
        cur <- isolate(input$mol_assoc_tabset)
        if (!is.null(cur) && cur %in% .redundant_method_tabs) {
          updateTabsetPanel(session, tab_id, selected = "Summary")
        }
        return()
      }

      if (cf$magma_gene) showTab(tab_id, "MAGMA") else hideTab(tab_id, "MAGMA")
      if (cf$twas) showTab(tab_id, "Expression - FUSION") else hideTab(tab_id, "Expression - FUSION")
      if (any(cf$pwas_panel_rosmap, cf$pwas_panel_banner))
        showTab(tab_id, "Protein - FUSION") else hideTab(tab_id, "Protein - FUSION")
      if (cf$smr_expression) showTab(tab_id, "Expression - SMR") else hideTab(tab_id, "Expression - SMR")
      if (cf$smr_protein_panel_rosmap) showTab(tab_id, "Protein - SMR") else hideTab(tab_id, "Protein - SMR")
    }

    observeEvent(config_flags(), apply_mol_assoc_visibility())
    if (!is.null(comparison_mode)) {
      observeEvent(comparison_mode(), apply_mol_assoc_visibility())
    }

    ########
    # Update selectInputs when summary data changes
    ########

    observe({
      all_func_res <- mol_assoc_summary_data()
      req(all_func_res)

      methods <- unique(all_func_res$Method)
      updateSelectInput(session, "selected_methods_mol", choices = methods, selected = methods)

      expr_panels <- unique(all_func_res$Panel[all_func_res$Type == 'Expr.' | all_func_res$Type == 'Splice'])
      updateSelectInput(session, "selected_expr_panels_mol", choices = expr_panels, selected = expr_panels)

      protein_panels <- unique(all_func_res$Panel[all_func_res$Type == 'Protein'])
      updateSelectInput(session, "selected_protein_panels_mol", choices = protein_panels, selected = protein_panels)

      res_group <- paste0(all_func_res$Method, '\n', all_func_res$Type)
      res_group[res_group == 'SNP\nFine-mapping\n'] <- 'SuSiE'
      res_group[res_group == 'MAGMA\n'] <- 'MAGMA'
      res_group[res_group == 'Nearest\nGene\n'] <- 'Nearest\nGene'

      hc_groups <- c('SuSiE', 'FUSION\nExpr.', 'FUSION\nSplice', 'SMR\nExpr.', 'FUSION\nProtein', 'SMR\nProtein')
      hc_groups <- hc_groups[hc_groups %in% res_group]

      updateSelectInput(session, "selected_group_hc_mol", choices = hc_groups, selected = hc_groups)
    })

    ########
    # Individual method tables
    ########

    output$mol_assoc_magma_table <- renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )

      tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/magma")

      datatable(tmp, rownames = F, options = list(
        rowCallback = JS(js),
        columnDefs = list(list(className = 'dt-center', targets = 0:5))))
    })

    output$mol_assoc_fusion_expr_table <- renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "}"
      )

      tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/exp/fusion")$res
      tmp$TWAS.Z <- round(tmp$TWAS.Z, 3)
      tmp$`High Confidence` <- tmp$TWAS.P.FDR < 0.05 & tmp$COLOC_logical
      tmp$COLOC_logical <- NULL
      tmp[is.na(tmp)] <- 'NA'

      names(tmp)[names(tmp) == 'TWAS.Z'] <- 'Z'
      names(tmp)[names(tmp) == 'TWAS.P'] <- 'P'
      names(tmp)[names(tmp) == 'TWAS.P.FDR'] <- 'P.FDR'

      tmp <- tmp[, c("PANEL","CHR","P0","P1","Ensembl ID","Gene Symbol","Z","P","P.FDR","COLOC.PP3","COLOC.PP4","High Confidence"), with = F]

      datatable(tmp, rownames = F, options = list(
        rowCallback = JS(js),
        columnDefs = list(list(className = 'dt-center', targets = 0:11))))
    })

    output$mol_assoc_fusion_protein_table <- renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "}"
      )

      tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/protein/fusion")$res
      tmp$pwas_all.Z <- round(tmp$pwas_all.Z, 3)
      tmp$`High Confidence` <- tmp$pwas_all.P.FDR < 0.05 & tmp$COLOC_logical
      tmp$COLOC_logical <- NULL
      tmp[is.na(tmp)] <- 'NA'

      names(tmp)[names(tmp) == 'pwas_all.Z'] <- 'Z'
      names(tmp)[names(tmp) == 'pwas_all.P'] <- 'P'
      names(tmp)[names(tmp) == 'pwas_all.P.FDR'] <- 'P.FDR'

      tmp <- tmp[, c("PANEL","CHR","P0","P1","Ensembl ID","Gene Symbol","Z","P","P.FDR","COLOC.PP3","COLOC.PP4","High Confidence"), with = F]

      datatable(tmp, rownames = F, options = list(
        rowCallback = JS(js),
        columnDefs = list(list(className = 'dt-center', targets = 0:11))))
    })

    output$mol_assoc_smr_expr_table <- renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "  var x = data[9];",
        "  $('td:eq(9)', row).html(x.toExponential(2));",
        "}"
      )

      tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/exp/smr")$res
      tmp$`High Confidence` <- tmp$p_SMR.FDR < 0.05 & tmp$p_HEIDI > 0.05
      tmp <- tmp[, c("PANEL","CHR","BP","Ensembl ID","Gene Symbol","b_SMR","se_SMR","p_SMR","p_SMR.FDR","p_HEIDI","High Confidence"), with = F]
      tmp$b_SMR <- round(tmp$b_SMR, 3)
      tmp$se_SMR <- round(tmp$se_SMR, 3)

      names(tmp)[names(tmp) == 'b_SMR'] <- 'BETA'
      names(tmp)[names(tmp) == 'se_SMR'] <- 'SE'
      names(tmp)[names(tmp) == 'p_SMR'] <- 'P'
      names(tmp)[names(tmp) == 'p_SMR.FDR'] <- 'P.FDR'
      names(tmp)[names(tmp) == 'p_HEIDI'] <- "P (HEIDI)"

      datatable(tmp, rownames = F, options = list(
        rowCallback = JS(js),
        columnDefs = list(list(className = 'dt-center', targets = 0:10))))
    })

    output$mol_assoc_smr_protein_table <- renderDataTable({
      req(gwas_data(), selected_gwas())
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "  var y = data[9];",
        "  $('td:eq(9)', row).html(y.toExponential(2));",
        "}"
      )

      tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/protein/smr")$res
      tmp$`High Confidence` <- tmp$p_SMR.FDR < 0.05 & tmp$p_HEIDI > 0.05
      tmp <- tmp[, c("PANEL","CHR","BP","Ensembl ID","Gene Symbol","b_SMR","se_SMR","p_SMR","p_SMR.FDR","p_HEIDI"), with = F]
      tmp$b_SMR <- round(tmp$b_SMR, 3)
      tmp$se_SMR <- round(tmp$se_SMR, 3)

      names(tmp)[names(tmp) == 'b_SMR'] <- 'BETA'
      names(tmp)[names(tmp) == 'se_SMR'] <- 'SE'
      names(tmp)[names(tmp) == 'p_SMR'] <- 'P'
      names(tmp)[names(tmp) == 'p_SMR.FDR'] <- 'P.FDR'
      names(tmp)[names(tmp) == 'p_HEIDI'] <- "P (HEIDI)"

      datatable(tmp, rownames = F, options = list(
        rowCallback = JS(js),
        columnDefs = list(list(className = 'dt-center', targets = 0:9))))
    })

    ########
    # Panel info table
    ########

    output$panel_info_table <- renderDataTable({
      req(gwas_data(), selected_gwas(), config_flags())
      cf <- config_flags()

      all_panel_info <- NULL

      if (cf$twas) {
        tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/exp/fusion")$panels
        all_panel_info <- rbind(all_panel_info, tmp)
      }

      if (cf$smr_expression) {
        tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/exp/smr")$panels
        all_panel_info <- rbind(all_panel_info, tmp)
      }

      if (any(cf$pwas_panel_rosmap, cf$pwas_panel_banner)) {
        tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/protein/fusion")$panels
        all_panel_info <- rbind(all_panel_info, tmp)
      }

      if (cf$smr_protein_panel_rosmap) {
        tmp <- gd_read(gwas_data(), selected_gwas(), "mol_assoc/protein/smr")$panels
        all_panel_info <- rbind(all_panel_info, tmp)
      }

      all_panel_info <- all_panel_info[order(all_panel_info$Software, all_panel_info$Type, all_panel_info$Panel), ]

      names(all_panel_info)[names(all_panel_info) == 'N_indiv'] <- 'N Individuals'
      names(all_panel_info)[names(all_panel_info) == 'N_gene'] <- 'N Genes'

      datatable(all_panel_info, rownames = F, options = list(
        columnDefs = list(list(className = 'dt-center', targets = '_all'))))
    })

    ########
    # Summary data reactive
    ########

    mol_assoc_summary_data <- reactive({
      req(gwas_data(), selected_gwas(), config_flags())
      build_mol_assoc_data(gwas_data(), selected_gwas(), config_flags())
    })

    # Gene -> locus lookup for the per-locus layout. Positions are resolved
    # through the canonical reference (with a per-method harvest fallback) and
    # scoped to the DISPLAYED genes only, so loci are tight clusters of the
    # shown features rather than genome-wide bands.
    mol_locus_map <- reactive({
      req(gwas_data(), selected_gwas(), config_flags(), mol_assoc_summary_data_filtered())
      ids <- unique(as.character(mol_assoc_summary_data_filtered()$ID))
      ids <- ids[ids != "Placeholder"]
      pos <- resolve_feature_positions(ids, gwas_data(), selected_gwas(), config_flags())
      assign_loci(pos)
    })

    ########
    # Filtered summary data with debounce
    ########

    mol_assoc_summary_data_filtered <- reactive({
      req(mol_assoc_summary_data(), input$selected_methods_mol)
      all_func_res <- mol_assoc_summary_data()

      # Determine HC genes from the FULL, unfiltered data using the
      # "high-confidence" method selection. This way a gene identified by
      # SuSiE fine-mapping still counts as HC even if SNP fine-mapping is
      # not in the "Include results from these methods:" selection.
      hc_genes <- NULL
      if (input$conf_only_mol) {
        res_group <- paste0(all_func_res$Method, '\n', all_func_res$Type)
        res_group[res_group == 'SNP\nFine-mapping\n'] <- 'SuSiE'
        res_group[res_group == 'MAGMA\n']             <- 'MAGMA'
        res_group[res_group == 'Nearest\nGene\n']     <- 'Nearest\nGene'
        hc_genes <- all_func_res$ID[which(
          ((all_func_res$Sig == T & all_func_res$Coloc == T) |
            all_func_res$Panel == "SuSie (L=1)") &
          res_group %in% input$selected_group_hc_mol)]
      }

      # Filter results table by user specified methods (for DISPLAY only —
      # HC status was already determined above from the full data).
      all_func_res <- all_func_res[all_func_res$Method %in% input$selected_methods_mol, ]

      # Filter results table by user specified expression and protein panels
      if (any(all_func_res$Type == 'Expr.' | all_func_res$Type == 'Splice')) {
        all_func_res <- all_func_res[!((all_func_res$Type == 'Expr.' | all_func_res$Type == 'Splice') & !(all_func_res$Panel %in% input$selected_expr_panels_mol)), ]
      }
      if (any(all_func_res$Type == 'Protein')) {
        all_func_res <- all_func_res[!(all_func_res$Type == 'Protein' & !(all_func_res$Panel %in% input$selected_protein_panels_mol)), ]
      }

      # Insert NA rows for all panels and methods so when filtering by gene, all selected panels and methods remain
      na_rows <- all_func_res[!(duplicated(paste0(all_func_res$Panel, all_func_res$Method))), ]
      na_rows$ID <- 'Placeholder'
      na_rows$Z <- NA
      na_rows$Sig <- NA
      na_rows$Coloc <- NA

      all_func_res <- rbind(na_rows, all_func_res)

      # Restrict rows to HC genes (if the toggle was on).
      if (input$conf_only_mol) {
        all_func_res <- all_func_res[all_func_res$ID %in% c(hc_genes, 'Placeholder'), ]
      }

      input_genes <- unlist(strsplit(input$geneInput_mol, "[, ]"))
      selected_genes <- input_genes[input_genes != ""]

      if (length(selected_genes) > 0) {
        if (sum(grepl(paste(selected_genes, collapse = '|'), all_func_res$ID, ignore.case = T)) > 0) {
          selected_genes <- c(selected_genes, 'Placeholder')
          all_func_res <- all_func_res[grepl(paste(selected_genes, collapse = '|'), all_func_res$ID, ignore.case = T) & !is.na(all_func_res$ID), ]
        } else {
          all_func_res <- data.frame(matrix(nrow = 0, ncol = 5))
        }
      }

      return(all_func_res)
    }) %>% debounce(1000)

    ########
    # Plot dimensions
    ########

    plot_dim_mol <- reactive({
      all_func_res <- mol_assoc_summary_data_filtered()
      # Filtering to a gene with no HC data returns a 0-row placeholder
      # frame with no Method/Type columns; short-circuit before touching them.
      if (is.null(all_func_res) || nrow(all_func_res) == 0) {
        return(list(height = 100, width = 100, panel_h_pt = 0, left_pad_pt = 0))
      }
      all_func_res$Group <- make_group_labels(all_func_res$Method, all_func_res$Type)
      # Similar to plot_dim_drug, but with smaller mins: mol_assoc has a
      # single Z-score colourbar (drug stacks two), so it needs less
      # reserved panel height for the legend.
      dims <- calc_plot_dims(all_func_res, y_col = "ID", x_col = "Panel", facet_col = "Group",
                             facet_order = get_group_order(),
                             font_size = input$plot_font_size %||% 14,
                             min_height = 200, min_panel_h_pt = 120)

      # Locus mode adds a right-hand Locus strip plus panel spacing between
      # locus bands (4pt gaps).
      if (identical(input$mol_layout, "locus")) {
        genes <- unique(all_func_res$ID[all_func_res$ID != 'Placeholder'])
        lmap <- mol_locus_map()
        n_loci <- length(unique(lmap$Locus[lmap$ID %in% genes]))
        if (any(!(genes %in% lmap$ID))) n_loci <- n_loci + 1
        dims$height <- dims$height + 4 * max(n_loci - 1, 0)
        dims$width  <- dims$width + 120
      }
      dims
    })

    ########
    # Heatmap plot
    ########

    output$mol_assoc_plot <- renderPlot({
      all_func_res <- mol_assoc_summary_data_filtered()
      if (plot_dim_mol()[['height']] >= 10000 || nrow(all_func_res) == 0) return(NULL)
      gt <- build_mol_assoc_gtable(
        all_func_res = all_func_res,
        locus_mode = identical(input$mol_layout, "locus") && nrow(mol_locus_map()) > 0,
        locus_map = mol_locus_map(),
        font_size = input$plot_font_size %||% 14,
        point_size = input$plot_point_size %||% 4,
        theme_fn = mol_theme_fn(input$plot_theme),
        title = input$plot_title %||% "",
        panel_h_pt = plot_dim_mol()[['panel_h_pt']],
        left_pad_pt = plot_dim_mol()[['left_pad_pt']]
      )
      if (!is.null(gt)) grid::grid.draw(gt)
    })

    ########
    # UI outputs for plot and messages
    ########

    output$mol_assoc_plot.ui <- renderUI({
      ns <- session$ns
      filtered <- mol_assoc_summary_data_filtered()
      has_real_genes <- any(filtered$ID != 'Placeholder')
      if (plot_dim_mol()[['height']] < 10000 && nrow(filtered) > 0 && has_real_genes) {
        legend_items <- list(
          "Rows / columns" = "Each row is a gene; each column is an expression or protein reference panel, grouped by method (shown in the facet headers).",
          "Colour" = "Association Z-score: blue = the molecular feature is lower with the trait-increasing allele, red = higher, white is approximately no association.",
          "Black outline" = "The association is FDR-significant in that panel and method.",
          "Black square" = "FDR-significant AND supported by colocalisation.",
          "Green" = "Gene highlighted by SuSiE fine-mapping or as the nearest gene to a lead variant."
        )
        if (identical(input$mol_layout, "locus")) {
          legend_items <- c(list(
            "Locus bands" = "Genes are grouped into loci (features within 500 kb, labelled by coordinate range) and ordered by genomic position; genes with no known position appear in a trailing 'Unplaced' band."
          ), legend_items)
        }
        tagList(
          plotOutput(ns("mol_assoc_plot"), height = plot_dim_mol()[['height']], width = plot_dim_mol()[['width']]),
          gd_legend(legend_items, heading = "How to read this plot")
        )
      } else {
        NULL
      }
    })

    output$message_too_large_mol <- renderUI({
      if (plot_dim_mol()[['height']] > 10000) {
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "Plot is too large. Restrict to high-confidence genes or specify a list of genes."
        ))
      }
    })

    output$message_no_positions_mol <- renderUI({
      req(config_flags())
      if (identical(input$mol_layout, "locus") && nrow(mol_locus_map()) == 0) {
        HTML(paste0(
          "<div style='color: #666; font-style: italic; margin: 10px 0;'>",
          "No genomic positions are available for these results (position-bearing ",
          "methods such as MAGMA, FUSION or SMR were not run), so the plot is shown ",
          "alphabetically instead of by locus.",
          "</div>"
        ))
      }
    })

    output$message_no_genes_mol <- renderUI({
      filtered <- mol_assoc_summary_data_filtered()
      real_genes <- filtered$ID[filtered$ID != 'Placeholder']
      if (length(real_genes) == 0 && as.logical(input$conf_only_mol)) {
        HTML(paste0(
          "<div style='color: #666; font-style: italic; margin: 10px 0;'>",
          "No high-confidence genes were identified. ",
          "Set <b>Show high-confidence only</b> to <b>False</b>, ",
          "and specify gene symbols using <b>Enter gene symbols</b> to visualise results for those genes.",
          "</div>"
        ))
      } else if (nrow(filtered) == 0) {
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "No genes are present."
        ))
      }
    })

    ########
    # Plot download: keep width/height in sync with the on-screen plot so
    # the default download matches what the user is currently looking at.
    # User overrides get replaced whenever the on-screen dims change.
    ########

    observeEvent(plot_dim_mol(), {
      dims <- plot_dim_mol()
      if (dims$width > 100 && dims$height < 10000) {
        updateNumericInput(session, "dl_width",  value = round(dims$width  / 72, 1))
        updateNumericInput(session, "dl_height", value = round(dims$height / 72, 1))
      }
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        sprintf("mol_assoc_summary_%s.%s",
                format(Sys.time(), "%Y%m%d_%H%M%S"),
                input$dl_format %||% "png")
      },
      content = function(file) {
        all_func_res <- mol_assoc_summary_data_filtered()
        if (is.null(all_func_res) || nrow(all_func_res) == 0) {
          # Write a small placeholder so the file isn't empty on error.
          grDevices::png(file, width = 4, height = 1, units = "in", res = 96)
          grid::grid.text("No data to plot.")
          grDevices::dev.off()
          return(invisible())
        }
        gt <- build_mol_assoc_gtable(
          all_func_res = all_func_res,
          locus_mode = identical(input$mol_layout, "locus") && nrow(mol_locus_map()) > 0,
          locus_map = mol_locus_map(),
          font_size = input$plot_font_size %||% 14,
          point_size = input$plot_point_size %||% 4,
          theme_fn = mol_theme_fn(input$plot_theme),
          title = input$plot_title %||% "",
          panel_h_pt = plot_dim_mol()[['panel_h_pt']],
          left_pad_pt = plot_dim_mol()[['left_pad_pt']]
        )
        w <- input$dl_width  %||% 12
        h <- input$dl_height %||% 8
        fmt <- input$dl_format %||% "png"
        switch(fmt,
          png = grDevices::png(file, width = w, height = h, units = "in",
                               res = input$dl_dpi %||% 300),
          pdf = grDevices::pdf(file, width = w, height = h),
          svg = grDevices::svg(file, width = w, height = h)
        )
        if (!is.null(gt)) grid::grid.draw(gt)
        grDevices::dev.off()
      }
    )
  })
}
