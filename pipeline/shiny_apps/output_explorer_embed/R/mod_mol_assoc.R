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
        fluidPage(
          sidebarPanel(
            selectInput(ns("selected_methods_mol"), "Select methods", "", multiple = T),
            selectInput(ns("selected_expr_panels_mol"), "Select expression panels", "", multiple = T),
            selectInput(ns("selected_protein_panels_mol"), "Select protein panels", "", multiple = T),
            radioButtons(ns("conf_only_mol"), "Show high-confidence only :",
                         choices = c("True" = T, "False" = F), selected = T),
            radioButtons(ns("mol_layout"), "Feature ordering:",
                         choices = c("Alphabetical" = "alphabetical",
                                     "By locus (genomic position)" = "locus"),
                         selected = "alphabetical"),
            textInput(ns("geneInput_mol"), "Enter gene symbols (whitespace- or comma-seperated):"),
            hr(),
            h5('Select high confidence criteria:'),
            selectInput(ns("selected_group_hc_mol"), "Select methods", "", multiple = T)
          ),
          mainPanel(
            uiOutput(ns("message_too_large_mol")),
            uiOutput(ns("message_no_genes_mol")),
            uiOutput(ns("message_no_positions_mol")),
            uiOutput(ns("mol_assoc_plot.ui"))
          )
        )
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
molAssocServer <- function(id, gwas_data, selected_gwas, config_flags) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ########
    # Show/hide tabs based on config
    ########

    observeEvent(config_flags(), {
      cf <- config_flags()
      tab_id <- "mol_assoc_tabset"
      if (!cf$magma_gene) hideTab(tab_id, "MAGMA")
      if (!cf$twas) hideTab(tab_id, "Expression - FUSION")
      if (!any(cf$pwas_panel_rosmap, cf$pwas_panel_banner)) hideTab(tab_id, "Protein - FUSION")
      if (!cf$smr_expression) hideTab(tab_id, "Expression - SMR")
      if (!cf$smr_protein_panel_rosmap) hideTab(tab_id, "Protein - SMR")
    })

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

      # Filter results table by user specified methods
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

      # Filter results table if user specifies high confidence genes only
      if (input$conf_only_mol) {
        # Create group variable
        res_group <- paste0(all_func_res$Method, '\n', all_func_res$Type)
        res_group[res_group == 'SNP\nFine-mapping\n'] <- 'SuSiE'
        res_group[res_group == 'MAGMA\n'] <- 'MAGMA'
        res_group[res_group == 'Nearest\nGene\n'] <- 'Nearest\nGene'

        hc_genes <- all_func_res$ID[which(((all_func_res$Sig == T & all_func_res$Coloc == T) | all_func_res$Panel == "SuSie (L=1)") & res_group %in% input$selected_group_hc_mol)]
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
      all_func_res$Group <- make_group_labels(all_func_res$Method, all_func_res$Type)
      dims <- calc_plot_dims(all_func_res, y_col = "ID", x_col = "Panel", facet_col = "Group",
                             facet_order = get_group_order())

      # Per-locus layout adds a horizontal band (facet row) per locus plus a
      # right-hand locus strip, so allow extra height/width for that chrome.
      if (identical(input$mol_layout, "locus")) {
        genes <- unique(all_func_res$ID[all_func_res$ID != 'Placeholder'])
        lmap <- mol_locus_map()
        n_loci <- length(unique(lmap$Locus[lmap$ID %in% genes]))
        if (any(!(genes %in% lmap$ID))) n_loci <- n_loci + 1  # trailing "Unplaced" band
        dims$height <- dims$height + 30 * max(n_loci, 1)
        dims$width  <- dims$width + 120
      }
      dims
    })

    ########
    # Heatmap plot
    ########

    output$mol_assoc_plot <- renderPlot({

      all_func_res <- mol_assoc_summary_data_filtered()

      if (plot_dim_mol()[['height']] < 10000 & nrow(all_func_res) > 0) {

        # Insert missing data
        all_func_res_all <- NULL
        for (i in unique(all_func_res$Panel)) {
          for (j in unique(all_func_res$Method[all_func_res$Panel == i])) {

            all_func_res_panel <- all_func_res[all_func_res$Panel == i & all_func_res$Method == j, ]
            all_func_res_other <- all_func_res[!(all_func_res$Panel %in% all_func_res_panel$Panel) & !(all_func_res$Method %in% all_func_res_panel$Method), ]
            all_func_res_other <- all_func_res_other[!(all_func_res_other$ID %in% all_func_res_panel$ID), ]
            all_func_res_other <- unique(all_func_res_other$ID)

            if (length(all_func_res_other) > 0) {
              all_func_res_panel_missing <- data.frame(ID = all_func_res_other)

              all_func_res_panel_missing$Panel <- i
              all_func_res_panel_missing$ID <- all_func_res_other
              all_func_res_panel_missing$Z <- NA
              all_func_res_panel_missing$Sig <- 0
              all_func_res_panel_missing$Coloc <- 0
              all_func_res_panel_missing$Method <- j
              all_func_res_panel_missing$Type <- all_func_res_panel$Type[1]
              all_func_res_panel_missing$Group <- all_func_res_panel$Group[1]

              all_func_res_panel_missing <- all_func_res_panel_missing[, names(all_func_res_panel)]

              all_func_res_all <- rbind(all_func_res_all, all_func_res_panel_missing)
            }

            all_func_res_all <- rbind(all_func_res_all, all_func_res_panel)
          }
        }

        # Now remove the NA rows
        all_func_res_all <- all_func_res_all[all_func_res_all$ID != 'Placeholder', ]

        all_func_res_all$Group <- paste0(all_func_res_all$Method, '\n', all_func_res_all$Type)
        all_func_res_all$Group[all_func_res_all$Group == 'SNP\nFine-mapping\n'] <- 'SuSiE'
        all_func_res_all$Group[all_func_res_all$Group == 'MAGMA\n'] <- 'MAGMA'
        all_func_res_all$Group[all_func_res_all$Group == 'Nearest\nGene\n'] <- 'Nearest\nGene'

        groups <- c('SuSiE', 'FUSION\nExpr.', 'FUSION\nSplice', 'SMR\nExpr.', 'FUSION\nProtein', 'SMR\nProtein', 'MAGMA', 'Nearest\nGene')
        groups <- groups[groups %in% all_func_res_all$Group]

        all_func_res_all$Group <- factor(all_func_res_all$Group, levels = groups)

        # Per-locus layout only if positions are available; otherwise fall back
        # to alphabetical (a note explains this via message_no_positions_mol).
        locus_mode <- identical(input$mol_layout, "locus") && nrow(mol_locus_map()) > 0

        if (locus_mode) {
          # Join locus assignments; genes without a position go to a trailing
          # "Unplaced" band ordered alphabetically at the bottom.
          lmap <- mol_locus_map()
          all_func_res_all <- merge(
            all_func_res_all, lmap[, c("ID", "BP", "Locus", "locus_order")],
            by = "ID", all.x = TRUE)
          unplaced <- is.na(all_func_res_all$locus_order)
          max_order <- suppressWarnings(max(all_func_res_all$locus_order, na.rm = TRUE))
          if (!is.finite(max_order)) max_order <- 0
          all_func_res_all$Locus[unplaced] <- "Unplaced"
          all_func_res_all$locus_order[unplaced] <- max_order + 1
          all_func_res_all$BP[unplaced] <- Inf

          # Gene order: by locus, then position within locus, then symbol.
          gene_order <- unique(all_func_res_all[order(all_func_res_all$locus_order,
                                                      all_func_res_all$BP,
                                                      as.character(all_func_res_all$ID)), "ID"])
          all_func_res_all$ID <- factor(all_func_res_all$ID, levels = rev(gene_order))
          locus_levels <- unique(all_func_res_all$Locus[order(all_func_res_all$locus_order)])
          all_func_res_all$Locus <- factor(all_func_res_all$Locus, levels = locus_levels)
        } else {
          all_func_res_all <- all_func_res_all[order(as.character(all_func_res_all$ID)), ]
          all_func_res_all$ID <- factor(all_func_res_all$ID, levels = rev(unique(as.character(all_func_res_all$ID))))
        }

        x <- c(-max(abs(all_func_res_all$Z), na.rm = T), 0, max(abs(all_func_res_all$Z), na.rm = T))
        x <- (x - min(x)) / (max(x) - min(x))

        all_func_res_all <- data.table(all_func_res_all)

        # Extra left plot.margin to absorb 45° tick label overflow past the
        # y-axis area (kept in sync with calc_plot_dims's width budget).
        left_pad_pt <- plot_dim_mol()[["left_pad_pt"]]

        # Create the heatmap plot (geoms/scale/theme shared across layouts)
        heatmap <- ggplot(data = all_func_res_all, aes(x = Panel, y = ID)) +
          theme_bw() +
          geom_point(data = all_func_res_all[all_func_res_all$Sig == T, ], aes(x = Panel, y = ID), colour = 'black', size = 5) +
          geom_point(data = all_func_res_all[all_func_res_all$Coloc == T & all_func_res_all$Sig == T, ], aes(x = Panel, y = ID), colour = 'black', shape = 15, size = 6) +
          geom_point(data = all_func_res_all[(all_func_res_all$Method == 'SNP\nFine-mapping' | all_func_res_all$Method == 'Nearest\nGene') & !(is.na(all_func_res_all$Z)), ], aes(x = Panel, y = ID), colour = '#00FF00', size = 5) +
          geom_point(data = all_func_res_all[all_func_res_all$Method != 'SNP\nFine-mapping' & all_func_res_all$Method != 'Nearest\nGene', ], aes(colour = Z), size = 4) +
          scale_colour_gradientn(colours = c("#0066FF", "#0099FF", "#FFFFFF", "#FF6666", "#FF0000"), na.value = NA, name = "Z-score", limits = c(-max(abs(all_func_res_all$Z), na.rm = T), max(abs(all_func_res_all$Z), na.rm = T)), values = x) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
          labs(x = '', y = '') +
          theme(text = element_text(size = 14),
                plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = left_pad_pt, unit = "pt"))

        if (locus_mode) {
          # Locus bands (facet rows) x method groups (facet columns). space="free"
          # sizes each band by its gene count and each column by its panel count,
          # so the manual gtable width hack used below is not needed here.
          heatmap <- heatmap +
            facet_grid(Locus ~ Group, scales = "free", space = "free") +
            theme(strip.text.y = element_text(angle = 0),
                  panel.spacing.y = grid::unit(4, "pt"))
          print(heatmap)
        } else {
          # Absolute per-facet widths matching calc_plot_dims's budget, so
          # unaffected facets (e.g. MAGMA) keep the same width when panels are
          # added to TWAS/SMR facets.
          group_widths <- vapply(groups, function(g) {
            panel_names <- unique(all_func_res_all$Panel[all_func_res_all$Group == g])
            facet_target_width_pt(g, panel_names)
          }, numeric(1))

          heatmap <- heatmap +
            facet_wrap(~ Group, nrow = 1, scales = "free_x") +
            scale_y_discrete(limits = unique(rev(all_func_res_all$ID)))

          gt <- ggplot_gtable(ggplot_build(heatmap))
          for (i in seq_along(groups)) {
            panel_l <- gt$layout$l[grepl(paste0('^panel-', i, '-1$'), gt$layout$name)]
            gt$widths[panel_l] <- grid::unit(group_widths[i], "pt")
          }
          grid.draw(gt)
        }
      } else {
        NULL
      }
    })

    ########
    # UI outputs for plot and messages
    ########

    output$mol_assoc_plot.ui <- renderUI({
      ns <- session$ns
      filtered <- mol_assoc_summary_data_filtered()
      has_real_genes <- any(filtered$ID != 'Placeholder')
      if (plot_dim_mol()[['height']] < 10000 & nrow(filtered) > 0 & has_real_genes) {
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
  })
}
