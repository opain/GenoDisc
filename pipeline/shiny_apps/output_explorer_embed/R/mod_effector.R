# Effector-gene prioritisation (sub-tabs of Molecular Associations): CALDERA
# per-credible-set causal-gene probabilities and credible-set coding annotation
# (AlphaMissense + consequence calls). Reads the packaged blocks `caldera` and
# `credset_coding` (a list of $summary + $variants). PoPS no longer has its own
# tab — each gene's PoPS score is shown in the CALDERA gene table (the `pops`
# column). CALDERA is computed for both L=1 and L=10 SuSiE credible sets (the
# `L` column); the user picks which to view.

# ---- shared UI helpers -------------------------------------------------------

# Collapsible "How to read this" intro (gd-details pattern from mod_mol_assoc.R).
.eff_howto <- function(..., title = "How to read this") {
  tags$details(
    class = "gd-details",
    style = "margin-bottom: 16px;",   # breathing room before the content below
    tags$summary(title),
    tags$div(class = "gd-details-body",
             tags$p(class = "gd-details-intro", ...))
  )
}

# A figure block: a titled plot with a small PNG download button beneath it.
.eff_figure <- function(ns, plot_id, height = "420px", max_width = 900) {
  tags$div(
    style = sprintf("max-width: %dpx;", max_width),
    plotOutput(ns(plot_id), height = height),
    div(style = "margin: 4px 0 2px;",
        downloadButton(ns(paste0(plot_id, "_dl")), "Download PNG",
                       class = "btn btn-default btn-sm"))
  )
}

# Width-capped wrapper for a table (mirrors the Manhattan / MAGMA result tables).
.eff_table <- function(ns, table_id, max_width = 900) {
  tags$div(style = sprintf("max-width: %dpx;", max_width),
           DT::dataTableOutput(ns(table_id)))
}

# L=1 / L=10 radio for CALDERA views. SuSiE L=1 assumes one causal signal per
# locus (covers every fine-mapped locus); L=10 allows up to ten (only loci with
# a multi-signal decomposition appear).
.eff_L_radio <- function(ns, id) {
  radioButtons(ns(id), "SuSiE credible sets (L):",
               choices = c("L = 1 (one signal/locus)" = "1",
                           "L = 10 (up to ten signals)" = "10"),
               selected = "1", inline = TRUE)
}

# DT headerCallback that attaches hover tooltips (title attr) to named columns.
.eff_header_tips <- function(tips) {
  obj <- paste0("{", paste(sprintf("'%s':'%s'", names(tips),
                                   gsub("'", "\\\\'", unlist(tips))), collapse = ","), "}")
  DT::JS(
    "function(thead, data, start, end, display){",
    paste0("  var tips = ", obj, ";"),
    "  $(thead).find('th').each(function(){",
    "    var t = tips[$(this).text()]; if (t) { $(this).attr('title', t); }",
    "  });",
    "}")
}

# Map a bcftools-csq / AlphaMissense consequence to a broad class (for the
# stacked PIP-composition plot) and to a tidy per-variant display label.
.CSQ_CLASS <- c(
  missense = "Protein-altering", stop_gained = "Protein-altering",
  stop_lost = "Protein-altering", start_lost = "Protein-altering",
  frameshift = "Protein-altering", inframe_deletion = "Protein-altering",
  inframe_insertion = "Protein-altering", splice_acceptor = "Splice region",
  splice_donor = "Splice region", splice_region = "Splice region",
  `3_prime_utr` = "UTR", `5_prime_utr` = "UTR", intron = "Intronic",
  non_coding = "Non-coding transcript", synonymous = "Synonymous")
.CSQ_LEVELS <- c("Protein-altering", "Splice region", "UTR", "Synonymous",
                 "Intronic", "Non-coding transcript", "Intergenic / unannotated")
.CSQ_COLS <- c("Protein-altering" = "#c0392b", "Splice region" = "#e67e22",
               "UTR" = "#f1c40f", "Synonymous" = "#95a5a6", "Intronic" = "#5dade2",
               "Non-coding transcript" = "#2980b9", "Intergenic / unannotated" = "#d5dbe3")
.csq_class <- function(x) {
  cl <- unname(.CSQ_CLASS[as.character(x)])
  cl[is.na(cl)] <- "Intergenic / unannotated"   # ".", NA, not_annotated_non_snv
  factor(cl, levels = .CSQ_LEVELS)
}
.csq_tidy <- function(x) {
  x <- as.character(x); x[is.na(x) | x == "."] <- "intergenic"
  x[x == "non_coding"] <- "non-coding transcript"
  x[x == "not_annotated_non_snv"] <- "unannotated (non-SNV)"
  gsub("_", " ", x)
}
# Our AlphaMissense table stores short class labels; show the canonical ones.
.am_display <- function(x) {
  m <- c(benign = "likely_benign", pathogenic = "likely_pathogenic",
         ambiguous = "ambiguous")[as.character(x)]
  ifelse(is.na(m), as.character(x), m)
}

effectorUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Effector genes",
    br(),
    p("Locus-level effector-gene prioritisation and credible-set coding annotation ",
      "for the selected GWAS. ", strong("CALDERA"), " ranks the candidate causal gene ",
      "at each fine-mapped locus (its PoPS score is shown per gene); ",
      strong("Credible-set coding"), " flags loci whose signal is carried by a ",
      "protein-altering variant. Read them together with the fine-mapping evidence."),
    hr(),
    tabsetPanel(
      id = ns("effector_tabs"),

      # ------------------------------------------------------------- CALDERA
      tabPanel(
        title = "CALDERA", value = "caldera",
        br(),
        .eff_howto(
          strong("CALDERA"), " gives each gene in a credible set a probability of ",
          "being the causal gene, using distance to the lead variant, the gene's PoPS ",
          "score, non-synonymous credible-set PIP and the local gene count. ",
          "Probabilities are ", strong("relative within a locus"), " (they sum to 1 per ",
          "credible set), so they rank genes at that locus — they are not an absolute ",
          "measure of causality. The ", strong("overview"), " table lists every credible ",
          "set; click a row to drill into its per-gene ranking. The ", strong("nearest gene"),
          " is a transparency baseline. ", tags$em("Note: "), "CALDERA's nearest gene is ",
          "from a GENCODE protein-coding set, so it can differ from the SNP-Associations ",
          "nearest gene (a broader Ensembl set; may show an older symbol for the same gene)."),
        .eff_L_radio(ns, "caldera_L"),
        h4("All loci — overview"),
        tags$p(tags$em("Click a row to show that credible set's gene ranking below.")),
        DT::dataTableOutput(ns("caldera_overview")),   # full screen width
        uiOutput(ns("caldera_msg")),
        gd_legend(list(
          "locus / cs"           = "Lead variant (locus) and credible-set index.",
          "top_gene / top_prob"  = "The highest-probability gene in the credible set and its CALDERA probability.",
          "nearest_gene"         = "Gene nearest the lead variant (transparency baseline).",
          "nearest_agree"        = "TRUE if the top CALDERA gene is also the nearest gene.",
          "n_genes"              = "Number of genes CALDERA scored in this credible set.",
          "flag_protein_altering" = "TRUE if the credible set carries a meaningful share of PIP on protein-altering variants.",
          "pip_protein_altering / frac_protein_altering" = "Summed / fractional credible-set PIP on protein-altering variants; top_coding_* name the leading coding variant's gene, consequence and AlphaMissense class."
        ), heading = "Overview columns"),
        hr(),
        h4("Per-locus gene ranking"),
        uiOutput(ns("caldera_sel_label")),
        tags$details(class = "gd-details", style = "margin-top: 14px; margin-bottom: 12px;",
          tags$summary("Plot options"),
          tags$div(class = "gd-details-body",
            fluidRow(
              column(4,
                sliderInput(ns("cal_font"), "Font size (pt):", min = 8, max = 24,
                            value = 14, step = 1),
                selectInput(ns("cal_theme"), "Theme:",
                            choices = c("Black & white" = "bw", "Minimal" = "minimal",
                                        "Classic" = "classic", "Light" = "light"),
                            selected = "bw")),
              column(4,
                selectInput(ns("cal_dl_format"), "Download format:",
                            choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg")),
                numericInput(ns("cal_dl_dpi"), "Resolution (DPI, PNG only):",
                             value = 300, min = 72, max = 600, step = 1)),
              column(4,
                numericInput(ns("cal_dl_w"), "Width (inches):", value = 9,
                             min = 2, max = 30, step = 0.5),
                numericInput(ns("cal_dl_h"), "Height (inches, auto):", value = 6,
                             min = 1, max = 60, step = 0.5))),
            downloadButton(ns("caldera_bar_dl"), "Download plot",
                           class = "btn btn-default btn-sm")
          )),
        tags$div(style = "max-width: 1200px; margin-bottom: 20px;",
                 plotOutput(ns("caldera_bar"), height = "auto")),
        DT::dataTableOutput(ns("caldera_gene_table")),
        gd_legend(list(
          "caldera"      = "CALDERA causal-gene probability for this gene in this credible set (0–1, sums to 1 per locus).",
          "Unnormalised probability" = "The model's raw per-gene score before it is normalised to sum to 1 within the credible set; useful for comparing absolute confidence across loci.",
          "rank"         = "Rank of this gene within the credible set by CALDERA probability (1 = top).",
          "nearest_gene" = "TRUE if this is the gene nearest the lead variant (baseline comparator).",
          "pops"         = "The gene's PoPS score; dist = bp distance to the lead variant; coding = summed non-synonymous credible-set PIP."
        )), br()
      ),

      # ------------------------------------------------- Credible-set coding
      tabPanel(
        title = "Credible-set coding", value = "credset_coding",
        br(),
        .eff_howto(
          "Each credible-set variant is annotated with its predicted consequence; the plot ",
          "shows, per credible set, how its fine-mapping PIP is split across consequence classes. ",
          tags$br(), tags$b("What PIP means here: "), "the per-credible-set posterior (SuSiE alpha), ",
          "which sums to the coverage level — not the marginal PIP. ",
          tags$br(), tags$b("What a flag means: "), "coding variants are only ~1–2% of all variants, ",
          "so even a modest protein-altering share is notable; but low PIP on a missense variant can ",
          "simply mean it is in LD with a regulatory causal variant. ",
          tags$br(), tags$b("Link to CALDERA: "), "non-synonymous PIP is one of CALDERA's four features, ",
          "so a coding flag directly raises that gene's CALDERA score — the two tabs are ",
          strong("not"), " independent evidence. ",
          tags$br(), "AlphaMissense calling a common missense variant ", strong("likely_benign"),
          " is expected and does ", strong("not"), " mean it is non-functional."),
        .eff_L_radio(ns, "coding_L"),
        h4("PIP composition by consequence class"),
        tags$details(class = "gd-details", style = "margin-top: 10px; margin-bottom: 12px;",
          tags$summary("Plot options"),
          tags$div(class = "gd-details-body",
            fluidRow(
              column(4,
                sliderInput(ns("cod_font"), "Font size (pt):", min = 8, max = 24, value = 14, step = 1),
                selectInput(ns("cod_theme"), "Theme:",
                            choices = c("Black & white" = "bw", "Minimal" = "minimal",
                                        "Classic" = "classic", "Light" = "light"), selected = "bw")),
              column(4,
                selectInput(ns("cod_dl_format"), "Download format:",
                            choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg")),
                numericInput(ns("cod_dl_dpi"), "Resolution (DPI, PNG only):", value = 300, min = 72, max = 600, step = 1)),
              column(4,
                numericInput(ns("cod_dl_w"), "Width (inches):", value = 9, min = 2, max = 30, step = 0.5),
                numericInput(ns("cod_dl_h"), "Height (inches, auto):", value = 7, min = 1, max = 80, step = 0.5))),
            downloadButton(ns("coding_bar_dl"), "Download plot", class = "btn btn-default btn-sm")
          )),
        tags$div(style = "max-width: 1200px; margin-bottom: 20px;",
                 plotOutput(ns("coding_bar"), height = "auto")),
        uiOutput(ns("coding_bar_msg")),
        h4("Per credible set"),
        DT::dataTableOutput(ns("coding_summary_table")),
        gd_legend(list(
          "coding_flag" = "Graded call on the protein-altering PIP share: ≥0.5 'likely coding mechanism', 0.1–0.5 'possible coding contribution', else 'none'.",
          "pip_protein_altering / frac_protein_altering" = "Summed / fractional credible-set PIP on protein-altering variants.",
          "caldera_gene / caldera_agree" = "CALDERA's top gene for the same credible set (at the selected L) and whether it matches the top coding gene.",
          "top_coding_am_class / top_coding_am_score" = "AlphaMissense class + pathogenicity score of the top coding variant."
        ), heading = "Per-credible-set columns"),
        br(), h4("Per variant"),
        DT::dataTableOutput(ns("coding_variants_table")),
        div(style = "margin: 6px 0;",
            downloadButton(ns("coding_variants_dl"),
                           "Download variants (CSV, incl. A1/A2)",
                           class = "btn btn-default btn-sm")),
        uiOutput(ns("coding_msg")),
        gd_legend(list(
          "PIP (per-cs) / marginal PIP" = "Per-credible-set posterior (SuSiE alpha, used by the summaries above) and the marginal PIP combining all effects.",
          "consequence"   = "Predicted consequence (bcftools csq / AlphaMissense); '.' is shown as intergenic, 'non_coding' as a non-coding transcript.",
          "REF / ALT"     = "Genome reference / alternate alleles. The GWAS A1/A2 alleles are included in the download only.",
          "am_class / am_pathogenicity" = "AlphaMissense class (likely_benign / ambiguous / likely_pathogenic) + score for missense variants."
        )), br()
      )
    )
  )
}

effectorServer <- function(id, gwas_data, selected_gwas, config_flags) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    read_block <- function(block) {
      gd <- tryCatch(gwas_data(), error = function(e) NULL)
      g  <- tryCatch(selected_gwas(), error = function(e) NULL)
      if (is.null(gd) || is.null(g) || !nzchar(g)) return(NULL)
      tryCatch(gd_read(gd, g, block), error = function(e) NULL)
    }

    empty_dt <- function() DT::datatable(data.frame(), rownames = FALSE,
                                         selection = "none", options = list(dom = "t"))
    msg <- function(text) div(class = "gd-legend", em(text))
    show_dt <- function(df, header_tips = NULL, selection = "none") {
      df <- as.data.frame(df)
      # Round probability/score columns to 3 dp (CALDERA tables otherwise show
      # many decimal places). Integer-valued columns (BP, dist) are unaffected.
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) round(x, 3))
      df[is.na(df)] <- "NA"
      opts <- list(pageLength = 15, scrollX = TRUE,
                   columnDefs = list(list(className = "dt-center",
                                          targets = seq_len(ncol(df)) - 1)))
      if (!is.null(header_tips)) opts$headerCallback <- .eff_header_tips(header_tips)
      DT::datatable(df, rownames = FALSE, filter = "top", selection = selection, options = opts)
    }

    caldera_data <- reactive(read_block("caldera"))
    coding_data  <- reactive(read_block("credset_coding"))

    # CALDERA rows for the chosen L (falls back to all rows if there is no L
    # column, e.g. a bundle built before the L split).
    caldera_at <- function(Lval) {
      d <- caldera_data()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      d <- as.data.frame(d)
      if ("L" %in% names(d)) {
        d <- d[suppressWarnings(as.integer(d$L)) == as.integer(Lval), , drop = FALSE]
      }
      if (nrow(d) == 0) return(NULL)
      d$caldera <- suppressWarnings(as.numeric(d$caldera))
      d
    }
    caldera_L_sel    <- reactive(caldera_at(input$caldera_L %||% "1"))

    reg_download <- function(plot_id, build, w = 10, h = 6) {
      output[[paste0(plot_id, "_dl")]] <- downloadHandler(
        filename = function() sprintf("effector_%s.png", plot_id),
        content = function(file) {
          p <- tryCatch(build(), error = function(e) NULL)
          open_plot_device(file, fmt = "png", w = w, h = h, units = "in", dpi = 150)
          if (!is.null(p)) print(p) else grid::grid.text("No data to plot.")
          grDevices::dev.off()
        }
      )
    }

    # Draw a CALDERA bar from a gene-level subset `sub` (one credible set).
    draw_caldera_bar <- function(sub, title, font_size = 14,
                                 theme_fn = ggplot2::theme_bw) {
      if (is.null(sub) || nrow(sub) == 0) return(NULL)
      sub$nearest_gene <- as.logical(sub$nearest_gene)
      sub$gene <- factor(sub$gene, levels = unique(sub$gene[order(sub$caldera)]))
      sub$fill_grp <- ifelse(sub$rank == 1, "Top-ranked", "Other")
      ggplot2::ggplot(sub, ggplot2::aes(x = caldera, y = gene, fill = fill_grp)) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::geom_point(data = sub[sub$nearest_gene %in% TRUE, , drop = FALSE],
                            ggplot2::aes(x = caldera, y = gene),
                            shape = 23, size = 3, fill = "white", colour = "black",
                            inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("Top-ranked" = "#1a7f5a", "Other" = "#9aa7c7"),
                                   name = NULL) +
        ggplot2::labs(x = "CALDERA causal-gene probability", y = NULL,
                      title = title, subtitle = "◇ = nearest gene (baseline)") +
        theme_fn(base_size = font_size) + ggplot2::theme(legend.position = "bottom")
    }

    # ======================================================================
    # CALDERA — overview (all loci) + per-locus gene ranking
    # ======================================================================

    # One row per (locus, cs): top gene + probability, nearest-gene agreement,
    # and the credible-set-level coding columns (moved off the gene table).
    caldera_overview_data <- reactive({
      d <- caldera_L_sel(); if (is.null(d)) return(NULL)
      cs_cols <- intersect(c("pip_protein_altering", "frac_protein_altering",
                             "top_coding_gene", "top_coding_consequence",
                             "top_coding_am_class", "flag_protein_altering"), names(d))
      parts <- lapply(split(d, list(d$locus, d$cs), drop = TRUE), function(g) {
        top  <- g[which.max(g$caldera), , drop = FALSE]
        near_idx <- which(as.logical(g$nearest_gene) %in% TRUE)
        near <- if (length(near_idx)) as.character(g$gene[near_idx[1]]) else NA_character_
        # Canonicalise the displayed symbols and decide nearest-agreement on
        # Ensembl identity (falling back to canonical symbol) so a symbol-vintage
        # difference between the two does not read as a disagreement.
        top_ens  <- if ("ensgid" %in% names(g)) as.character(top$ensgid[1]) else NA_character_
        near_ens <- if (length(near_idx) && "ensgid" %in% names(g))
          as.character(g$ensgid[near_idx[1]]) else NA_character_
        h_top  <- harmonise_genes(as.character(top$gene[1]))
        h_near <- harmonise_genes(near)
        top_id  <- if (!is.na(top_ens)  && top_ens  != "") top_ens  else h_top$ensembl_id
        near_id <- if (!is.na(near_ens) && near_ens != "") near_ens else h_near$ensembl_id
        agree <- (!is.na(top_id) && !is.na(near_id) && top_id == near_id) ||
                 (!is.na(near) && !is.na(h_top$symbol) && identical(h_top$symbol, h_near$symbol))
        top_disp  <- if (!is.na(h_top$symbol)  && !grepl("^ENSG", h_top$symbol))  h_top$symbol  else as.character(top$gene[1])
        near_disp <- if (!is.na(near) && !is.na(h_near$symbol) && !grepl("^ENSG", h_near$symbol)) h_near$symbol else near
        row <- data.frame(
          locus = g$locus[1], cs = g$cs[1],
          top_gene = top_disp,
          top_prob = round(top$caldera[1], 3),
          nearest_gene = near_disp,
          nearest_agree = isTRUE(agree),
          n_genes = nrow(g), stringsAsFactors = FALSE)
        for (cc in cs_cols) row[[cc]] <- g[[cc]][1]
        row
      })
      ov <- do.call(rbind, parts)
      if ("top_coding_gene" %in% names(ov)) {
        h <- harmonise_genes(as.character(ov$top_coding_gene))
        ov$top_coding_gene <- ifelse(!is.na(h$symbol) & !grepl("^ENSG", h$symbol),
                                     h$symbol, as.character(ov$top_coding_gene))
      }
      # Sort here (not in the render) so the row index from a table click maps to
      # the same row the drill-down uses.
      ov[order(-ov$top_prob), , drop = FALSE]
    })

    # Single-row-selectable overview: clicking a row drives the drill-down.
    output$caldera_overview <- DT::renderDataTable({
      d <- caldera_overview_data()
      if (is.null(d) || nrow(d) == 0) return(empty_dt())
      show_dt(d, selection = "single", header_tips = list(
        "nearest_agree" = "TRUE if CALDERA's top gene is also the nearest gene.",
        "top_prob"      = "CALDERA probability of the top-ranked gene in this credible set.",
        "frac_protein_altering" = "Fraction of credible-set PIP on protein-altering variants."))
    })
    output$caldera_msg <- renderUI({
      if (is.null(caldera_L_sel()))
        return(msg("No CALDERA results for this GWAS at the selected L."))
      NULL
    })

    # The selected overview row (defaults to the first row = highest-probability
    # credible set when nothing is clicked). Returns one (locus, cs) row or NULL.
    caldera_sel <- reactive({
      ov <- caldera_overview_data(); if (is.null(ov) || nrow(ov) == 0) return(NULL)
      i <- input$caldera_overview_rows_selected
      if (is.null(i) || length(i) == 0 || i < 1 || i > nrow(ov)) i <- 1L
      ov[i, , drop = FALSE]
    })

    output$caldera_sel_label <- renderUI({
      r <- caldera_sel(); if (is.null(r)) return(NULL)
      div(class = "gd-legend",
          "Showing credible set: ", tags$strong(as.character(r$locus)),
          sprintf(" — cs %s (L = %s)", r$cs, input$caldera_L %||% "1"),
          tags$span(style = "color:#888;", "  (click another overview row to change)"))
    })

    # One credible set (the selected locus + cs) at a time.
    caldera_locus_sub <- reactive({
      r <- caldera_sel(); d <- caldera_L_sel()
      if (is.null(r) || is.null(d)) return(NULL)
      sub <- d[as.character(d$locus) == as.character(r$locus) & d$cs == r$cs, , drop = FALSE]
      if (nrow(sub) == 0) NULL else sub
    })
    cal_n_genes  <- reactive({ s <- caldera_locus_sub(); if (is.null(s)) 0L else nrow(s) })
    # ~26 px/gene (roughly the Molecular-Associations heatmap spacing) + margins.
    cal_bar_h_px <- reactive(max(220, 96 + 26 * cal_n_genes()))

    build_caldera_bar <- function() {
      sub <- caldera_locus_sub(); if (is.null(sub)) return(NULL)
      r <- caldera_sel()
      ttl <- if (!is.null(r)) sprintf("CALDERA at %s (cs %s)", r$locus, r$cs) else "CALDERA"
      draw_caldera_bar(sub, ttl, font_size = input$cal_font %||% 14,
                       theme_fn = mol_theme_fn(input$cal_theme %||% "bw"))
    }
    output$caldera_bar <- renderPlot({ p <- build_caldera_bar(); if (!is.null(p)) p },
                                     height = function() cal_bar_h_px())

    # Default the download height to the gene count on each new selection.
    observeEvent(caldera_locus_sub(), {
      updateNumericInput(session, "cal_dl_h",
                         value = round(max(3, 1.6 + 0.3 * cal_n_genes()), 1))
    })
    output$caldera_bar_dl <- downloadHandler(
      filename = function() {
        r <- caldera_sel()
        tag <- if (!is.null(r)) paste0(r$locus, "_cs", r$cs) else "plot"
        sprintf("caldera_%s.%s", tag, input$cal_dl_format %||% "png")
      },
      content = function(file) {
        p <- tryCatch(build_caldera_bar(), error = function(e) NULL)
        open_plot_device(file, fmt = input$cal_dl_format %||% "png",
                         w = input$cal_dl_w %||% 9, h = input$cal_dl_h %||% 6,
                         units = "in", dpi = input$cal_dl_dpi %||% 300)
        if (!is.null(p)) print(p) else grid::grid.text("No data to plot.")
        grDevices::dev.off()
      })

    # Gene-level table (credible-set-level coding columns live in the overview).
    output$caldera_gene_table <- DT::renderDataTable({
      sub <- caldera_locus_sub(); if (is.null(sub)) return(empty_dt())
      cols <- intersect(c("locus", "cs", "gene", "ensgid", "caldera", "multi",
                          "rank", "nearest_gene", "pops", "dist", "coding"), names(sub))
      g <- sub[order(-sub$caldera), cols, drop = FALSE]
      if ("gene" %in% names(g)) {
        key <- if ("ensgid" %in% names(g))
          ifelse(!is.na(g$ensgid) & g$ensgid != "", as.character(g$ensgid), as.character(g$gene))
          else as.character(g$gene)
        h <- harmonise_genes(key)
        g$gene <- ifelse(!is.na(h$symbol) & !grepl("^ENSG", h$symbol), h$symbol, as.character(g$gene))
      }
      names(g)[names(g) == "multi"] <- "Unnormalised probability"
      show_dt(g, header_tips = list(
        "Unnormalised probability" = "Raw per-gene CALDERA score before normalising to sum to 1 within the credible set.",
        "pops" = "The gene's PoPS score (genome-wide prioritisation feature)."))
    })

    # ======================================================================
    # Credible-set coding
    # ======================================================================
    # L-filtered coding summary / variants (coding_L radio).
    coding_part <- function(which) {
      d <- coding_data(); if (is.null(d) || is.null(d[[which]])) return(NULL)
      p <- as.data.frame(d[[which]]); if (nrow(p) == 0) return(NULL)
      if ("L" %in% names(p))
        p <- p[suppressWarnings(as.integer(p$L)) == as.integer(input$coding_L %||% "1"), , drop = FALSE]
      if (nrow(p) == 0) NULL else p
    }
    coding_summary  <- reactive(coding_part("summary"))
    coding_variants <- reactive(coding_part("variants"))

    # CALDERA's top gene per (locus, cs) at the selected L (for the agreement column).
    coding_caldera_top <- reactive({
      d <- caldera_data(); if (is.null(d) || nrow(d) == 0) return(NULL)
      d <- as.data.frame(d)
      if ("L" %in% names(d))
        d <- d[suppressWarnings(as.integer(d$L)) == as.integer(input$coding_L %||% "1"), , drop = FALSE]
      if (nrow(d) == 0) return(NULL)
      d$caldera <- suppressWarnings(as.numeric(d$caldera))
      do.call(rbind, lapply(split(d, list(d$locus, d$cs), drop = TRUE), function(x) {
        r <- x[which.max(x$caldera), ]
        data.frame(locus = as.character(r$locus), cs = r$cs,
                   caldera_gene = as.character(r$gene),
                   caldera_ensgid = if ("ensgid" %in% names(r)) as.character(r$ensgid) else NA_character_,
                   stringsAsFactors = FALSE)
      }))
    })

    # ---- stacked PIP-by-consequence-class bar ----
    coding_bar_data <- reactive({
      v <- coding_variants(); if (is.null(v)) return(NULL)
      v$PIP <- suppressWarnings(as.numeric(v$PIP))
      v <- v[is.finite(v$PIP), , drop = FALSE]; if (nrow(v) == 0) return(NULL)
      v$class <- .csq_class(v$consequence)
      v$label <- paste0(v$locus, " (cs", v$cs, ")")
      v
    })
    cod_n_sets   <- reactive({ v <- coding_bar_data(); if (is.null(v)) 0L else length(unique(v$label)) })
    cod_bar_h_px <- reactive(max(240, 90 + 22 * cod_n_sets()))

    build_coding_bar <- function() {
      v <- coding_bar_data(); if (is.null(v)) return(NULL)
      agg <- aggregate(PIP ~ label + class, data = v, FUN = sum)
      # Order by protein-altering share. Compute PA-PIP over the full variant set
      # (zeroing non-PA) so its names align 1:1 with `tot` (subset-indexing can
      # drop/misalign names and collapse labels into a spurious NA bar).
      pav <- v$PIP; pav[v$class != "Protein-altering"] <- 0
      tot <- tapply(v$PIP, v$label, sum)
      pa  <- tapply(pav,   v$label, sum)
      frac <- as.numeric(pa) / as.numeric(tot)
      agg$label <- factor(agg$label, levels = names(tot)[order(frac)])  # highest PA share at top
      agg <- agg[!is.na(agg$label), , drop = FALSE]
      ggplot2::ggplot(agg, ggplot2::aes(x = PIP, y = label, fill = class)) +
        ggplot2::geom_col(width = 0.76) +
        ggplot2::scale_fill_manual(values = .CSQ_COLS, drop = FALSE, name = "Consequence class") +
        ggplot2::labs(x = "Credible-set PIP", y = NULL,
                      title = "PIP composition by consequence class",
                      subtitle = "One bar per credible set; sorted by protein-altering share (top = highest)") +
        mol_theme_fn(input$cod_theme %||% "bw")(base_size = input$cod_font %||% 14) +
        ggplot2::theme(legend.position = "bottom")
    }
    output$coding_bar <- renderPlot({ p <- build_coding_bar(); if (!is.null(p)) p },
                                    height = function() cod_bar_h_px())
    output$coding_bar_msg <- renderUI({
      if (is.null(coding_summary()))
        return(msg("No credible-set coding annotation for this GWAS at the selected L."))
      NULL
    })
    observeEvent(coding_bar_data(), {
      updateNumericInput(session, "cod_dl_h", value = round(max(3, 1.5 + 0.28 * cod_n_sets()), 1))
    })
    output$coding_bar_dl <- downloadHandler(
      filename = function() sprintf("coding_composition_L%s.%s",
                 input$coding_L %||% "1", input$cod_dl_format %||% "png"),
      content = function(file) {
        p <- tryCatch(build_coding_bar(), error = function(e) NULL)
        open_plot_device(file, fmt = input$cod_dl_format %||% "png",
                         w = input$cod_dl_w %||% 9, h = input$cod_dl_h %||% 7,
                         units = "in", dpi = input$cod_dl_dpi %||% 300)
        if (!is.null(p)) print(p) else grid::grid.text("No data to plot.")
        grDevices::dev.off()
      })

    # ---- per-credible-set summary table (graded flag + CALDERA agreement + AM score) ----
    output$coding_summary_table <- DT::renderDataTable({
      s <- coding_summary(); if (is.null(s)) return(empty_dt())
      for (cc in c("pip_total","pip_protein_altering","frac_protein_altering","top_coding_pip"))
        if (cc %in% names(s)) s[[cc]] <- suppressWarnings(as.numeric(s[[cc]]))
      v <- coding_variants()
      if (!is.null(v) && "top_coding_snp" %in% names(s) && "am_pathogenicity" %in% names(v)) {
        am <- v[!duplicated(v$SNP), c("SNP","am_pathogenicity")]
        s$top_coding_am_score <- round(suppressWarnings(as.numeric(
          am$am_pathogenicity[match(s$top_coding_snp, am$SNP)])), 3)
      }
      if ("top_coding_am_class" %in% names(s)) s$top_coding_am_class <- .am_display(s$top_coding_am_class)
      ct <- coding_caldera_top()
      if (!is.null(ct)) {
        i <- match(paste(s$locus, s$cs), paste(ct$locus, ct$cs))
        caldera_gene_raw <- ct$caldera_gene[i]
        caldera_ensgid   <- ct$caldera_ensgid[i]
        # Decide agreement on canonical gene identity, not a raw symbol compare,
        # so a GENCODE-vs-csq symbol-vintage mismatch (e.g. CFAP410 vs C21orf2)
        # for the same gene still counts as agreement. Prefer the Ensembl id;
        # fall back to canonical symbol, then the old raw-string compare.
        h_cal  <- harmonise_genes(caldera_gene_raw)
        h_code <- harmonise_genes(as.character(s$top_coding_gene))
        cal_id  <- ifelse(!is.na(caldera_ensgid) & caldera_ensgid != "",
                          caldera_ensgid, h_cal$ensembl_id)
        code_id <- h_code$ensembl_id
        agree_id  <- !is.na(cal_id) & !is.na(code_id) & cal_id == code_id
        agree_sym <- !is.na(h_cal$symbol) & nzchar(as.character(h_code$symbol)) &
          h_cal$symbol == h_code$symbol
        agree_raw <- !is.na(caldera_gene_raw) & nzchar(as.character(s$top_coding_gene)) &
          (as.character(caldera_gene_raw) == as.character(s$top_coding_gene))
        agree <- ifelse(!is.na(cal_id) & !is.na(code_id), agree_id, agree_sym | agree_raw)
        agree[is.na(agree)] <- FALSE
        # Display canonical symbols on both sides for consistency with other tabs.
        s$caldera_gene <- ifelse(!is.na(h_cal$symbol) & !grepl("^ENSG", h_cal$symbol),
                                 h_cal$symbol, caldera_gene_raw)
        s$top_coding_gene <- ifelse(!is.na(h_code$symbol) & !grepl("^ENSG", h_code$symbol),
                                    h_code$symbol, as.character(s$top_coding_gene))
        s$caldera_agree <- agree
      }
      cols <- intersect(c("locus","cs","L","n_variants","pip_total","pip_protein_altering",
                          "frac_protein_altering","coding_flag","top_coding_gene",
                          "caldera_gene","caldera_agree","top_coding_consequence",
                          "top_coding_am_class","top_coding_am_score","top_coding_snp"), names(s))
      s <- s[order(-s$frac_protein_altering), cols, drop = FALSE]
      show_dt(s, header_tips = list(
        "caldera_agree" = "TRUE if CALDERA's top gene matches the top coding gene for this credible set.",
        "coding_flag"   = "Graded protein-altering call on the PIP share."))
    })

    # ---- per-variant table (REF/ALT default; per-cs alpha + marginal PIP) ----
    output$coding_variants_table <- DT::renderDataTable({
      v <- coding_variants(); if (is.null(v)) return(empty_dt())
      v$consequence <- .csq_tidy(v$consequence)
      if ("am_class" %in% names(v)) v$am_class <- .am_display(v$am_class)
      disp <- intersect(c("locus","cs","SNP","CHR","BP","REF","ALT","PIP","marginal_pip",
                          "consequence","gene","am_class","am_pathogenicity","protein_altering"), names(v))
      v2 <- v[order(v$locus, v$cs, -suppressWarnings(as.numeric(v$PIP))), disp, drop = FALSE]
      names(v2)[names(v2) == "PIP"] <- "PIP (per-cs)"
      names(v2)[names(v2) == "marginal_pip"] <- "marginal PIP"
      show_dt(v2, header_tips = list(
        "PIP (per-cs)" = "Per-credible-set posterior (SuSiE alpha) — the value the summaries use.",
        "marginal PIP" = "Marginal PIP combining all SuSiE effects (can exceed the per-cs value)."))
    })
    output$coding_variants_dl <- downloadHandler(
      filename = function() sprintf("credset_coding_variants_L%s.csv", input$coding_L %||% "1"),
      content = function(file) {
        v <- coding_variants()
        if (is.null(v)) v <- data.frame()
        utils::write.csv(v, file, row.names = FALSE)   # full columns incl. A1/A2
      })
    output$coding_msg <- renderUI({
      d <- coding_data()
      if (is.null(d) || (is.null(d$summary) && is.null(d$variants)))
        return(msg("No credible-set coding annotation for this GWAS."))
      NULL
    })

  })
}
