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
          "This flags credible sets whose signal is carried by ", strong("protein-altering"),
          " variants (missense / stop / frameshift / splice). The bar below shows, per ",
          "credible set, the share of fine-mapping PIP that sits on protein-altering ",
          "variants. A high share points to a coding mechanism and names the gene directly. ",
          "Two caveats: the flag measures ", strong("PIP share"), ", not pathogenicity; and ",
          "AlphaMissense labelling a common missense variant ", strong("\"benign\""),
          " is expected and does ", strong("not"), " mean it is non-functional."),
        .eff_figure(ns, "coding_bar"),
        uiOutput(ns("coding_bar_msg")),
        br(), h4("Per credible set"),
        DT::dataTableOutput(ns("coding_summary_table")),
        br(), h4("Per variant"),
        DT::dataTableOutput(ns("coding_variants_table")),
        uiOutput(ns("coding_msg")),
        gd_legend(list(
          "pip_protein_altering" = "Summed PIP on protein-altering (missense/stop/frameshift/splice/…) variants in the credible set.",
          "frac_protein_altering" = "Fraction of the credible set's total PIP carried by protein-altering variants.",
          "am_class"             = "AlphaMissense class for a missense variant: likely_benign / ambiguous / likely_pathogenic.",
          "am_pathogenicity"     = "AlphaMissense pathogenicity score (0–1).",
          "consequence"          = "Predicted coding consequence (bcftools csq / AlphaMissense)."
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
        near <- g$gene[as.logical(g$nearest_gene) %in% TRUE]
        near <- if (length(near)) near[1] else NA_character_
        row <- data.frame(
          locus = g$locus[1], cs = g$cs[1],
          top_gene = as.character(top$gene[1]),
          top_prob = round(top$caldera[1], 3),
          nearest_gene = as.character(near),
          nearest_agree = !is.na(near) && identical(as.character(top$gene[1]), as.character(near)),
          n_genes = nrow(g), stringsAsFactors = FALSE)
        for (cc in cs_cols) row[[cc]] <- g[[cc]][1]
        row
      })
      ov <- do.call(rbind, parts)
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
      names(g)[names(g) == "multi"] <- "Unnormalised probability"
      show_dt(g, header_tips = list(
        "Unnormalised probability" = "Raw per-gene CALDERA score before normalising to sum to 1 within the credible set.",
        "pops" = "The gene's PoPS score (genome-wide prioritisation feature)."))
    })

    # ======================================================================
    # Credible-set coding
    # ======================================================================
    build_coding_bar <- function() {
      d <- coding_data()
      if (is.null(d) || is.null(d$summary) || nrow(d$summary) == 0) return(NULL)
      s <- as.data.frame(d$summary)
      for (cc in c("pip_total", "pip_protein_altering", "frac_protein_altering"))
        if (cc %in% names(s)) s[[cc]] <- suppressWarnings(as.numeric(s[[cc]]))
      s$label <- paste0(s$locus, " (cs", s$cs, if ("L" %in% names(s)) paste0(", L", s$L) else "", ")")
      s$other <- pmax(s$pip_total - s$pip_protein_altering, 0)
      s <- s[order(s$frac_protein_altering, s$pip_protein_altering), ]
      s$label <- factor(s$label, levels = s$label)
      long <- rbind(
        data.frame(label = s$label, part = "Protein-altering", pip = s$pip_protein_altering),
        data.frame(label = s$label, part = "Other",            pip = s$other))
      long$part <- factor(long$part, levels = c("Other", "Protein-altering"))
      gene_lab <- s[as.logical(s$flag_protein_altering) %in% TRUE,
                    c("label", "pip_total", "top_coding_gene"), drop = FALSE]
      ggplot2::ggplot(long, ggplot2::aes(x = label, y = pip, fill = part)) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::scale_fill_manual(values = c("Protein-altering" = "#c0392b",
                                              "Other" = "#cfd6e6"), name = NULL) +
        { if (nrow(gene_lab) > 0)
            ggplot2::geom_text(data = gene_lab,
              ggplot2::aes(x = label, y = pip_total, label = top_coding_gene),
              hjust = -0.1, size = 3, inherit.aes = FALSE) } +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Credible-set PIP",
                      title = "Credible-set PIP carried by protein-altering variants",
                      subtitle = "Red = protein-altering share; labelled sets are flagged") +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
    }
    output$coding_bar <- renderPlot({ p <- build_coding_bar(); if (!is.null(p)) p })
    output$coding_bar_msg <- renderUI({
      d <- coding_data()
      if (is.null(d) || is.null(d$summary) || nrow(d$summary) == 0)
        return(msg("No credible-set coding annotation for this GWAS."))
      NULL
    })
    reg_download("coding_bar", build_coding_bar, w = 9, h = 7)

    output$coding_summary_table <- DT::renderDataTable({
      d <- coding_data()
      if (is.null(d) || is.null(d$summary) || nrow(d$summary) == 0) return(empty_dt())
      show_dt(d$summary)
    })
    output$coding_variants_table <- DT::renderDataTable({
      d <- coding_data()
      if (is.null(d) || is.null(d$variants) || nrow(d$variants) == 0) return(empty_dt())
      show_dt(d$variants)
    })
    output$coding_msg <- renderUI({
      d <- coding_data()
      if (is.null(d) || (is.null(d$summary) && is.null(d$variants)))
        return(msg("No credible-set coding annotation for this GWAS."))
      NULL
    })

  })
}
