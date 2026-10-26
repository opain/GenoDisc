# Effector-gene prioritisation tab: PoPS gene scores, CALDERA per-credible-set
# causal-gene probabilities, and credible-set coding annotation (AlphaMissense +
# consequence calls). Reads the packaged blocks mol_assoc/pops, caldera and
# credset_coding (the last is a list of $summary + $variants). Pattern mirrors
# R/mod_mol_assoc.R (renderPlot + gd_legend + downloadHandler/open_plot_device);
# null-safe per sub-tab. Each sub-tab opens with a collapsible "How to read this"
# intro stating the method's scope + main caveat, then a figure, then the table.

# ---- shared UI helpers -------------------------------------------------------

# Collapsible "How to read this" intro (gd-details pattern from mod_mol_assoc.R).
.eff_howto <- function(..., title = "How to read this") {
  tags$details(
    class = "gd-details",
    tags$summary(title),
    tags$div(class = "gd-details-body",
             tags$p(class = "gd-details-intro", ...))
  )
}

# A figure block: a titled plot with a small PNG download button beneath it.
.eff_figure <- function(ns, plot_id, height = "420px") {
  tagList(
    plotOutput(ns(plot_id), height = height),
    div(style = "margin: 4px 0 2px;",
        downloadButton(ns(paste0(plot_id, "_dl")), "Download PNG",
                       class = "btn btn-default btn-sm"))
  )
}

effectorUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Effector genes",
    br(),
    p("Locus-level effector-gene prioritisation and credible-set coding annotation ",
      "for the selected GWAS. PoPS, CALDERA and the coding annotation answer three ",
      "different questions — ", strong("which gene"), " at a locus (CALDERA), ",
      strong("how plausible"), " a gene is genome-wide (PoPS), and whether a signal ",
      "is ", strong("carried by a coding variant"), " (coding) — so read them together."),
    hr(),
    tabsetPanel(
      id = ns("effector_tabs"),

      # ---------------------------------------------------------------- PoPS
      tabPanel(
        title = "PoPS", value = "pops",
        br(),
        .eff_howto(
          strong("PoPS (Polygenic Priority Score)"), " ranks every gene genome-wide ",
          "by combining its MAGMA gene association with a large published functional-",
          "genomics feature matrix. It is a ", strong("genome-wide prioritisation"),
          ", not a locus-specific or causal call: a high-PoPS gene near a GWAS signal ",
          "is a strong candidate, but PoPS alone does not prove causality or assign a ",
          "gene to a locus. The track below places each gene at its genomic position; ",
          "peaks near your associated loci are the ones to follow up."),
        .eff_figure(ns, "pops_track"),
        uiOutput(ns("pops_track_msg")),
        br(), h4("All genes"),
        DT::dataTableOutput(ns("pops_table")),
        uiOutput(ns("pops_msg")), br()
      ),

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
          "measure of causality. The ", strong("nearest gene"), " is marked as a ",
          "transparency baseline: CALDERA agreeing with it is reassuring, disagreeing ",
          "is where it adds value. Pick a locus to see its genes ranked."),
        uiOutput(ns("caldera_locus_ui")),
        .eff_figure(ns, "caldera_bar"),
        uiOutput(ns("caldera_bar_msg")),
        br(), h4("All loci (table)"),
        DT::dataTableOutput(ns("caldera_table")),
        uiOutput(ns("caldera_msg")),
        gd_legend(list(
          "caldera"      = "CALDERA causal-gene probability for this gene in this credible set (0–1, sums to 1 per locus).",
          "rank"         = "Rank of this gene within the credible set by CALDERA probability (1 = top).",
          "nearest_gene" = "TRUE if this is the gene nearest the lead variant (baseline comparator).",
          "pops"         = "PoPS score feature; dist = bp distance to the lead variant; coding = summed non-synonymous credible-set PIP."
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
      ),

      # --------------------------------------------------- Locus report card
      tabPanel(
        title = "Locus report card", value = "report",
        br(),
        .eff_howto(
          "For one locus, this pulls the three lines of evidence together: CALDERA's ",
          "gene ranking, each gene's genome-wide PoPS score, and whether the credible ",
          "set carries a protein-altering variant. The ", strong("verdict"),
          " line summarises whether they converge on a single gene and by what ",
          "mechanism (coding vs regulatory). Use it to decide the most likely effector ",
          "gene at a locus at a glance."),
        uiOutput(ns("report_locus_ui")),
        uiOutput(ns("report_verdict")),
        .eff_figure(ns, "report_caldera"),
        uiOutput(ns("report_msg")),
        br(), h4("Credible-set variants at this locus"),
        DT::dataTableOutput(ns("report_variants")),
        br()
      )
    )
  )
}

effectorServer <- function(id, gwas_data, selected_gwas, config_flags) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Read a block for the currently-selected GWAS; NULL on any problem.
    read_block <- function(block) {
      gd <- tryCatch(gwas_data(), error = function(e) NULL)
      g  <- tryCatch(selected_gwas(), error = function(e) NULL)
      if (is.null(gd) || is.null(g) || !nzchar(g)) return(NULL)
      tryCatch(gd_read(gd, g, block), error = function(e) NULL)
    }

    empty_dt <- function() DT::datatable(data.frame(), rownames = FALSE,
                                         options = list(dom = "t"))
    msg <- function(text) div(class = "gd-legend", em(text))
    show_dt <- function(df, sci_cols = NULL) {
      df <- as.data.frame(df); df[is.na(df)] <- "NA"
      DT::datatable(df, rownames = FALSE, filter = "top",
                    options = list(pageLength = 15, scrollX = TRUE,
                                   columnDefs = list(list(className = "dt-center",
                                                          targets = seq_len(ncol(df)) - 1))))
    }

    # Reactive block readers (cached per flush by Shiny).
    pops_data   <- reactive(read_block("mol_assoc/pops"))
    caldera_data <- reactive(read_block("caldera"))
    coding_data  <- reactive(read_block("credset_coding"))

    # Register a PNG download for a renderPlot-style builder. `build()` returns a
    # ggplot (or NULL); width/height in inches.
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

    # ======================================================================
    # PoPS
    # ======================================================================
    build_pops_track <- function() {
      d <- pops_data()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      d <- as.data.frame(d)
      d$PoPS_Score <- suppressWarnings(as.numeric(d$PoPS_Score))
      d <- d[is.finite(d$PoPS_Score), , drop = FALSE]
      if (nrow(d) == 0) return(NULL)
      pos <- tryCatch(
        resolve_feature_positions(d$Gene, gwas_data(), selected_gwas(), config_flags()),
        error = function(e) NULL)
      if (is.null(pos) || nrow(pos) == 0) {
        # Fallback: ranked lollipop of the top genes (no positions available).
        top <- d[order(-d$PoPS_Score), ][seq_len(min(25, nrow(d))), ]
        top$Gene <- factor(top$Gene, levels = rev(top$Gene))
        return(
          ggplot2::ggplot(top, ggplot2::aes(x = Gene, y = PoPS_Score)) +
            ggplot2::geom_segment(ggplot2::aes(xend = Gene, y = 0, yend = PoPS_Score),
                                  colour = "#9aa7c7") +
            ggplot2::geom_point(colour = "#2b5bd7", size = 2.6) +
            ggplot2::coord_flip() +
            ggplot2::labs(x = NULL, y = "PoPS score",
                          title = "Top PoPS-prioritised genes (positions unavailable)") +
            ggplot2::theme_bw())
      }
      m <- merge(d, pos, by.x = "Gene", by.y = "ID")
      m$CHR <- suppressWarnings(as.integer(m$CHR))
      m <- m[is.finite(m$CHR) & is.finite(m$BP), , drop = FALSE]
      if (nrow(m) == 0) return(NULL)
      # Cumulative x across chromosomes.
      m <- m[order(m$CHR, m$BP), ]
      chr_max <- tapply(m$BP, m$CHR, max)
      offset <- cumsum(as.numeric(c(0, head(chr_max, -1))))
      names(offset) <- names(chr_max)
      m$xpos <- m$BP + offset[as.character(m$CHR)]
      axis_df <- data.frame(CHR = names(chr_max),
                            centre = offset + chr_max / 2)
      m$chr_parity <- factor(m$CHR %% 2)
      lab <- m[order(-m$PoPS_Score), ][seq_len(min(15, nrow(m))), ]
      ggplot2::ggplot(m, ggplot2::aes(x = xpos, y = PoPS_Score, colour = chr_parity)) +
        ggplot2::geom_point(size = 1.5) +
        ggplot2::scale_colour_manual(values = c("0" = "#2b5bd7", "1" = "#7aa0f0"),
                                     guide = "none") +
        ggplot2::scale_x_continuous(breaks = axis_df$centre, labels = axis_df$CHR,
                                    expand = c(0.01, 0)) +
        ggrepel::geom_text_repel(data = lab,
                                 ggplot2::aes(label = Gene), size = 3,
                                 max.overlaps = 20, colour = "black",
                                 segment.colour = "grey60") +
        ggplot2::labs(x = "Chromosome", y = "PoPS score",
                      title = "Genome-wide PoPS scores (top genes labelled)") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank())
    }
    output$pops_track <- renderPlot({ p <- build_pops_track(); if (!is.null(p)) p })
    output$pops_track_msg <- renderUI({
      if (is.null(pops_data()) || nrow(pops_data()) == 0)
        return(msg("No PoPS results for this GWAS."))
      NULL
    })
    reg_download("pops_track", build_pops_track, w = 12, h = 5)

    output$pops_table <- DT::renderDataTable({
      d <- pops_data()
      if (is.null(d) || nrow(d) == 0) return(empty_dt())
      d <- as.data.frame(d)
      if ("PoPS_Score" %in% names(d))
        d <- d[order(-suppressWarnings(as.numeric(d$PoPS_Score))), ]
      show_dt(d)
    })
    output$pops_msg <- renderUI({
      if (is.null(pops_data()) || nrow(pops_data()) == 0)
        return(msg("No PoPS results for this GWAS."))
      NULL
    })

    # ======================================================================
    # CALDERA
    # ======================================================================
    output$caldera_locus_ui <- renderUI({
      d <- caldera_data()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      loci <- sort(unique(as.character(d$locus)))
      selectInput(ns("caldera_locus"), "Locus (lead variant):",
                  choices = loci, selected = loci[1], width = "320px")
    })

    build_caldera_bar <- function() {
      d <- caldera_data()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      d <- as.data.frame(d)
      loc <- input$caldera_locus %||% sort(unique(as.character(d$locus)))[1]
      sub <- d[as.character(d$locus) == loc, , drop = FALSE]
      if (nrow(sub) == 0) return(NULL)
      sub$caldera <- suppressWarnings(as.numeric(sub$caldera))
      sub$nearest_gene <- as.logical(sub$nearest_gene)
      sub$gene <- factor(sub$gene, levels = sub$gene[order(sub$caldera)])
      sub$fill_grp <- ifelse(sub$rank == 1, "Top-ranked", "Other")
      multi_cs <- length(unique(sub$cs)) > 1
      g <- ggplot2::ggplot(sub, ggplot2::aes(x = gene, y = caldera, fill = fill_grp)) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::geom_point(data = sub[sub$nearest_gene %in% TRUE, , drop = FALSE],
                            ggplot2::aes(x = gene, y = caldera),
                            shape = 23, size = 3, fill = "white", colour = "black",
                            inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("Top-ranked" = "#1a7f5a", "Other" = "#9aa7c7"),
                                   name = NULL) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "CALDERA causal-gene probability",
                      title = paste0("CALDERA at ", loc),
                      subtitle = "◇ = nearest gene (baseline)") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")
      if (multi_cs) g <- g + ggplot2::facet_wrap(~ cs, scales = "free_y",
                                                 labeller = ggplot2::label_both)
      g
    }
    output$caldera_bar <- renderPlot({ p <- build_caldera_bar(); if (!is.null(p)) p })
    output$caldera_bar_msg <- renderUI({
      if (is.null(caldera_data()) || nrow(caldera_data()) == 0)
        return(msg("No CALDERA results for this GWAS."))
      NULL
    })
    reg_download("caldera_bar", build_caldera_bar, w = 9, h = 6)

    output$caldera_table <- DT::renderDataTable({
      d <- caldera_data()
      if (is.null(d) || nrow(d) == 0) return(empty_dt())
      show_dt(d)
    })
    output$caldera_msg <- renderUI({
      if (is.null(caldera_data()) || nrow(caldera_data()) == 0)
        return(msg("No CALDERA results for this GWAS."))
      NULL
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
      s$label <- paste0(s$locus, " (cs", s$cs, ")")
      s$other <- pmax(s$pip_total - s$pip_protein_altering, 0)
      s <- s[order(s$frac_protein_altering, s$pip_protein_altering), ]
      s$label <- factor(s$label, levels = s$label)
      long <- rbind(
        data.frame(label = s$label, part = "Protein-altering", pip = s$pip_protein_altering),
        data.frame(label = s$label, part = "Other",            pip = s$other))
      long$part <- factor(long$part, levels = c("Other", "Protein-altering"))
      flagged <- s$label[as.logical(s$flag_protein_altering) %in% TRUE]
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
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")
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

    # ======================================================================
    # Locus report card
    # ======================================================================
    output$report_locus_ui <- renderUI({
      d <- caldera_data()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      loci <- sort(unique(as.character(d$locus)))
      selectInput(ns("report_locus"), "Locus (lead variant):",
                  choices = loci, selected = loci[1], width = "320px")
    })

    # Assemble the per-locus evidence for the chosen locus. Returns a list or NULL.
    report_bits <- reactive({
      d <- caldera_data()
      if (is.null(d) || nrow(d) == 0) return(NULL)
      d <- as.data.frame(d)
      loc <- input$report_locus %||% sort(unique(as.character(d$locus)))[1]
      cal <- d[as.character(d$locus) == loc, , drop = FALSE]
      if (nrow(cal) == 0) return(NULL)
      cal$caldera <- suppressWarnings(as.numeric(cal$caldera))
      top <- cal[which.max(cal$caldera), , drop = FALSE]
      near <- cal[as.logical(cal$nearest_gene) %in% TRUE, , drop = FALSE]
      near_gene <- if (nrow(near) > 0) near$gene[1] else NA_character_
      # Coding flag for this locus (from credset_coding$summary).
      cd <- coding_data()
      coding_gene <- NA_character_; coding_frac <- NA_real_
      if (!is.null(cd) && !is.null(cd$summary) && nrow(cd$summary) > 0) {
        cs <- as.data.frame(cd$summary)
        cs <- cs[as.character(cs$locus) == loc &
                   (as.logical(cs$flag_protein_altering) %in% TRUE), , drop = FALSE]
        if (nrow(cs) > 0) {
          cs <- cs[order(-suppressWarnings(as.numeric(cs$frac_protein_altering))), ]
          coding_gene <- cs$top_coding_gene[1]
          coding_frac <- suppressWarnings(as.numeric(cs$frac_protein_altering[1]))
        }
      }
      list(loc = loc, cal = cal, top_gene = as.character(top$gene[1]),
           top_prob = top$caldera[1], near_gene = as.character(near_gene),
           coding_gene = coding_gene, coding_frac = coding_frac)
    })

    output$report_verdict <- renderUI({
      b <- report_bits()
      if (is.null(b)) return(msg("No CALDERA results for this GWAS."))
      if (!is.na(b$coding_gene) && nzchar(b$coding_gene)) {
        verdict <- sprintf(
          "Coding mechanism: the credible set carries %.0f%% of its PIP on protein-altering variants in %s. %s is the strongest effector-gene candidate.",
          100 * (b$coding_frac %||% NA), b$coding_gene,
          if (identical(b$coding_gene, b$top_gene)) paste0(b$coding_gene, " (also CALDERA-top)") else b$coding_gene)
      } else if (!is.na(b$near_gene) && identical(b$near_gene, b$top_gene)) {
        verdict <- sprintf(
          "CALDERA top-ranks %s (p = %.2f), which is also the nearest gene — consistent prioritisation, with no coding signal (regulatory mechanism likely).",
          b$top_gene, b$top_prob)
      } else {
        verdict <- sprintf(
          "CALDERA top-ranks %s (p = %.2f) over the nearest gene%s — prioritised on functional grounds rather than proximity; no coding signal.",
          b$top_gene, b$top_prob,
          if (!is.na(b$near_gene)) paste0(" (", b$near_gene, ")") else "")
      }
      div(class = "gd-legend",
          tags$strong(paste0("Verdict for ", b$loc, ": ")), verdict)
    })

    build_report_caldera <- function() {
      b <- report_bits(); if (is.null(b)) return(NULL)
      sub <- b$cal
      sub$nearest_gene <- as.logical(sub$nearest_gene)
      sub$gene <- factor(sub$gene, levels = sub$gene[order(sub$caldera)])
      sub$fill_grp <- ifelse(sub$rank == 1, "Top-ranked", "Other")
      g <- ggplot2::ggplot(sub, ggplot2::aes(x = gene, y = caldera, fill = fill_grp)) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::geom_point(data = sub[sub$nearest_gene %in% TRUE, , drop = FALSE],
                            ggplot2::aes(x = gene, y = caldera),
                            shape = 23, size = 3, fill = "white", colour = "black",
                            inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("Top-ranked" = "#1a7f5a", "Other" = "#9aa7c7"),
                                   name = NULL) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "CALDERA causal-gene probability",
                      title = paste0("Gene ranking at ", b$loc),
                      subtitle = "◇ = nearest gene (baseline)") +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
      if (length(unique(sub$cs)) > 1)
        g <- g + ggplot2::facet_wrap(~ cs, scales = "free_y", labeller = ggplot2::label_both)
      g
    }
    output$report_caldera <- renderPlot({ p <- build_report_caldera(); if (!is.null(p)) p })
    reg_download("report_caldera", build_report_caldera, w = 9, h = 6)
    output$report_msg <- renderUI({
      if (is.null(caldera_data()) || nrow(caldera_data()) == 0)
        return(msg("No CALDERA results for this GWAS."))
      NULL
    })

    output$report_variants <- DT::renderDataTable({
      b <- report_bits(); cd <- coding_data()
      if (is.null(b) || is.null(cd) || is.null(cd$variants) || nrow(cd$variants) == 0)
        return(empty_dt())
      v <- as.data.frame(cd$variants)
      v <- v[as.character(v$locus) == b$loc, , drop = FALSE]
      if (nrow(v) == 0) return(empty_dt())
      show_dt(v)
    })
  })
}
