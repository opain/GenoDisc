
# Friendly labels + a stable, reasonably colour-blind-safe palette for the QC
# removal reasons recorded by the GenoUtils cleaner (gwas_qc$manhattan_raw$reason
# and gwas_qc$lost_signals). Keyed by the reason code; used by both the
# unfiltered-Manhattan legend and the lost-signals tables.
GD_REMOVED_REASON_LABELS <- c(
  not_in_reference  = "Not in reference panel",
  not_biallelic_snp = "Not a biallelic SNP",
  reported_maf      = "Reported MAF below threshold",
  reference_maf     = "Reference MAF below threshold",
  maf_discordance   = "MAF discordant with reference",
  bad_p             = "Out-of-bounds P-value",
  duplicate_id      = "Duplicate variant ID",
  n_outlier         = "Sample-size outlier",
  se_zero           = "Standard error = 0"
)
GD_REMOVED_REASON_COLS <- c(
  not_in_reference  = "#e41a1c",
  not_biallelic_snp = "#377eb8",
  reported_maf      = "#4daf4a",
  reference_maf     = "#984ea3",
  maf_discordance   = "#ff7f00",
  bad_p             = "#a65628",
  duplicate_id      = "#f781bf",
  n_outlier         = "#17becf",
  se_zero           = "#999999"
)
gd_reason_label <- function(x) {
  lab <- GD_REMOVED_REASON_LABELS[as.character(x)]
  ifelse(is.na(lab), as.character(x), unname(lab))
}

# Unfiltered (raw-input, pre-reference) Manhattan plot with QC-removed variants
# highlighted by reason. Mirrors build_manhattan_plot() in mod_snp_assoc.R but
# colours by removed/reason instead of clump category. `mr` is
# gwas_qc$manhattan_raw: columns CHR, BP, P, SNP, removed (logical), reason
# (char, NA for retained). Offsets are computed here (the raw data is not
# pre-decorated with cum_bp/neglogP the way manhattan_data is).
build_raw_manhattan_plot <- function(mr,
                                     lost_signals        = NULL,
                                     sig_threshold       = 5e-8,
                                     removed_p_thresh    = 1e-5,
                                     display_p_threshold = 1e-3,
                                     point_size          = 0.9,
                                     font_size           = 14,
                                     theme_fn            = ggplot2::theme_bw,
                                     title               = "",
                                     col_chr_odd         = "#8c8c8c",
                                     col_chr_even        = "#4d4d4d") {
  if (is.null(mr) || !NROW(mr)) return(NULL)
  ss <- data.table::as.data.table(data.table::copy(mr))
  if (!all(c("CHR", "BP", "P") %in% names(ss))) return(NULL)
  ss[, CHR := suppressWarnings(as.numeric(CHR))]
  ss <- ss[!is.na(CHR) & !is.na(BP) & !is.na(P) & P > 0]
  if (!nrow(ss)) return(NULL)
  data.table::setorder(ss, CHR, BP)

  # Cumulative-BP offsets from the (thinned) data — same construction as
  # scripts/plot_manhattan.R, approximate chromosome ends are fine for an overview.
  chr_lens <- ss[, .(chr_max = max(BP)), by = CHR][order(CHR)]
  chr_lens[, cum_offset := cumsum(as.numeric(chr_max)) - as.numeric(chr_max)]
  chr_lens[, label_pos  := cum_offset + as.numeric(chr_max) / 2]
  ss <- merge(ss, chr_lens[, .(CHR, cum_offset)], by = "CHR")
  ss[, cum_bp  := BP + cum_offset]
  ss[, neglogP := -log10(P)]
  ss[, chr_bg  := ifelse(CHR %% 2 == 1, col_chr_odd, col_chr_even)]

  if (!"removed" %in% names(ss)) ss[, removed := FALSE]
  ss[, removed := as.logical(removed)]
  ss[is.na(removed), removed := FALSE]

  # Only plot the signal region (P < display_p_threshold), matching the
  # SNP-associations Manhattan. Offsets above were computed from the full thinned
  # data so the x-axis still spans each chromosome; the raw data also carries a
  # background subsample for that purpose, which we drop from the plot here.
  ss <- ss[P < display_p_threshold]
  if (!nrow(ss)) return(NULL)

  rem <- if ("reason" %in% names(ss)) ss[removed == TRUE & !is.na(reason)] else ss[0]
  if (nrow(rem) > 0) {
    # Only the reasons actually present, in canonical order, so the legend
    # doesn't list the full vocabulary when just a couple apply.
    present <- intersect(names(GD_REMOVED_REASON_LABELS), unique(as.character(rem$reason)))
    rem[, reason := factor(reason, levels = present)]
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = ss[removed == FALSE],
      ggplot2::aes(x = cum_bp, y = neglogP, colour = chr_bg),
      size = point_size) +
    ggplot2::scale_colour_identity()

  if (nrow(rem) > 0) {
    lv <- levels(rem$reason)   # present reasons only (subset values/labels so the legend matches)
    p <- p +
      ggplot2::geom_point(
        data = rem,
        ggplot2::aes(x = cum_bp, y = neglogP, fill = reason),
        shape = 21, colour = "grey20", stroke = 0.3, size = point_size * 2.4) +
      ggplot2::scale_fill_manual(
        values = GD_REMOVED_REASON_COLS[lv], labels = GD_REMOVED_REASON_LABELS[lv],
        name = "Removed — reason", drop = TRUE, na.translate = FALSE)
  }

  p +
    ggplot2::geom_hline(yintercept = -log10(sig_threshold),
                        colour = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(removed_p_thresh),
                        colour = "blue", linetype = "dotted") +
    ggplot2::scale_x_continuous(breaks = chr_lens$label_pos, labels = chr_lens$CHR,
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::coord_cartesian(ylim = c(-log10(display_p_threshold), NA), clip = "off") +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](italic(P))),
                  title = if (nzchar(title)) title else NULL) +
    theme_fn(base_size = font_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position    = "bottom")
}

# Total variants removed at each QC step, parsed from the cleaner log's running
# "... N variants remain." counts (drop = previous remaining - current). Covers
# ALL removed variants (the per-variant removed.gz only keeps P < threshold).
# The pre-reference losses (not-a-SNP + not-in-panel) are reported together as
# `not_in_reference`, matching how the cleaner logs them.
gd_parse_removed_totals <- function(log_lines) {
  keys <- c("not_in_reference", "info", "reported_maf", "reference_maf",
            "maf_discordance", "bad_p", "duplicate_id", "n_outlier", "se_zero")
  specs <- c(
    not_in_reference = "After matching variants to the reference, ([0-9]+) variants remain",
    info             = "After removal of SNPs with INFO < [^,]*, ([0-9]+) variants remain",
    reported_maf     = "After removal of SNPs with reported MAF < [^,]*, ([0-9]+) variants remain",
    reference_maf    = "After removal of SNPs with reference MAF < [^,]*, ([0-9]+) variants remain",
    maf_discordance  = "After removal of SNPs with absolute (?:folded-)?MAF difference [^,]*, ([0-9]+) variants remain",
    bad_p            = "After removal of SNPs with out-of-bound P values, ([0-9]+) variants remain",
    duplicate_id     = "After removal of SNPs with duplicate IDs, ([0-9]+) variants remain",
    n_outlier        = "After removal of SNPs with N > [^,]*, ([0-9]+) variants remain",
    se_zero          = "After removal of SNPs with SE == 0, ([0-9]+) variants remain")
  getnum <- function(re) {
    hits <- grep(re, log_lines, perl = TRUE, value = TRUE)
    if (!length(hits)) return(NA_real_)
    m <- regmatches(hits[1], regexec(re, hits[1], perl = TRUE))[[1]]
    suppressWarnings(as.numeric(m[2]))
  }
  res  <- stats::setNames(rep(NA_real_, length(keys)), keys)
  prev <- getnum("GWAS contains ([0-9]+) variants")   # input total
  for (k in keys) {
    rem <- getnum(specs[[k]])
    if (!is.na(rem)) {
      res[k] <- if (!is.na(prev)) max(0, prev - rem) else NA_real_
      prev <- rem               # step ran: advance the running count
    } else {
      res[k] <- 0               # step's line absent (column not present) => nothing removed
    }
  }
  res
}

# Build the "removed by reason" summary: total removed (from the cleaner log) and
# the subset with P < 5e-8 (from lost_signals$counts_by_reason, folding the
# not-a-SNP split back into not_in_reference). Returns a display data.frame.
gd_build_removed_counts_table <- function(cleaner_log, lost_signals) {
  if (is.null(cleaner_log) || !length(cleaner_log)) return(NULL)
  tot <- gd_parse_removed_totals(cleaner_log)
  sig <- stats::setNames(rep(0L, length(tot)), names(tot))
  cb  <- lost_signals$counts_by_reason
  if (!is.null(cb) && nrow(cb) > 0) {
    for (i in seq_len(nrow(cb))) {
      k <- as.character(cb$reason[i])
      if (identical(k, "not_biallelic_snp")) k <- "not_in_reference"
      if (k %in% names(sig)) sig[k] <- sig[k] + as.integer(cb$n[i])
    }
  }
  keep <- names(tot)[(!is.na(tot) & tot > 0) | sig > 0]
  if (!length(keep)) return(NULL)
  labs <- gd_reason_label(keep)
  labs[keep == "not_in_reference"] <- "Not in reference panel / not a SNP"
  fmt <- function(x) ifelse(is.na(x), "—", formatC(as.integer(x), format = "d", big.mark = ","))
  out <- data.frame(
    Reason = labs,
    `Total removed` = fmt(tot[keep]),
    `Removed (P < 5e-8)` = fmt(sig[keep]),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  total_all <- sum(tot[keep], na.rm = TRUE)
  out <- rbind(out, data.frame(
    Reason = "Total",
    `Total removed` = fmt(total_all),
    `Removed (P < 5e-8)` = fmt(sum(sig[keep])),
    check.names = FALSE, stringsAsFactors = FALSE))
  out
}

gwasQcUI <- function(id) {
  ns <- NS(id)
  # Per-tab GWAS picker — rendered via `uiOutput` so the choice list is
  # baked in at widget-creation time. updateSelectInput sent to a
  # selectize-backed selectInput before the client binding is ready
  # silently drops the choices update. No wrapping div here so that a
  # NULL return in single-GWAS mode leaves zero residual space.
  gwas_picker <- function(output_id) uiOutput(ns(output_id))
  tabPanel(
    title = "GWAS QC",
    br(),
    p("This tab shows key quality-control statistics for the loaded GWAS. ",
      "The QC Summary table covers every selected GWAS; the other sub-tabs ",
      "are per-GWAS and expose a GWAS selector at the top so you can flick ",
      "between them without changing your bundle-wide selection."),
    hr(),
    tabsetPanel(
      id = ns("gwas_qc_tabs"),
      tabPanel(
        title = "QC Summary",
        value = "qc_summary",
        br(),
        div(style = "max-width: 700px;", tableOutput(ns("qc_table"))),
        uiOutput(ns("qc_legend"))
      ),
      tabPanel(
        title = "Allele Frequency Plot",
        value = "maf_plot",
        br(),
        gwas_picker("maf_gwas_ui"),
        uiOutput(ns("maf_plot_ui"))
      ),
      tabPanel(
        title = "QQ Plot",
        value = "qq_plot",
        br(),
        gwas_picker("qq_gwas_ui"),
        uiOutput(ns("qq_plot_ui"))
      ),
      tabPanel(
        title = "Removed Variants",
        value = "removed_variants",
        br(),
        gwas_picker("rawmh_gwas_ui"),
        uiOutput(ns("raw_manhattan_ui"))
      ),
      # SNP-h² and rG have moved to their own top-level tab
      # (SNP-h² & rG). Removing them keeps this tab focused on QC.
      tabPanel(
        title = "Sumstat QC Log",
        value = "cleaner_log",
        br(),
        gwas_picker("log_gwas_ui"),
        uiOutput(ns("cleaner_log_ui"))
      )
    )
  )
}

gwasQcServer <- function(id, gwas_data, selected_gwas, gwas_list, config_flags,
                          selected_gwas_multi = NULL,
                          comparison_mode = NULL) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # Per-tab GWAS pickers rendered via renderUI so the choice list is
    # baked in at widget-creation time. Preserves the user's per-tab
    # pick across bundle-selection changes when the pick is still valid.
    .render_per_gwas_picker <- function(output_id, input_id, label = "GWAS:") {
      output[[output_id]] <- renderUI({
        gwas_vec <- selected_gwas_multi()
        req(length(gwas_vec) >= 1)
        # Suppress the picker for single-GWAS bundles — nothing to pick
        # between, and the widget is just noise.
        if (length(gwas_vec) < 2L) return(NULL)
        cur <- isolate(input[[input_id]])
        keep <- if (!is.null(cur) && cur %in% gwas_vec) cur else gwas_vec[1L]
        choices <- gwas_choices(gwas_data(), gwas_vec)
        widest <- max(nchar(names(choices)), 12L, na.rm = TRUE)
        w_px   <- max(320L, as.integer(round(widest * 8.5)) + 60L)
        div(style = sprintf("width: %dpx; margin-bottom: 10px;", w_px),
          selectInput(session$ns(input_id), label,
                       choices = choices, selected = keep,
                       multiple = FALSE, width = "100%"))
      })
    }
    if (!is.null(selected_gwas_multi)) {
      .render_per_gwas_picker("maf_gwas_ui",   "maf_gwas")
      .render_per_gwas_picker("qq_gwas_ui",    "qq_gwas")
      .render_per_gwas_picker("rawmh_gwas_ui", "rawmh_gwas")
      .render_per_gwas_picker("log_gwas_ui",   "log_gwas")
    }

    # Plain-text explanations for the QC metrics, used by the legend under the table.
    qc_help <- list(
      "N variants pre-QC" = "Number of variants in the input GWAS before quality control.",
      "Identified genome build" = paste0(
        "Genome build (e.g. GRCh37/GRCh38) inferred from the variant positions; ",
        "shown as 'Unknown' if it could not be determined."),
      "N variants post-QC" = "Number of variants remaining after quality-control filtering.",
      "Lambda GC" = paste0(
        "Genomic inflation factor: median chi-square divided by 0.455. About 1.0 means the ",
        "test statistics are well calibrated; much above 1.1 can indicate confounding (or, ",
        "for a well-powered polygenic trait, genuine widespread signal)."),
      "Max. chi^2" = paste0(
        "The strongest single-variant association in the study, as a chi-square statistic ",
        "derived from the p-value. Larger values mean a stronger top signal."),
      "N genome-wide significant variants" = paste0(
        "Number of variants reaching genome-wide significance (p < 5e-8) after QC. ",
        "Zero is common for underpowered GWAS."),
      "LDSC intercept (SE)" = paste0(
        "LD-score regression intercept. About 1.0 indicates little confounding or ",
        "sample overlap; values above 1 suggest residual confounding or population ",
        "stratification rather than polygenic signal (contrast with Lambda GC). ",
        "Standard error in brackets.")
    )

    # Which GWAS to summarise: the scalar selected GWAS in single mode,
    # the full multi-select vector in compare mode.
    qc_gwas_vec <- reactive({
      req(gwas_data())
      if (!is.null(comparison_mode) && isTRUE(comparison_mode())) {
        selected_gwas_multi()
      } else {
        selected_gwas()
      }
    })

    # Parameter rows shown in the QC Summary table, in display order.
    # (LDSC SNP-h² has moved to the SNP-h² & rG tab; the LDSC intercept
    # stays here as a QC-relevant metric — values close to 1 imply
    # polygenic signal, > 1 suggests confounding.)
    .qc_row_labels <- c(
      "GWAS Label"                         = "label",
      "N variants pre-QC"                  = "n_var_orig",
      "Identified genome build"            = "build",
      "N variants post-QC"                 = "n_snp_final",
      "Lambda GC"                          = "lambda_gc",
      "Max. chi^2"                         = "max_chi2",
      "N genome-wide significant variants" = "n_sig_snp",
      "LDSC intercept (SE)"                = "int_paren"
    )

    .fmt_qc_cell <- function(x, kind) {
      if (all(is.na(x))) return("—")
      switch(kind,
        n_var_orig  = ,
        n_snp_final = ,
        n_sig_snp   = ifelse(is.na(x), "—", formatC(as.integer(x), format = "d", big.mark = ",")),
        lambda_gc   = ,
        max_chi2    = ifelse(is.na(x), "—", sprintf("%.3f", as.numeric(x))),
        # label / build / int_paren / anything else: as-is with NA -> em-dash
        ifelse(is.na(x) | !nzchar(as.character(x)), "—", as.character(x))
      )
    }

    # Wide table (single row per GWAS from the data-layer helper) reshaped
    # into the display shape: rows are parameters, columns are GWAS.
    qc_val <- reactive({
      req(gwas_data(), config_flags())
      req(length(qc_gwas_vec()) >= 1)
      cf <- config_flags()

      long <- build_qc_summary_long(gd = gwas_data(),
                                      gwas_vec = qc_gwas_vec(),
                                      gwas_list = gwas_list())
      if (nrow(long) == 0) return(NULL)

      long <- long[match(qc_gwas_vec(), gwas)]

      # LDSC intercept lives on the QC tab as a QC-relevant stat
      # (values > 1 suggest confounding vs. polygenic signal). Compose
      # the "estimate (SE)" string once so the display renders directly.
      long[, int_paren := ifelse(is.na(int_est), NA_character_,
                                   sprintf("%.3f (%.3f)", int_est, int_se))]

      # Drop the intercept row if LDSC wasn't run — no need to show a
      # column of em-dashes.
      row_labels <- .qc_row_labels
      if (!isTRUE(cf$ldsc)) {
        row_labels <- row_labels[row_labels != "int_paren"]
      }

      cols <- lapply(seq_len(nrow(long)), function(i) {
        vapply(names(row_labels), function(lab) {
          field <- row_labels[[lab]]
          .fmt_qc_cell(long[[field]][i], field)
        }, character(1))
      })
      out <- data.frame(
        Parameter = names(row_labels),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      for (i in seq_along(cols)) out[[long$gwas[i]]] <- cols[[i]]
      out
    })

    # Small static summary — plain shiny::renderTable (avoids DT client-side
    # rendering quirks seen with this table). Parameter column is bold.
    output$qc_table <- renderTable({
      dat <- qc_val()
      req(dat)
      dat$Parameter <- paste0("<strong>", dat$Parameter, "</strong>")
      # colnames only meaningful in compare mode (single mode = "GWAS_NAME"
      # column header is redundant with the label row). Show colnames in
      # both modes for consistency.
      dat
    }, sanitize.text.function = function(x) x,
       colnames = TRUE,
       align = "l",
       rownames = FALSE)

    output$qc_legend <- renderUI({
      req(config_flags())
      cf <- config_flags()
      items <- list(
        "Lambda GC" = qc_help[["Lambda GC"]],
        "Max. chi^2" = qc_help[["Max. chi^2"]],
        "N genome-wide significant variants" = qc_help[["N genome-wide significant variants"]]
      )
      if (isTRUE(cf$ldsc)) {
        items[["LDSC intercept (SE)"]] <- qc_help[["LDSC intercept (SE)"]]
      }
      gd_legend(items)
    })

    # Falls back to selected_gwas() when the per-tab picker hasn't been
    # populated yet (e.g. first render before the observe fires).
    .pick_gwas <- function(input_id) {
      v <- input[[input_id]]
      if (is.null(v) || !nzchar(v)) selected_gwas() else v
    }

    # MAF plot rendering
    output$maf_plot_ui <- renderUI({
      req(gwas_data())
      g <- .pick_gwas("maf_gwas")
      req(g)
      b64 <- gd_read(gwas_data(), g, "gwas_qc")$maf_plot_base64

      if (!is.null(b64)) {
        tagList(
          tags$img(
            src = paste0("data:image/png;base64,", b64),
            style = "max-height: 430px; width: auto;"
          ),
          tags$p(
            style = "font-size: 0.85em; color: #6c757d; margin-top: 8px; max-width: 800px;",
            "This is a quality-control check. It compares the allele frequency reported in your ",
            "GWAS summary statistics against the allele frequency of the same variant in the ",
            "reference panel (1000 Genomes). Large disagreements usually indicate allele ",
            "mis-coding, strand issues, an ancestry mismatch between your GWAS and the reference, ",
            "or data errors, so these variants are removed before downstream analysis."
          ),
          gd_legend(list(
            "X-axis" = "Reference-panel allele frequency.",
            "Y-axis" = "Allele frequency reported in your GWAS summary statistics.",
            "Diagonal line" = "y = x: points on this line agree between the GWAS and the reference.",
            "Points shown" = paste0(
              "Only variants whose two frequencies differ by more than 0.2 are plotted — these ",
              "are the discordant variants removed during QC. An empty or nearly-empty plot is a ",
              "good sign (few or no frequency mismatches).")
          ), heading = "How to read this plot")
        )
      } else {
        div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          icon("chart-bar", style = "font-size: 2em;"),
          br(), br(),
          tags$strong("MAF Plot Unavailable"),
          br(),
          "Allele frequency data was not provided in the input GWAS."
        )
      }
    })

    # QQ plot rendering
    output$qq_plot_ui <- renderUI({
      req(gwas_data())
      g <- .pick_gwas("qq_gwas")
      req(g)
      b64 <- gd_read(gwas_data(), g, "gwas_qc")$qq_plot_base64

      if (!is.null(b64)) {
        tagList(
          tags$img(
            src = paste0("data:image/png;base64,", b64),
            style = "max-height: 500px; width: auto;"
          ),
          gd_legend(list(
            "Diagonal line" = "Expected p-value distribution under the null hypothesis of no association.",
            "Points on the line" = "Test statistics are well calibrated (no inflation).",
            "Early upward departure" = paste0(
              "Points lifting above the line across most of the range suggests inflation or ",
              "confounding (compare with Lambda GC on the QC Summary tab)."),
            "Upward tail at the far right only" = paste0(
              "The strongest associations rising above the line, while the rest follows it, ",
              "is the expected signature of true genetic signal.")
          ), heading = "How to read this plot")
        )
      } else {
        div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          icon("chart-line", style = "font-size: 2em;"),
          br(), br(),
          tags$strong("QQ Plot Unavailable"),
          br(),
          "This results package was produced before the QQ plot rule was added."
        )
      }
    })

    # Cleaner log rendering
    output$cleaner_log_ui <- renderUI({
      req(gwas_data())
      g <- .pick_gwas("log_gwas")
      req(g)
      log_lines <- gd_read(gwas_data(), g, "gwas_qc")$cleaner_dat$log

      if (is.null(log_lines) || length(log_lines) == 0) {
        return(p("Log not available."))
      }

      pre(
        style = "background-color: #1e1e1e; color: #00ff00; font-family: 'Courier New', monospace; font-size: 12px; padding: 15px; border-radius: 4px; white-space: pre-wrap;",
        paste(log_lines, collapse = "\n")
      )
    })

    # ------------------------------------------------------------------
    # Removed Variants: unfiltered (raw-input) Manhattan + lost-signal summary
    # ------------------------------------------------------------------
    # The whole gwas_qc block for the picked GWAS; manhattan_raw / lost_signals
    # are NULL for bundles produced before this feature (handled below).
    rawmh_qc <- reactive({
      req(gwas_data())
      g <- .pick_gwas("rawmh_gwas"); req(g)
      gd_read(gwas_data(), g, "gwas_qc")
    })

    output$raw_manhattan <- renderPlot({
      qc <- rawmh_qc(); req(qc)
      p <- build_raw_manhattan_plot(
        qc$manhattan_raw, qc$lost_signals,
        removed_p_thresh = qc$lost_signals$removed_p_thresh %||% 1e-5)
      # NB: qualify shiny::validate/need — the app also loads jsonlite, whose
      # validate() masks shiny's and errors on a `need` object.
      shiny::validate(shiny::need(!is.null(p), "No plottable raw Manhattan data."))
      p
    })

    output$removed_counts_tbl <- renderTable({
      qc <- rawmh_qc(); req(qc)
      tbl <- gd_build_removed_counts_table(qc$cleaner_dat$log, qc$lost_signals)
      req(!is.null(tbl))
      # bold the Total row
      tbl$Reason[nrow(tbl)] <- paste0("<strong>", tbl$Reason[nrow(tbl)], "</strong>")
      tbl
    }, sanitize.text.function = function(x) x,
       colnames = TRUE, align = "lrr", rownames = FALSE)

    output$lost_loci <- DT::renderDT({
      loci <- rawmh_qc()$lost_signals$loci
      req(!is.null(loci) && nrow(loci) > 0)
      d <- as.data.frame(loci)
      d$reason      <- gd_reason_label(d$reason)
      d$lead_P      <- formatC(d$lead_P, format = "e", digits = 2)
      d$lost_signal <- ifelse(isTRUE_vec(d$lost_signal), "Yes", "No")
      names(d) <- c("CHR", "Lead SNP", "Lead BP", "Lead P", "Reason",
                    "# removed", "# retained (±500kb)",
                    "# retained GW-sig", "Signal lost?")
      dt <- DT::datatable(
        d, rownames = FALSE, escape = TRUE,
        options = list(dom = "tp", paging = TRUE, pageLength = 10,
                       order = list(list(3, "asc"))))
      DT::formatStyle(dt, "Signal lost?", target = "row",
        backgroundColor = DT::styleEqual("Yes", "#fdecea"))
    })

    output$raw_manhattan_ui <- renderUI({
      req(gwas_data())
      g <- .pick_gwas("rawmh_gwas"); req(g)
      qc <- gd_read(gwas_data(), g, "gwas_qc")
      mr <- qc$manhattan_raw
      ls <- qc$lost_signals

      if (is.null(mr)) {
        return(div(
          style = "background-color: #e9ecef; border-radius: 8px; padding: 60px 20px; text-align: center; color: #6c757d;",
          icon("chart-area", style = "font-size: 2em;"),
          br(), br(),
          tags$strong("Unfiltered Manhattan Unavailable"),
          br(),
          "This results package was produced before removed-variant reporting was added, ",
          "or the input GWAS had no chromosome/position columns."
        ))
      }

      gw     <- ls$gw_threshold %||% 5e-8
      rpt    <- ls$removed_p_thresh %||% 1e-5
      n_sig  <- ls$n_removed_sig %||% 0
      n_lost <- ls$n_lost_signal_loci %||% 0
      n_loci <- if (!is.null(ls$loci)) nrow(ls$loci) else 0

      headline <- if (n_sig == 0) {
        div(style = "padding: 10px 14px; border-radius: 6px; background: #e6f4ea; color: #1e4620; margin-bottom: 12px; max-width: 1000px;",
          icon("circle-check"),
          sprintf(" No genome-wide significant variants (P < %s) were removed during QC — no lost signals.",
                  format(gw, scientific = TRUE)))
      } else {
        div(style = "padding: 10px 14px; border-radius: 6px; background: #fdecea; color: #611a15; margin-bottom: 12px; max-width: 1000px;",
          icon("triangle-exclamation"),
          sprintf(" %s genome-wide significant variant(s) were removed during QC across %s approximate loci; %s locus/loci lost their entire genome-wide significant signal.",
                  format(n_sig, big.mark = ","), n_loci, n_lost))
      }

      tagList(
        p(style = "color: var(--gd-text-mute); max-width: 1000px;",
          "Manhattan plot of the ", tags$strong("raw input GWAS"),
          " (before reference harmonisation), showing variants with P < 1e-3. Grey points are ",
          "variants carried through; coloured points were removed during QC/harmonisation, coloured ",
          "by reason. This reveals associations that were dropped before downstream analysis — the ",
          "filtered Manhattan on the SNP Associations tab only shows variants kept in the reference panel."),
        headline,
        div(style = "max-width: 1200px;", plotOutput(ns("raw_manhattan"), height = "480px")),
        div(style = "margin: 6px 0 14px;",
            downloadButton(ns("rawmh_dl"), "Download plot (PNG)", class = "btn-sm")),
        gd_legend(list(
          "Grey points"      = "Variants with P < 1e-3 in the raw input GWAS (before reference harmonisation) that passed QC, shaded alternately by chromosome.",
          "Coloured points"  = "Variants removed during QC/harmonisation, coloured by the reason they were dropped.",
          "Red dashed line"  = sprintf("Genome-wide significance (P = %s).", format(gw, scientific = TRUE)),
          "Blue dotted line" = sprintf("Reporting threshold (P = %s): removed variants below this are recorded individually.", format(rpt, scientific = TRUE)),
          "“Signal lost”"      = "A ~1 Mb locus where every genome-wide significant variant was removed and no retained variant there still reaches significance."
        ), heading = "How to read this plot"),
        h5("Variants removed during QC, by reason"),
        p(style = "color: var(--gd-text-mute); max-width: 1000px;",
          "Number of variants removed at each QC step (across all SNPs), and the subset ",
          "that were genome-wide significant (P < 5e-8). Totals are from the cleaner log; ",
          "pre-reference losses (variants absent from the panel and non-SNPs) are grouped together."),
        div(style = "max-width: 640px;", tableOutput(ns("removed_counts_tbl"))),
        if (n_sig > 0) tagList(
          h5("Approximate affected loci"),
          div(style = "max-width: 1100px;", DT::DTOutput(ns("lost_loci")))
        ) else NULL
      )
    })

    output$rawmh_dl <- downloadHandler(
      filename = function() {
        g <- .pick_gwas("rawmh_gwas")
        sprintf("%s_unfiltered_manhattan.png", g %||% "gwas")
      },
      content = function(file) {
        qc <- rawmh_qc()
        p <- build_raw_manhattan_plot(
          qc$manhattan_raw, qc$lost_signals,
          removed_p_thresh = qc$lost_signals$removed_p_thresh %||% 1e-5)
        if (is.null(p)) {
          grDevices::png(file, width = 4, height = 1, units = "in", res = 96)
          grid::grid.text("No raw Manhattan data in this bundle.")
          grDevices::dev.off()
          return(invisible())
        }
        open_plot_device(file, fmt = "png", w = 12, h = 5, units = "in", dpi = 300)
        print(p)
        grDevices::dev.off()
      }
    )

  })
}

# Vectorised, NA-safe isTRUE for a logical column (lost_signal may contain NA).
isTRUE_vec <- function(x) !is.na(x) & x
