# Overview tab — cross-GWAS summary
#
# Two DT tables:
#   1. QC row: N, lambda_GC, h2 (SE), N sig SNPs, N SNPs final
#   2. Yield row: sig counts across selected sig basis + threshold
#
# Visible only when >= 2 GWAS are selected (gated at the app.R hide/show layer).

overviewUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Overview",
    br(),
    tags$details(class = "gd-details",
      tags$summary("Filter data"),
      tags$div(class = "gd-details-body",
        fluidRow(
          column(4,
            radioButtons(ns("sig_basis"), "Significance basis (yield counts):",
                          choices = c("FDR" = "fdr", "P" = "p"),
                          selected = "fdr", inline = TRUE)
          ),
          column(4,
            numericInput(ns("sig_threshold"), "Significance threshold:",
                          value = 0.05, min = 1e-12, max = 1, step = 0.01)
          ),
          column(4,
            selectInput(ns("gwas_sort"), "Order GWAS by:",
                         choices = c("As selected"   = "as_selected",
                                      "Alphabetical" = "alphabetical",
                                      "Sample size"  = "n",
                                      "N sig SNPs"   = "n_sig_snp",
                                      "SNP-h²"       = "h2"),
                         selected = "as_selected")
          )
        )
      )
    ),
    br(),
    h4("Per-GWAS QC and power"),
    p(
      style = "color: var(--gd-text-mute);",
      "Sample size, genomic inflation (λ", tags$sub("GC"),
      "), SNP-h² (LDSC observed-scale), and genome-wide significant SNP count for each selected GWAS."
    ),
    DT::DTOutput(ns("qc_tbl")),
    br(),
    h4("Yield: counts of significant entities per GWAS"),
    p(
      style = "color: var(--gd-text-mute);",
      "Counts respect the significance basis and threshold above. ",
      "ATC counts use best-per-cell across TWAS-GSEA panels."
    ),
    DT::DTOutput(ns("yield_tbl")),
    br(),
    gd_legend(list(
      "N" = "Approximate per-GWAS sample size, taken as the maximum N across clumped lead SNPs.",
      "λ GC" = "Genomic inflation factor from the sumstat cleaner log.",
      "h² (SE)" = "LDSC observed-scale SNP heritability with its standard error. Blank when LDSC was not run.",
      "N sig SNPs" = "Count of SNPs passing p < 5×10⁻⁸ in the cleaned sumstats.",
      "Loci" = "Number of independent clumped loci.",
      "Sig Genes" = "MAGMA gene-level significant entities under the chosen basis + threshold.",
      "Sig Tissues" = "GTEx MAGMA tissue-specific significant tissues.",
      "Sig ATC" = "Union of significant ATC classes across MAGMA and TWAS-GSEA (best-per-cell)."
    ))
  )
}

overviewServer <- function(id, gwas_data, selected_gwas_multi,
                            comparison_mode, comparison_long) {
  moduleServer(id, function(input, output, session) {

    qc_dat <- reactive({
      req(gwas_data(), comparison_mode())
      build_overview_qc(gwas_data(), selected_gwas_multi())
    })

    ordered_gwas <- reactive({
      req(comparison_mode())
      sort_mode <- if (is.null(input$gwas_sort)) "as_selected" else input$gwas_sort
      order_gwas(selected_gwas_multi(), sort_mode, qc_dat())
    })

    output$qc_tbl <- DT::renderDT({
      req(comparison_mode())
      d <- qc_dat()
      ord <- ordered_gwas()
      d <- d[match(ord, gwas)]
      out <- data.frame(
        GWAS        = d$gwas,
        N           = formatC(d$N, format = "d", big.mark = ","),
        `lambda GC` = ifelse(is.na(d$lambda_gc), "—", sprintf("%.3f", d$lambda_gc)),
        `h2 (SE)`   = ifelse(is.na(d$h2), "—",
                             sprintf("%.3f (%.3f)", d$h2, d$h2_se)),
        `N sig SNPs`  = ifelse(is.na(d$n_sig_snp), "—",
                                formatC(d$n_sig_snp, format = "d", big.mark = ",")),
        `N SNPs final` = ifelse(is.na(d$n_snp_final), "—",
                                 formatC(d$n_snp_final, format = "d", big.mark = ",")),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      # NA sample-size shows as "NA"; hide it.
      out$N[out$N == "NA"] <- "—"
      DT::datatable(
        out, rownames = FALSE, escape = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE)
      )
    })

    yield_dat <- reactive({
      req(comparison_mode(), comparison_long())
      build_overview_yield(
        comparison_long(), gwas_data(), selected_gwas_multi(),
        sig_basis     = if (is.null(input$sig_basis)) "fdr" else input$sig_basis,
        sig_threshold = if (is.null(input$sig_threshold)) 0.05 else as.numeric(input$sig_threshold)
      )
    })

    output$yield_tbl <- DT::renderDT({
      req(comparison_mode())
      d <- yield_dat()
      ord <- ordered_gwas()
      d <- d[match(ord, gwas)]
      fmt <- function(v) ifelse(is.na(v), "—", formatC(v, format = "d", big.mark = ","))
      out <- data.frame(
        GWAS          = d$gwas,
        Loci          = fmt(d$n_loci),
        `Sig Genes`   = fmt(d$n_genes_sig),
        `Sig Tissues` = fmt(d$n_tissues_sig),
        `Sig ATC`     = fmt(d$n_atc_sig),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(
        out, rownames = FALSE, escape = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE)
      )
    })
  })
}
