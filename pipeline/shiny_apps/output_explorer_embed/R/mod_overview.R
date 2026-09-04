# Overview tab — bundle-level summary (works for single and multi-GWAS bundles)
#
# Two DT tables:
#   1. QC row: N, lambda_GC, h2 (SE), N sig SNPs, N SNPs post-QC
#   2. Yield row: sig counts across selected sig basis + threshold
#
# Visible whenever a bundle is loaded, regardless of how many GWAS are selected.

overviewUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Overview",
    br(),
    tags$details(class = "gd-details",
      tags$summary("Filter data"),
      tags$div(class = "gd-details-body",
        fluidRow(
          column(6,
            numericInput(ns("sig_threshold"), "FDR significance threshold:",
                          value = 0.05, min = 1e-12, max = 1, step = 0.01)
          ),
          column(6,
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
    tags$div(style = "max-width: 1100px;",
      h4("QC and power"),
      p(
        style = "color: var(--gd-text-mute);",
        "Sample size, genomic inflation (λ", tags$sub("GC"),
        "), SNP-h² (LDSC observed-scale), and genome-wide significant SNP count."
      ),
      DT::DTOutput(ns("qc_tbl")),
      br(),
      h4("Yield: counts of significant entities"),
      p(
        style = "color: var(--gd-text-mute);",
        "Counts of high-confidence entities under the FDR threshold above. ",
        "Molecular counts require both FDR-significance and colocalisation (FUSION) or HEIDI support (SMR)."
      ),
      DT::DTOutput(ns("yield_tbl"))
    ),
    br(),
    gd_legend(list(
      "N" = "Approximate sample size, taken as the maximum N across clumped lead SNPs.",
      "λ GC" = "Genomic inflation factor from the sumstat cleaner log.",
      "h² (SE)" = "LDSC observed-scale SNP heritability with its standard error. Blank when LDSC was not run.",
      "N sig SNPs" = "Count of SNPs passing p < 5×10⁻⁸ in the cleaned sumstats.",
      "N SNPs post-QC" = "Number of SNPs remaining after quality-control filtering in the sumstat cleaner.",
      "Loci" = "Number of independent clumped loci.",
      "Sig TWAS" = "TWAS-FUSION genes that are FDR-significant AND colocalised (COLOC-supported), unique across panels.",
      "Sig PWAS" = "PWAS-FUSION genes that are FDR-significant AND colocalised (COLOC-supported), unique across panels.",
      "Sig SMR-eQTL" = "SMR-expression genes that are FDR-significant AND pass HEIDI (p_HEIDI > 0.05), unique across panels.",
      "Sig SMR-pQTL" = "SMR-protein genes that are FDR-significant AND pass HEIDI (p_HEIDI > 0.05), unique across panels.",
      "Sig SuSiE" = "Unique genes containing a full 95% credible set reported by SuSiE fine-mapping.",
      "Sig Tissues" = "GTEx MAGMA tissue-specific FDR-significant tissues.",
      "Sig ATC" = "Union of FDR-significant ATC classes across MAGMA and TWAS-GSEA (best-per-cell).",
      "Sig Drugs" = "Union of FDR-significant drugs across MAGMA and TWAS-GSEA (best-per-cell)."
    ))
  )
}

overviewServer <- function(id, gwas_data, selected_gwas_multi,
                            comparison_mode, comparison_long) {
  moduleServer(id, function(input, output, session) {

    qc_dat <- reactive({
      req(gwas_data())
      build_overview_qc(gwas_data(), selected_gwas_multi())
    })

    ordered_gwas <- reactive({
      req(gwas_data())
      sort_mode <- if (is.null(input$gwas_sort)) "as_selected" else input$gwas_sort
      order_gwas(selected_gwas_multi(), sort_mode, qc_dat())
    })

    output$qc_tbl <- DT::renderDT({
      d <- qc_dat()
      ord <- ordered_gwas()
      d <- d[match(ord, gwas)]
      out <- data.frame(
        GWAS        = d$gwas,
        N           = formatC(d$N, format = "d", big.mark = ","),
        `lambda GC` = ifelse(is.na(d$lambda_gc), "—", sprintf("%.3f", d$lambda_gc)),
        `h2 (SE)`   = ifelse(is.na(d$h2), "—",
                             sprintf("%.3f (%.3f)", d$h2, d$h2_se)),
        `N sig SNPs`   = ifelse(is.na(d$n_sig_snp), "—",
                                 formatC(d$n_sig_snp, format = "d", big.mark = ",")),
        `N SNPs post-QC` = ifelse(is.na(d$n_snp_final), "—",
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
      req(gwas_data(), comparison_long())
      build_overview_yield(
        comparison_long(), gwas_data(), selected_gwas_multi(),
        sig_threshold = if (is.null(input$sig_threshold)) 0.05 else as.numeric(input$sig_threshold)
      )
    })

    output$yield_tbl <- DT::renderDT({
      d <- yield_dat()
      ord <- ordered_gwas()
      d <- d[match(ord, gwas)]
      fmt <- function(v) ifelse(is.na(v), "—", formatC(v, format = "d", big.mark = ","))
      out <- data.frame(
        GWAS            = d$gwas,
        Loci            = fmt(d$n_loci),
        `Sig TWAS`      = fmt(d$n_twas_hc),
        `Sig PWAS`      = fmt(d$n_pwas_hc),
        `Sig SMR-eQTL`  = fmt(d$n_smr_expr_hc),
        `Sig SMR-pQTL`  = fmt(d$n_smr_prot_hc),
        `Sig SuSiE`     = fmt(d$n_susie_hc),
        `Sig Tissues`   = fmt(d$n_tissues_sig),
        `Sig ATC`       = fmt(d$n_atc_sig),
        `Sig Drugs`     = fmt(d$n_drugs_sig),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(
        out, rownames = FALSE, escape = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE)
      )
    })
  })
}
