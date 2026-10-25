# Effector-gene prioritisation tab: PoPS gene scores, CALDERA per-credible-set
# causal-gene probabilities, and credible-set coding annotation (AlphaMissense +
# consequence calls). Reads the packaged blocks mol_assoc/pops, caldera and
# credset_coding (the last is a list of $summary + $variants). Pattern mirrors
# R/mod_mol_assoc.R (DT tables + gd_legend); null-safe per sub-tab.

effectorUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Effector genes",
    br(),
    p("Locus-level effector-gene prioritisation and credible-set coding annotation ",
      "for the selected GWAS."),
    hr(),
    tabsetPanel(
      id = ns("effector_tabs"),
      tabPanel(
        title = "PoPS", value = "pops",
        br(), p("Polygenic Priority Score (PoPS) per gene — genome-wide gene ",
                "prioritisation from MAGMA gene associations + functional features. ",
                "Higher scores indicate stronger prioritisation."),
        DT::dataTableOutput(ns("pops_table")),
        uiOutput(ns("pops_msg")), br()
      ),
      tabPanel(
        title = "CALDERA", value = "caldera",
        br(), p("CALDERA causal-gene probabilities per credible set (one row per gene ",
                "in each locus). ", strong("caldera"), " is the per-locus probability ",
                "(sums to 1 within a credible set); ", strong("rank"), " orders genes ",
                "within a locus; ", strong("nearest_gene"), " flags the closest gene as a ",
                "transparency baseline."),
        DT::dataTableOutput(ns("caldera_table")),
        uiOutput(ns("caldera_msg")),
        gd_legend(list(
          "caldera"      = "CALDERA causal-gene probability for this gene in this credible set (0–1, sums to 1 per locus).",
          "rank"         = "Rank of this gene within the credible set by CALDERA probability (1 = top).",
          "nearest_gene" = "TRUE if this is the gene nearest the lead variant (baseline comparator).",
          "pops"         = "PoPS score feature; dist = bp distance to the lead variant; coding = summed non-synonymous credible-set PIP."
        )), br()
      ),
      tabPanel(
        title = "Credible-set coding", value = "credset_coding",
        br(), p("Per-credible-set coding summary, then the per-variant annotation. ",
                strong("flag_protein_altering"), " marks sets carrying a meaningful share ",
                "of PIP on protein-altering variants."),
        h4("Per credible set"),
        DT::dataTableOutput(ns("coding_summary_table")),
        br(), h4("Per variant"),
        DT::dataTableOutput(ns("coding_variants_table")),
        uiOutput(ns("coding_msg")),
        gd_legend(list(
          "pip_protein_altering" = "Summed PIP on protein-altering (missense/stop/frameshift/splice/…) variants in the credible set.",
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

    # Read a block for the currently-selected GWAS; NULL on any problem.
    read_block <- function(block) {
      gd <- tryCatch(gwas_data(), error = function(e) NULL)
      g  <- tryCatch(selected_gwas(), error = function(e) NULL)
      if (is.null(gd) || is.null(g) || !nzchar(g)) return(NULL)
      tryCatch(gd_read(gd, g, block), error = function(e) NULL)
    }

    empty_dt <- function() DT::datatable(data.frame(), rownames = FALSE,
                                         options = list(dom = "t"))
    msg <- function(text) renderUI(div(class = "gd-legend", em(text)))
    show_dt <- function(df, sci_cols = NULL) {
      df <- as.data.frame(df); df[is.na(df)] <- "NA"
      DT::datatable(df, rownames = FALSE, filter = "top",
                    options = list(pageLength = 15, scrollX = TRUE,
                                   columnDefs = list(list(className = "dt-center",
                                                          targets = seq_len(ncol(df)) - 1))))
    }

    # ---- PoPS ----
    output$pops_table <- DT::renderDataTable({
      d <- read_block("mol_assoc/pops")
      if (is.null(d) || nrow(d) == 0) return(empty_dt())
      show_dt(d)
    })
    output$pops_msg <- renderUI({
      d <- read_block("mol_assoc/pops")
      if (is.null(d) || nrow(d) == 0)
        return(div(class = "gd-legend", em("No PoPS results for this GWAS.")))
      NULL
    })

    # ---- CALDERA ----
    output$caldera_table <- DT::renderDataTable({
      d <- read_block("caldera")
      if (is.null(d) || nrow(d) == 0) return(empty_dt())
      show_dt(d)
    })
    output$caldera_msg <- renderUI({
      d <- read_block("caldera")
      if (is.null(d) || nrow(d) == 0)
        return(div(class = "gd-legend", em("No CALDERA results for this GWAS."))); NULL
    })

    # ---- Credible-set coding (block is a list: $summary + $variants) ----
    output$coding_summary_table <- DT::renderDataTable({
      d <- read_block("credset_coding")
      if (is.null(d) || is.null(d$summary) || nrow(d$summary) == 0) return(empty_dt())
      show_dt(d$summary)
    })
    output$coding_variants_table <- DT::renderDataTable({
      d <- read_block("credset_coding")
      if (is.null(d) || is.null(d$variants) || nrow(d$variants) == 0) return(empty_dt())
      show_dt(d$variants)
    })
    output$coding_msg <- renderUI({
      d <- read_block("credset_coding")
      if (is.null(d) || (is.null(d$summary) && is.null(d$variants)))
        return(div(class = "gd-legend", em("No credible-set coding annotation for this GWAS."))); NULL
    })
  })
}
