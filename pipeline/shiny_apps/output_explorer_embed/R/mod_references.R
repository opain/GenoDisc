referencesUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title="References",
    br(),
    # Draft write-up downloads (auto-generated per run; absent in older bundles).
    uiOutput(ns("writeup_downloads")),
    p("Please be sure to cite the software and datasets used by this pipeline. Relevent citations are shown below:"),
    fluidRow(
      column(width=8,
             tableOutput(ns("reference_table")),
      )
    ),
  )
}

referencesServer <- function(id, gwas_data = NULL) {
  moduleServer(id, function(input, output, session) {

    # ---- draft write-up downloads ------------------------------------------
    # The pipeline writes a draft methods/results section + reference list into
    # the bundle at package/writeup/. Surface them as downloads here. gwas_data()
    # is the opened gd handle; gd$root_dir is the package directory. Bundles
    # produced before this feature have no writeup/ dir -> show a muted note.
    writeup_dir <- reactive({
      if (is.null(gwas_data)) return(NULL)
      gd <- tryCatch(gwas_data(), error = function(e) NULL)
      if (is.null(gd) || is.null(gd$root_dir) || is.na(gd$root_dir)) return(NULL)
      d <- file.path(gd$root_dir, "writeup")
      if (dir.exists(d)) d else NULL
    })

    output$writeup_downloads <- renderUI({
      ns <- session$ns
      d <- writeup_dir()
      if (is.null(d)) return(NULL)
      tagList(
        tags$div(
          style = "border:1px solid #d0d7de; border-radius:6px; padding:12px 14px; margin-bottom:16px;",
          tags$strong("Draft write-up (automatically generated)"),
          tags$p(style = "margin:6px 0 10px; color:#57606a;",
                 "A draft methods section, descriptive results section and reference list generated from this run. ",
                 tags$em("These are drafts — check and edit before use.")),
          downloadButton(ns("dl_methods"),  "Methods (.md)",     class = "btn-sm"),
          downloadButton(ns("dl_bib"),      "References (.bib)", class = "btn-sm"),
          downloadButton(ns("dl_all"),      "All files (.tar.gz)", class = "btn-sm")
        )
      )
    })

    output$dl_methods <- downloadHandler(
      filename = function() "genodisc_draft_methods.md",
      content  = function(file) {
        d <- writeup_dir(); req(d)
        file.copy(file.path(d, "methods.md"), file, overwrite = TRUE)
      }
    )
    output$dl_bib <- downloadHandler(
      filename = function() "genodisc_references.bib",
      content  = function(file) {
        d <- writeup_dir(); req(d)
        file.copy(file.path(d, "references.bib"), file, overwrite = TRUE)
      }
    )
    output$dl_all <- downloadHandler(
      filename = function() "genodisc_draft_writeup.tar.gz",
      content  = function(file) {
        d <- writeup_dir(); req(d)
        old <- setwd(dirname(d)); on.exit(setwd(old), add = TRUE)
        utils::tar(file, basename(d), compression = "gzip")
      }
    )


    # Small, static reference table - plain shiny::renderTable instead of a DT
    # widget, since this needs no sorting/searching/pagination (avoids the
    # DT/DataTables client-side rendering issues seen with this table).
    output$reference_table<-renderTable({
      ref_table<-NULL

      ref_table<-rbind(ref_table, data.frame(Name = 'LD Score Regression',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/25642630">25642630</a>',
                                             Website = '<a href="https://github.com/bulik/ldsc">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate SNP-h2 and LD-score intercept'))

      ref_table<-rbind(ref_table, data.frame(Name = 'SBayesRC',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/38689000">38689000</a>',
                                             Website = '<a href="https://github.com/zhilizheng/SBayesRC">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate SNP-h2, polygenicity and functional-annotation heritability enrichment'))

      ref_table<-rbind(ref_table, data.frame(Name = 'SumHer (LDAK)',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/30510236">30510236</a>',
                                             Website = '<a href="https://dougspeed.com/">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate SNP-h2 and functional-annotation heritability enrichment (BLD-LDAK model)'))

      ref_table<-rbind(ref_table, data.frame(Name = 'Manc-COJO',
                                             PMID = '<a href="https://doi.org/10.1038/s41588-026-02799-1">10.1038/s41588-026-02799-1</a>',
                                             Website = '<a href="https://github.com/light156/multi-ancestry-COJO">link</a>',
                                             Type = 'Software',
                                             Use = 'Conditional & joint analysis for independent association signals'))

      ref_table<-rbind(ref_table, data.frame(Name = 'SuSiE',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/35853082">35853082</a>',
                                             Website = '<a href="https://stephenslab.github.io/susie-paper">link</a>',
                                             Type = 'Software',
                                             Use = 'Variant-level finemapping'))

      ref_table<-rbind(ref_table, data.frame(Name = 'FUSION',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/26854917">26854917</a>',
                                             Website = '<a href="http://gusevlab.org/projects/fusion/">link</a>',
                                             Type = 'Software',
                                             Use = 'Perform TWAS/PWAS'))

      ref_table<-rbind(ref_table, data.frame(Name = 'COLOC',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/24830394">24830394</a>',
                                             Website = '<a href="https://chr1swallace.github.io/coloc/">link</a>',
                                             Type = 'Software',
                                             Use = 'Used for colocalisation within FUSION TWAS/PWAS analysis'))

      ref_table<-rbind(ref_table, data.frame(Name = 'SMR',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/29763751">29763751</a>',
                                             Website = '<a href="https://yanglab.westlake.edu.cn/software/smr/#Overview">link</a>',
                                             Type = 'Software',
                                             Use = 'Infer differential expression or protein levels'))

      ref_table<-rbind(ref_table, data.frame(Name = 'MAGMA',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/25885710">25885710</a>',
                                             Website = '<a href="https://ctg.cncr.nl/software/magma">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate gene associations and drug enrichment using DrugTargetor database'))

      ref_table<-rbind(ref_table, data.frame(Name = 'TWAS-GSEA',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/31230729">31230729</a>',
                                             Website = '<a href="https://github.com/opain/TWAS-GSEA">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate drug enrichment using TWAS and DrugTargetor'))

      ref_table<-rbind(ref_table, data.frame(Name = 'PsychENCODE DLPFC TWAS SNP-weights',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/30545856">30545856</a>',
                                             Website = '<a href="http://resource.psychencode.org/">link</a>',
                                             Type = 'Dataset',
                                             Use = 'TWAS'))

      ref_table<-rbind(ref_table, data.frame(Name = 'PsychENCODE DLPFC eQTL data for SMR',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/30545857">30545857</a>',
                                             Website = '<a href="http://cnsgenomics.com/data/SMR/">link</a>',
                                             Type = 'Dataset',
                                             Use = 'SMR'))

      ref_table<-rbind(ref_table, data.frame(Name = 'MetaBrain eQTL data for SMR',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/36823318">36823318</a>',
                                             Website = '<a href="https://www.metabrain.nl/cis-eqtls.html">link</a>',
                                             Type = 'Dataset',
                                             Use = 'SMR'))

      ref_table<-rbind(ref_table, data.frame(Name = 'ROSMAP/Banner DLPFC PWAS SNP-weights and eQTL in SMR results',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/33510477">33510477</a>',
                                             Website = '<a href="https://www.synapse.org/#!Synapse:syn23627957">link</a>',
                                             Type = 'Dataset',
                                             Use = 'PWAS and SMR'))

      ref_table<-rbind(ref_table, data.frame(Name = 'DrugTargetor',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/30517594">30517594</a>',
                                             Website = '<a href="https://drugtargetor.com/">link</a>',
                                             Type = 'Dataset',
                                             Use = 'Drug repurposing'))

      ref_table
    }, sanitize.text.function = function(x) x, rownames = FALSE)
  })
}
