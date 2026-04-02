referencesUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title="References",
    br(),
    p("Please be sure to cite the software and datasets used by this pipeline. Relevent citations are shown below:"),
    fluidRow(
      column(width=8,
             dataTableOutput(ns("reference_table")),
      )
    ),
  )
}

referencesServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    output$reference_table<-renderDataTable({
      ref_table<-NULL

      ref_table<-rbind(ref_table, data.frame(Name = 'LD Score Regression',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/25642630">25642630</a>',
                                             Website = '<a href="https://github.com/bulik/ldsc">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate SNP-h2 and LD-score intercept'))

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

      datatable(ref_table,
                rownames= FALSE,
                escape = -2:-3,
                selection = 'none',
                options = list(
                  paging = FALSE,
                  dom = 'lrt'
                )
                )
    })
  })
}
