# Load required libraries
library(shiny)
library(data.table)
library(DT)
library(grid)
library(ggplot2)
library(cowplot)
library(shinythemes)
library(shinycssloaders)
library(R.utils)
library(shinyjs)
library(sass)
library(memoise)
library(jsonlite)

# Load functions and modules
source('functions.R')
source('reader.R')
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)

options(shiny.maxRequestSize = 600 * 1024 * 1024)
options(shiny.legacy.datatable = TRUE)

# Define UI
ui <- fluidPage(

  shinyjs::useShinyjs(),

  tags$style(HTML("
    .custom-panel {
          padding: 20px;
          border: 1px solid #ccc;
    }
    .gd-legend {
          font-size: 0.85em;
          color: #6c757d;
          margin-top: 10px;
          max-width: 800px;
    }
    .gd-legend ul {
          margin: 4px 0 0 0;
          padding-left: 20px;
    }
    .gd-legend li {
          margin-bottom: 3px;
    }
    /* fileInput: give the filename text extra left-padding so long names
       don't creep under the Browse button, and make the upload progress
       bar tall enough that the percentage is readable. */
    .shiny-input-container .input-group .form-control {
          padding-left: 16px;
    }
    .shiny-file-input-progress {
          height: 24px !important;
          margin-top: 6px;
    }
    .shiny-file-input-progress .progress-bar {
          height: 24px !important;
          line-height: 24px;
          font-size: 13px;
    }
  ")),

  theme = shinythemes::shinytheme("paper"),

  tabsetPanel(id = "main_tabs",
    dataInputUI("data_input"),
    gwasQcUI("gwas_qc"),
    snpAssocUI("snp_assoc"),
    molAssocUI("mol_assoc"),
    enrichmentUI("enrichment"),
    referencesUI("references"),
    configurationUI("configuration")
  )
)

# Define server logic
server <- function(input, output, session) {

  # Hide all tabs except Data Input on startup
  observeEvent(session$clientData, {
    for (tab in c("GWAS QC", "SNP Associations", "Molecular Associations",
                   "Enrichment Analysis", "References", "Configuration")) {
      hideTab("main_tabs", tab)
    }
  }, once = TRUE)

  # Data input module returns shared state
  shared <- dataInputServer("data_input")

  # Parse config flags once
  config_flags <- reactive({
    req(shared$gwas_data())
    parse_config_flags(gd_config(shared$gwas_data())$flags_raw)
  })

  gwas_list <- reactive({
    req(shared$gwas_data())
    gd_config(shared$gwas_data())$gwas_list
  })

  # Show relevant tabs when data loads
  observeEvent(config_flags(), {
    cf <- config_flags()
    showTab("main_tabs", "GWAS QC")
    showTab("main_tabs", "References")
    showTab("main_tabs", "Configuration")

    toggle <- function(tab, cond) {
      if (cond) showTab("main_tabs", tab) else hideTab("main_tabs", tab)
    }
    toggle("SNP Associations", any(cf$clump, cf$cojo, cf$finemap))
    toggle("Molecular Associations", cf$mol_assoc)
    toggle("Enrichment Analysis", any(cf$magma_drugtargetor, cf$gcsc, cf$twas_gsea_drugtargetor, cf$twas_gsea_drugtargetor_nondirectional, cf$tissue_magma))
  })

  # Wire up all modules
  gwasQcServer("gwas_qc", shared$gwas_data, shared$selected_gwas, gwas_list, config_flags)
  snpAssocServer("snp_assoc", shared$gwas_data, shared$selected_gwas, config_flags)
  molAssocServer("mol_assoc", shared$gwas_data, shared$selected_gwas, config_flags)
  enrichmentServer("enrichment", shared$gwas_data, shared$selected_gwas, config_flags)
  referencesServer("references")
  configurationServer("configuration", shared$gwas_data)
}

# Run the Shiny app
shinyApp(ui, server)
