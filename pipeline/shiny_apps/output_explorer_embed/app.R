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

# Load functions and modules
source('functions.R')
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)

options(shiny.maxRequestSize = 600 * 1024 * 1024)

# Define UI
ui <- fluidPage(

  shinyjs::useShinyjs(),

  tags$style(HTML("
    .custom-panel {
          padding: 20px;
          border: 1px solid #ccc;
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
    parse_config_flags(shared$gwas_data()$configuration$config)
  })

  gwas_list <- reactive({
    req(shared$gwas_data())
    shared$gwas_data()$configuration$gwas_list
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
    toggle("Enrichment Analysis", any(cf$magma_drugtargetor, cf$gcsc, cf$twas_gsea_drugtargetor))
  })

  # Wire up all modules
  gwasQcServer("gwas_qc", shared$gwas_data, shared$selected_gwas, gwas_list)
  snpAssocServer("snp_assoc", shared$gwas_data, shared$selected_gwas)
  molAssocServer("mol_assoc", shared$gwas_data, shared$selected_gwas, config_flags)
  enrichmentServer("enrichment", shared$gwas_data, shared$selected_gwas, config_flags)
  referencesServer("references")
  configurationServer("configuration", shared$gwas_data)
}

# Run the Shiny app
shinyApp(ui, server)
