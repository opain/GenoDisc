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

  # GenoDisc-web design tokens (light mode) mirrored from the Django app's
  # _design_system.html so the two look congruent.
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = ""),
    tags$link(rel = "stylesheet",
              href = paste0("https://fonts.googleapis.com/css2?",
                            "family=Inter+Tight:wght@400;500;600;700&",
                            "family=Montserrat:wght@500;600&",
                            "display=swap"))
  ),

  tags$style(HTML("
    :root {
      --gd-bg:        #f6f7f9;
      --gd-bg-2:      #eceef2;
      --gd-panel:     #ffffff;
      --gd-panel-2:   #f1f3f7;
      --gd-border:    #dce0e7;
      --gd-border-2:  #c4cad5;
      --gd-text:      #1b2029;
      --gd-text-dim:  #454c5a;
      --gd-text-mute: #626a78;
      --gd-accent:    #0f766e;
      --gd-accent-h:  #0c655e;
      --gd-accent-2:  #a16207;
      --gd-accent-3:  #4f46e5;
      --gd-danger:    #dc2626;
      --gd-shadow-sm: 0 1px 2px rgba(0,0,0,.06);
      --gd-shadow:    0 8px 24px rgba(0,0,0,.10);
      --gd-r-btn: 10px;
      --gd-r-input: 6px;
      --gd-r-card: 16px;
      --gd-font-body: 'Inter Tight', -apple-system, BlinkMacSystemFont,
                     'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
      --gd-font-head: 'Montserrat', -apple-system, BlinkMacSystemFont,
                     'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
    }

    /* ===== Base ===== */
    body, .container-fluid, .container {
      font-family: var(--gd-font-body);
      font-size: 14.5px;
      line-height: 1.5;
      color: var(--gd-text);
      background: var(--gd-bg);
    }
    h1, h2, h3, h4, h5, h6 {
      font-family: var(--gd-font-head);
      font-weight: 600;
      letter-spacing: -0.01em;
      color: var(--gd-text);
    }
    p, li, label { color: var(--gd-text); }
    a, a:visited { color: var(--gd-accent); }
    a:hover, a:focus { color: var(--gd-accent-h); text-decoration: underline; }
    hr { border-top-color: var(--gd-border); }

    /* ===== Buttons ===== */
    .btn, .btn-default,
    .action-button,
    .shiny-download-link {
      background: transparent;
      border: 1px solid var(--gd-border-2);
      color: var(--gd-text);
      border-radius: var(--gd-r-btn);
      padding: 8px 14px;
      font-family: var(--gd-font-body);
      font-weight: 500;
      font-size: 13.5px;
      transition: background .12s ease, border-color .12s ease, color .12s ease;
    }
    .btn:hover, .btn-default:hover,
    .action-button:hover,
    .shiny-download-link:hover {
      background: var(--gd-panel-2);
      border-color: var(--gd-border-2);
      color: var(--gd-text);
    }
    .btn-primary,
    .btn-file,
    .shiny-download-link {
      background: var(--gd-accent);
      border-color: var(--gd-accent);
      color: #fff;
    }
    .btn-primary:hover,
    .btn-file:hover,
    .shiny-download-link:hover {
      background: var(--gd-accent-h);
      border-color: var(--gd-accent-h);
      color: #fff;
    }
    .btn:focus, .btn:active,
    .action-button:focus,
    .shiny-download-link:focus {
      outline: none;
      box-shadow: 0 0 0 3px rgba(15,118,110,.15);
    }

    /* ===== Form inputs ===== */
    .form-control,
    input[type='text'], input[type='number'], input[type='search'], textarea,
    .selectize-input {
      background: var(--gd-bg-2);
      border: 1px solid var(--gd-border-2);
      border-radius: var(--gd-r-input);
      color: var(--gd-text);
      font-family: var(--gd-font-body);
      font-size: 13.5px;
      padding: 8px 12px;
      box-shadow: none;
    }
    .form-control:focus,
    input[type='text']:focus, input[type='number']:focus, textarea:focus,
    .selectize-input.focus {
      border-color: var(--gd-accent);
      box-shadow: 0 0 0 3px rgba(15,118,110,.15);
      outline: none;
    }
    .selectize-dropdown {
      border: 1px solid var(--gd-border);
      border-radius: var(--gd-r-input);
      background: var(--gd-panel);
    }
    .selectize-dropdown .active {
      background: var(--gd-panel-2);
      color: var(--gd-text);
    }
    .selectize-control.multi .selectize-input > div {
      background: rgba(15,118,110,.10);
      color: var(--gd-accent-h);
      border-radius: 4px;
      padding: 2px 8px;
    }
    /* Slider (ionRangeSlider) — teal accent */
    .irs--shiny .irs-bar { background: var(--gd-accent); border-color: var(--gd-accent); }
    .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
      background: var(--gd-accent);
    }
    .irs--shiny .irs-handle {
      border-color: var(--gd-accent);
    }

    /* ===== Panels / cards ===== */
    .well, .sidebar-panel, .sidebarPanel {
      background: var(--gd-panel);
      border: 1px solid var(--gd-border);
      border-radius: var(--gd-r-card);
      box-shadow: var(--gd-shadow-sm);
      padding: 20px;
    }
    .custom-panel {
      background: var(--gd-panel);
      border: 1px solid var(--gd-border);
      border-radius: var(--gd-r-card);
      box-shadow: var(--gd-shadow-sm);
      padding: 20px;
    }

    /* ===== Collapsible <details> (Filter data / Plot options) =====
       Consolidated here so both mod_mol_assoc.R and mod_enrichment.R
       inherit consistent styling. */
    .gd-details > summary {
      cursor: pointer;
      padding: 8px 12px;
      background: var(--gd-panel-2);
      border: 1px solid var(--gd-border);
      border-radius: var(--gd-r-input);
      font-family: var(--gd-font-body);
      font-weight: 600;
      color: var(--gd-text);
      user-select: none;
      list-style: none;
      max-width: 1100px;
    }
    .gd-details > summary::-webkit-details-marker { display: none; }
    .gd-details > summary::before {
      content: '\\25B8';
      display: inline-block;
      width: 1em;
      margin-right: 4px;
      transition: transform 0.15s ease;
      font-size: 1.35em;
      line-height: 1;
      vertical-align: -0.05em;
      color: var(--gd-accent);
    }
    .gd-details[open] > summary::before { transform: rotate(90deg); }
    .gd-details > summary:hover { background: #e6e9ef; }
    .gd-details[open] > summary {
      border-radius: var(--gd-r-input) var(--gd-r-input) 0 0;
      border-bottom: none;
    }
    .gd-details + .gd-details { margin-top: 8px; }
    .gd-details-body {
      padding: 12px 15px;
      background: var(--gd-panel-2);
      border: 1px solid var(--gd-border);
      border-top: none;
      border-radius: 0 0 var(--gd-r-input) var(--gd-r-input);
      max-width: 1100px;
    }
    .gd-details-intro {
      color: var(--gd-text-mute);
      font-size: 0.9em;
      margin-bottom: 12px;
    }
    .gd-details .selectize-control.multi .selectize-input {
      max-height: 120px;
      overflow-y: auto;
    }

    /* ===== Tabs =====
       Bootstrap Paper (underlying theme) marks the active tab with a 2px
       blue bottom-border. Its :hover rule fights ours at equal specificity,
       causing a blue-outline flash + 1px height jump on hover of an active
       tab. !important on our tab rules guarantees no state flicker. */
    .nav-tabs { border-bottom: 1px solid var(--gd-border) !important; }
    .nav-tabs > li > a,
    .nav-tabs > li > a:hover,
    .nav-tabs > li > a:focus,
    .nav-tabs > li > a:active {
      border-radius: 8px 8px 0 0 !important;
      padding: 10px 15px !important;
      font-family: var(--gd-font-body) !important;
      font-weight: 500 !important;
      margin-right: 2px !important;
      transition: none !important;
      outline: none !important;
      box-shadow: none !important;
    }
    .nav-tabs > li > a {
      color: var(--gd-text-dim) !important;
      background: transparent !important;
      border: 1px solid transparent !important;
    }
    .nav-tabs > li > a:hover {
      background: var(--gd-panel-2) !important;
      color: var(--gd-text) !important;
      border-color: var(--gd-border) !important;
      border-bottom-color: transparent !important;
    }
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:hover,
    .nav-tabs > li.active > a:focus {
      color: var(--gd-accent) !important;
      background: var(--gd-panel) !important;
      border: 1px solid var(--gd-border) !important;
      border-bottom-color: transparent !important;
    }
    /* Keyboard-only focus ring — subtle teal glow (doesn't add layout). */
    .nav-tabs > li > a:focus-visible {
      box-shadow: 0 0 0 3px rgba(15,118,110,.15) !important;
    }

    /* ===== DataTables (DT) ===== */
    table.dataTable {
      font-family: var(--gd-font-body);
      font-size: 13.5px;
      border-collapse: separate;
    }
    table.dataTable thead th, table.dataTable thead td {
      border-bottom: 1px solid var(--gd-border);
      color: var(--gd-text-mute);
      font-size: 11.5px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.03em;
    }
    table.dataTable tbody td { padding: 10px 12px; }
    table.dataTable.hover tbody tr:hover,
    table.dataTable tbody tr:hover { background: rgba(0,0,0,.035); }
    .dataTables_wrapper .dataTables_paginate .paginate_button {
      border-radius: 6px;
      color: var(--gd-text) !important;
    }
    .dataTables_wrapper .dataTables_paginate .paginate_button.current,
    .dataTables_wrapper .dataTables_paginate .paginate_button.current:hover {
      background: var(--gd-accent) !important;
      border-color: var(--gd-accent) !important;
      color: #fff !important;
    }

    /* ===== Notifications ===== */
    .shiny-notification {
      border-radius: 8px;
      border: 1px solid var(--gd-border);
      background: var(--gd-panel);
      color: var(--gd-text);
      box-shadow: var(--gd-shadow-sm);
    }

    /* ===== Legend text (existing helper) ===== */
    .gd-legend {
      font-size: 0.85em;
      color: var(--gd-text-mute);
      margin-top: 10px;
      max-width: 800px;
    }
    .gd-legend ul { margin: 4px 0 0 0; padding-left: 20px; }
    .gd-legend li { margin-bottom: 3px; }

    /* ===== fileInput ===== */
    /* Extra padding so long filenames don't creep under the Browse button. */
    .shiny-input-container .input-group .form-control { padding-left: 16px; }
    /* Upload progress bar — teal, tall enough to read the % text. */
    .shiny-file-input-progress {
      height: 24px !important;
      margin-top: 6px;
      background: var(--gd-bg-2);
      border-radius: 6px;
    }
    .shiny-file-input-progress .progress-bar {
      height: 24px !important;
      line-height: 24px;
      font-size: 13px;
      background: var(--gd-accent);
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
