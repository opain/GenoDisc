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
library(colourpicker)

# Load functions and modules
source('functions.R')
source('reader.R')
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)

options(shiny.maxRequestSize = 600 * 1024 * 1024)
options(shiny.legacy.datatable = TRUE)

# Define UI
ui <- fluidPage(

  shinyjs::useShinyjs(),

  # GenoDisc-web design tokens (light + dark) mirrored from the Django app's
  # _design_system.html so the two look congruent. Dark is the default; users
  # opt into light and the choice persists via localStorage['genodisc_theme'],
  # matching Django's key and semantics.
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = ""),
    tags$link(rel = "stylesheet",
              href = paste0("https://fonts.googleapis.com/css2?",
                            "family=Inter+Tight:wght@400;500;600;700&",
                            "family=Montserrat:wght@500;600&",
                            "display=swap")),
    # Pre-paint: restore saved theme before body renders (no flash).
    tags$script(HTML("
      (function () {
        try {
          if (localStorage.getItem('genodisc_theme') === 'light') {
            document.documentElement.setAttribute('data-theme', 'light');
          }
        } catch (e) {}
      })();
    ")),
    # Wire up the toggle button once the DOM is ready.
    tags$script(HTML("
      document.addEventListener('DOMContentLoaded', function () {
        var KEY = 'genodisc_theme';
        var btn = document.getElementById('themeToggle');
        if (!btn) return;
        btn.addEventListener('click', function () {
          var light = document.documentElement.getAttribute('data-theme') === 'light';
          if (light) {
            document.documentElement.removeAttribute('data-theme');
          } else {
            document.documentElement.setAttribute('data-theme', 'light');
          }
          try { localStorage.setItem(KEY, light ? 'dark' : 'light'); } catch (e) {}
        });
      });
    "))
  ),

  tags$style(HTML("
    :root {
      /* Dark defaults — mirrors Django _design_system.html :root. */
      --gd-bg:        #0c0f14;
      --gd-bg-2:      #11151c;
      --gd-panel:     #161b24;
      --gd-panel-2:   #1c2330;
      --gd-border:    #242b3a;
      --gd-border-2:  #2e3647;
      --gd-text:      #e8ecf2;
      --gd-text-dim:  #c9d1e1;
      --gd-text-mute: #c9d1e1;
      --gd-accent:    #7df0c4;
      --gd-accent-h:  #66d9ab;
      --gd-accent-2:  #f3c969;
      --gd-accent-3:  #9eb6ff;
      --gd-danger:    #ff8985;
      --gd-shadow-sm: 0 1px 2px rgba(0,0,0,.4);
      --gd-shadow:    0 8px 24px rgba(0,0,0,.35);
      --gd-focus-ring: rgba(125,240,196,.20);
      --gd-pill-bg:    rgba(125,240,196,.10);
      --gd-pill-text:  #7df0c4;
      /* Text colour on an --gd-accent background. Accent is a bright mint in
         dark mode, so use a dark ink; light mode uses white on dark teal. */
      --gd-accent-text: #0c0f14;
      /* Subtle darker/lighter tint for hover states over --gd-panel-2 */
      --gd-hover-bg: #242b3a;
      /* Alt-row background for banded tables (noticeably lighter than panel). */
      --gd-band:    #21293a;
      --gd-r-btn: 10px;
      --gd-r-input: 6px;
      --gd-r-card: 16px;
      --gd-font-body: 'Inter Tight', -apple-system, BlinkMacSystemFont,
                     'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
      --gd-font-head: 'Montserrat', -apple-system, BlinkMacSystemFont,
                     'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
    }
    :root[data-theme='light'] {
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
      --gd-focus-ring: rgba(15,118,110,.15);
      --gd-pill-bg:    rgba(15,118,110,.10);
      --gd-pill-text:  #0c655e;
      --gd-accent-text: #ffffff;
      --gd-hover-bg:   #e6e9ef;
      --gd-band:       #f5f7fa;
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
      color: var(--gd-accent-text);
    }
    .btn-primary:hover,
    .btn-file:hover,
    .shiny-download-link:hover {
      background: var(--gd-accent-h);
      border-color: var(--gd-accent-h);
      color: var(--gd-accent-text);
    }
    .btn:focus, .btn:active,
    .action-button:focus,
    .shiny-download-link:focus {
      outline: none;
      box-shadow: 0 0 0 3px var(--gd-focus-ring);
    }

    /* ===== Form inputs ===== */
    .form-control,
    input[type='text'], input[type='number'], input[type='search'], textarea,
    .selectize-input,
    .selectize-control.single .selectize-input,
    .selectize-control.multi .selectize-input {
      background: var(--gd-bg-2) !important;
      border: 1px solid var(--gd-border-2) !important;
      border-radius: var(--gd-r-input);
      color: var(--gd-text) !important;
      font-family: var(--gd-font-body);
      font-size: 13.5px;
      padding: 8px 12px;
      box-shadow: none !important;
    }
    /* fileInput: hide the 'no file selected' text field so only the Browse
       button remains, styled as a complete standalone pill. Filename
       feedback still comes via the progress bar and the tab reveal after
       upload completes. */
    .input-group > .form-control { display: none !important; }
    .input-group > .input-group-btn > .btn,
    .input-group > .input-group-btn > label > .btn,
    .input-group > .input-group-btn > label.btn {
      padding: 8px 14px !important;
      font-size: 13.5px !important;
      line-height: 1.5 !important;
      border: 1px solid var(--gd-accent) !important;
      border-radius: var(--gd-r-input) !important;
    }
    .form-control:focus,
    input[type='text']:focus, input[type='number']:focus, textarea:focus,
    .selectize-input.focus {
      border-color: var(--gd-accent);
      box-shadow: 0 0 0 3px var(--gd-focus-ring);
      outline: none;
    }
    .selectize-dropdown, .selectize-dropdown-content {
      border: 1px solid var(--gd-border);
      border-radius: var(--gd-r-input);
      background: var(--gd-panel);
      color: var(--gd-text);
    }
    /* Dropdown option base colour — needs !important to defeat selectize's
       own stylesheet, which otherwise wins with equal specificity. */
    .selectize-dropdown .option,
    .selectize-dropdown [data-selectable] {
      color: var(--gd-text) !important;
      background: var(--gd-panel) !important;
    }
    /* Highlighted / hovered option — accent tint with readable text. The
       previous panel-2 background was almost indistinguishable from the
       base panel colour in both themes. */
    .selectize-dropdown .option:hover,
    .selectize-dropdown .option.active,
    .selectize-dropdown [data-selectable]:hover,
    .selectize-dropdown [data-selectable].active {
      background: var(--gd-accent) !important;
      color: var(--gd-accent-text) !important;
    }
    .selectize-input > input { color: var(--gd-text); }
    .selectize-control.multi .selectize-input > div {
      background: var(--gd-pill-bg);
      color: var(--gd-pill-text);
      border-radius: 4px;
      padding: 2px 8px;
    }
    /* Slider (ionRangeSlider) — teal accent */
    .irs--shiny .irs-bar { background: var(--gd-accent); border-color: var(--gd-accent); }
    .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
      background: var(--gd-accent);
      /* Text on --gd-accent needs the same colour as button text: dark in
         dark mode (where accent is a bright mint), white in light mode. */
      color: var(--gd-accent-text);
    }
    .irs--shiny .irs-from::before, .irs--shiny .irs-to::before,
    .irs--shiny .irs-single::before {
      border-top-color: var(--gd-accent);
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
    .gd-details > summary:hover {
      background: var(--gd-hover-bg);
      color: var(--gd-text);
    }
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
      box-shadow: 0 0 0 3px var(--gd-focus-ring) !important;
    }

    /* ===== DataTables (DT) ===== */
    table.dataTable, table.dataTable tbody, table.dataTable tbody tr,
    table.dataTable tbody td, table.dataTable thead th, table.dataTable thead td,
    .dataTables_wrapper {
      background: var(--gd-panel);
      color: var(--gd-text);
      font-family: var(--gd-font-body);
    }
    table.dataTable { font-size: 13.5px; border-collapse: separate; }
    table.dataTable thead th, table.dataTable thead td {
      border-bottom: 1px solid var(--gd-border);
      color: var(--gd-text-mute);
      font-size: 11.5px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.03em;
    }
    table.dataTable tbody td { padding: 10px 12px; }
    /* Zebra rows — apply to <td> cells directly, at high specificity, so we
       beat DataTables' built-in .stripe rules and any inline styles it may
       inject. :nth-child covers legacy shiny tables that don't get the
       .odd/.even classes; the .odd/.even rules cover DT's default markup. */
    .dataTables_wrapper table tbody tr:nth-child(even) > td,
    .dataTables_wrapper table tbody tr.even > td,
    table.dataTable tbody tr:nth-child(even) > td,
    table.dataTable tbody tr.even > td {
      background-color: var(--gd-panel) !important;
    }
    .dataTables_wrapper table tbody tr:nth-child(odd) > td,
    .dataTables_wrapper table tbody tr.odd > td,
    table.dataTable tbody tr:nth-child(odd) > td,
    table.dataTable tbody tr.odd > td {
      background-color: var(--gd-band) !important;
    }
    .dataTables_wrapper table tbody tr:hover > td,
    table.dataTable tbody tr:hover > td,
    table.dataTable.hover tbody tr:hover > td {
      background-color: var(--gd-panel-2) !important;
    }
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

    /* ===== First-load: hide downstream tabs until data is loaded =====
       Fires at CSS-parse time so there's no flash of all tabs while Shiny's
       observeEvent(session$clientData, hideTab(...)) is still waiting for
       the WebSocket. Once config_flags first resolves in the server we set
       data-app-ready on <html>, and this rule stops matching — the tabs
       Shiny showed via showTab() become visible; hideTab-hidden ones stay
       hidden thanks to the inline display:none Shiny set. */
    html:not([data-app-ready]) #main_tabs > li:not(:first-child) {
      display: none !important;
    }

    /* ===== Plot outputs =====
       Round the corners of every plotOutput visually via overflow:hidden on
       the wrapper. The underlying <img> stays rectangular, so right-click →
       Save Image As downloads the original unrounded PNG. */
    .shiny-plot-output {
      border-radius: var(--gd-r-card);
      overflow: hidden;
    }
    /* Same rounding for the pre-rendered PNG plots that some tabs load from
       the bundle (GWAS QC, SNP assoc). These are embedded as base64 data
       URLs — border-radius on <img> is a render-time effect only, so
       right-click Save Image As still gets the rectangular original. */
    img[src^='data:image/png'] {
      border-radius: var(--gd-r-card);
    }

    /* ===== Theme toggle button (fixed top-right) ===== */
    .gd-theme-toggle {
      position: fixed;
      top: 10px;
      right: 14px;
      z-index: 1050;
    }
    .gd-theme-toggle button {
      background: transparent;
      border: 1px solid var(--gd-border);
      color: var(--gd-text) !important;
      padding: 6px 10px !important;
      border-radius: 8px !important;
      line-height: 1;
      cursor: pointer;
    }
    .gd-theme-toggle button:hover {
      background: var(--gd-panel-2);
      border-color: var(--gd-border) !important;
    }
    .gd-theme-toggle svg { display: inline-block; vertical-align: middle; }
    .gd-theme-toggle .icon-moon { display: none; }
    :root[data-theme='light'] .gd-theme-toggle .icon-sun  { display: none; }
    :root[data-theme='light'] .gd-theme-toggle .icon-moon { display: inline-block; }
  ")),

  theme = shinythemes::shinytheme("paper"),

  tags$div(
    class = "gd-theme-toggle",
    tags$button(
      id = "themeToggle",
      type = "button",
      `aria-label` = "Toggle light or dark theme",
      title = "Toggle light/dark theme",
      HTML('<svg class="icon-sun"  width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.93 4.93l1.41 1.41M17.66 17.66l1.41 1.41M2 12h2M20 12h2M4.93 19.07l1.41-1.41M17.66 6.34l1.41-1.41"/></svg>'),
      HTML('<svg class="icon-moon" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>')
    )
  ),

  tabsetPanel(id = "main_tabs",
    dataInputUI("data_input"),
    overviewUI("overview"),
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

  # Hide all downstream tabs on startup (including the new Overview tab, which
  # only appears once the user selects >= 2 GWAS).
  observeEvent(session$clientData, {
    for (tab in c("Overview", "GWAS QC", "SNP Associations",
                   "Molecular Associations", "Enrichment Analysis",
                   "References", "Configuration")) {
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

    # Release the CSS-level hide of downstream tabs now that Shiny has
    # decided which of them should be visible for the loaded bundle.
    shinyjs::runjs("document.documentElement.setAttribute('data-app-ready', '1');")
  })

  # ---- Comparison mode plumbing ----

  # Show/hide the Overview tab based on how many GWAS the user has selected.
  # The user stays on whichever tab they were on; they can pick Overview
  # themselves when they want to see it.
  observe({
    if (shared$comparison_mode()) showTab("main_tabs", "Overview")
    else                          hideTab("main_tabs", "Overview")
  })

  # Cross-GWAS long tibble; rebuilds only when the bundle or the GWAS set
  # changes, not when per-view filter inputs change.
  comparison_long <- reactive({
    req(shared$gwas_data())
    req(length(shared$selected_gwas_multi()) >= 1)
    build_comparison_long(shared$gwas_data(),
                           shared$selected_gwas_multi(),
                           entity_types = c("tissue", "atc", "gene"))
  })

  # Wire up all modules
  overviewServer("overview", shared$gwas_data, shared$selected_gwas_multi,
                  shared$comparison_mode, comparison_long)
  gwasQcServer("gwas_qc", shared$gwas_data, shared$selected_gwas, gwas_list, config_flags)
  snpAssocServer("snp_assoc", shared$gwas_data, shared$selected_gwas, config_flags)
  molAssocServer("mol_assoc", shared$gwas_data, shared$selected_gwas, config_flags,
                  selected_gwas_multi = shared$selected_gwas_multi,
                  comparison_mode     = shared$comparison_mode,
                  comparison_long     = comparison_long)
  enrichmentServer("enrichment", shared$gwas_data, shared$selected_gwas, config_flags,
                    selected_gwas_multi = shared$selected_gwas_multi,
                    comparison_mode     = shared$comparison_mode,
                    comparison_long     = comparison_long)
  referencesServer("references")
  configurationServer("configuration", shared$gwas_data)
}

# Run the Shiny app
shinyApp(ui, server)
