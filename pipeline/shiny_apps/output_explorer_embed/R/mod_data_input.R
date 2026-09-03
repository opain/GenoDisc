dataInputUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title="Data Input",
    br(),
    p("This is an application for visualising the output of GenoDisc. To start, upload the 'bundle.tar.gz' file output by the GenoDisc pipeline (a legacy 'results_package.rds' also works), select one or more GWAS, and use the tabs to view interactive tables and plots of your results."),
    hr(),
    h5("Choose a bundle (.tar.gz) or legacy .rds file"),
    fileInput(ns("file"), NULL),
    h6('Or'),
    actionButton(ns("loadExample"), "Use example data"),
    tags$div(style = "font-size: 0.85em; color: #6c757d; margin-top: 6px;",
      "Example data: GenoDisc results for amyotrophic lateral sclerosis (ALS), ",
      "generated from the European-ancestry GWAS meta-analysis of ",
      tags$a("Van Rheenen et al. (2021, Nature Genetics)",
             href = "https://pubmed.ncbi.nlm.nih.gov/34873335/",
             target = "_blank", rel = "noopener noreferrer"), "."
    ),

    shinyjs::hidden(
      div(id = ns("gwas_selector_section"),
        hr(),
        h5("Select one or more GWAS"),
        tags$p(
          style = "font-size: 0.85em; color: var(--gd-text-mute); margin-bottom: 6px;",
          "Selecting more than one GWAS switches the enrichment tabs to cross-GWAS comparison views."
        ),
        selectizeInput(
          ns("gwas_selector"), NULL,
          choices = NULL, multiple = TRUE,
          options = list(plugins = list('remove_button'))
        ),
        br(), br(), br(), br(), br()
      )
    )
  )
}

dataInputServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    rds_path <- reactiveVal('')

    # Uploads land at a random datapath with no extension. gd_open dispatches
    # on file extension (tarball vs .rds), so we rename to preserve it.
    observeEvent(input$file, {
      orig <- input$file$name
      ext  <- if      (grepl("\\.tar\\.gz$", orig, ignore.case = TRUE)) ".tar.gz"
              else if (grepl("\\.tgz$",      orig, ignore.case = TRUE)) ".tar.gz"
              else if (grepl("\\.rds$",      orig, ignore.case = TRUE)) ".rds"
              else ""
      path <- input$file$datapath
      if (nzchar(ext) && !grepl(paste0(gsub("\\.", "\\\\.", ext), "$"), path)) {
        new_path <- paste0(path, ext)
        file.rename(path, new_path)
        path <- new_path
      }
      rds_path(path)
    })

    observeEvent(input$loadExample, {
      # Path is resolved against the module source file's directory (not
      # getwd()) so the button works whether the app was launched from the
      # app dir or not.
      here <- tryCatch(dirname(sys.frame(1L)$ofile), error = function(e) getwd())
      candidates <- c(
        file.path(here, "..", "data", "als_bundle.tar.gz"),
        file.path("data", "als_bundle.tar.gz")
      )
      hit <- candidates[file.exists(candidates)][1]
      if (!is.na(hit)) {
        rds_path(normalizePath(hit, winslash = "/"))
      } else {
        showNotification("Example data file not found (looked for data/als_bundle.tar.gz).", type = "error")
      }
    })

    gwas_data <- reactive({
      req(rds_path() != '')
      gd <- tryCatch(gd_open(rds_path()), error = function(e) {
        showNotification(paste0("Could not open file: ", conditionMessage(e)), type = "error")
        NULL
      })
      if (is.null(gd)) req(FALSE)
      if (length(gd_gwas(gd)) == 0L) {
        showNotification("Results package contains no GWAS.", type = "error")
        req(FALSE)
      }
      gd
    })

    # When a new bundle loads, default the selection to ALL GWAS so multi-GWAS
    # bundles land straight in comparison mode (single-GWAS bundles stay in the
    # existing single-view flow because length(input$gwas_selector) == 1).
    observeEvent(gwas_data(), {
      names_ <- gd_gwas(gwas_data())
      updateSelectizeInput(session, "gwas_selector", choices = names_, selected = names_)
      shinyjs::show("gwas_selector_section")
    })

    # Vector of currently-selected GWAS in the order the user chose (selectize
    # returns the picked order). Empty is coerced to the first GWAS available
    # so downstream reactives always have something to work with.
    selected_gwas_multi <- reactive({
      req(gwas_data())
      names_ <- gd_gwas(gwas_data())
      sel <- input$gwas_selector
      if (is.null(sel) || length(sel) == 0) return(names_[1])
      sel <- sel[sel %in% names_]
      if (length(sel) == 0) names_[1] else sel
    })

    # Scalar "current" GWAS for single-GWAS views. Kept for backward
    # compatibility with the existing per-GWAS module servers; equals the
    # first-selected GWAS.
    selected_gwas <- reactive({
      selected_gwas_multi()[1L]
    })

    comparison_mode <- reactive({
      length(selected_gwas_multi()) > 1L
    })

    list(
      gwas_data           = gwas_data,
      selected_gwas       = selected_gwas,
      selected_gwas_multi = selected_gwas_multi,
      comparison_mode     = comparison_mode
    )
  })
}
