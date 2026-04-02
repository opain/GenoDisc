dataInputUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title="Data Input",
    br(),
    p("This is an application for visualising the output of GenoDisc. To start, upload the 'results_package.rds' file output by the GenoDisc pipeline, select a GWAS, and use the tabs to view interactive tables and plots of your results."),
    p("Click ",a("here", href = "https://github.com/opain/GenoDisc"), " here to learn more about the pipeline. Please cite  ",a("our publication", href = "https://github.com/opain/GenoDisc"), "  and relevent software and datasets included in your analysis."),
    hr(),
    h5("Choose an .RDS file"),
    fileInput(ns("file"), NULL),
    h6('Or'),
    actionButton(ns("loadExample"), "Use example data"),

    hr(),

    h5("Select a GWAS"),
    selectInput(ns("gwas_selector"), NULL, ""),
    br(),
    br(),
    br(),
    br(),
    br()
  )
}

dataInputServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    rds_path <- reactiveVal('')

    observeEvent(input$file, {
      rds_path(input$file$datapath)
    })

    observeEvent(input$loadExample, {
      rds_path('example.rds')
    })

    gwas_data <- reactive({
      req(rds_path() != '')
      rds <- tryCatch(readRDS(rds_path()), error = function(e) NULL)
      if (is.null(rds)) {
        showNotification("Could not read the uploaded file. Please ensure it is a valid .rds file.", type = "error")
        req(FALSE)
      }
      err <- validate_rds(rds)
      if (!is.null(err)) {
        showNotification(err, type = "error")
        req(FALSE)
      }
      rds
    })

    observeEvent(gwas_data(), {
      rds <- gwas_data()
      gwas_names <- names(rds)[names(rds) != 'configuration']
      updateSelectInput(session, "gwas_selector", choices = gwas_names)
    })

    selected_gwas<-reactive({
      req(gwas_data(), input$gwas_selector)
      gwas_selected <- ifelse(input$gwas_selector %in% names(gwas_data()),
                              input$gwas_selector,
                              names(gwas_data())[names(gwas_data()) != 'configuration'][1])
      return(gwas_selected)
    })

    list(gwas_data = gwas_data, selected_gwas = selected_gwas)
  })
}
