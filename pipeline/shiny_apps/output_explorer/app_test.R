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

# Loads functions
source('functions.R')

options(shiny.maxRequestSize = 600 * 1024 * 1024)

# Define UI for the Shiny app
ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  tags$style(HTML("
    .navbar {
      min-height: 80px; /* Set the minimum height of the navbar */
    }
    .navbar .nav.navbar-nav li {
      margin: 15px 0; /* Adjust vertical spacing for buttons */
    }
    .footer {
      display: flex;
      justify-content: space-between;
      align-items: center;
      text-align: left;
      padding: 20px;
      background-color: #f5f5f5; /* Light grey background color */
    }
    .disabled-button {
        pointer-events: none;
        opacity: 0.4;
    }
    .progress {
        height: 22px;
        width: 500px;
    }
    .input-group .form-control {
        padding: 10px;
        font-size: 12px;
        width: 400px;
    }
    .load1 .loader {
        margin: 88px;
        position: relative;
    }
    .custom-panel {
          padding: 20px; /* Adjust the padding value as needed */
          border: 1px solid #ccc; /* Optional: Add a border for better visibility */
    }
  ")),
  
  navbarPage(
    title = div(img(src='logo/horizontal.png',
                    style="margin-top: -14px;
                               padding-right:10px;
                               padding-bottom:0px",
                    height = 80)),
    theme = shinythemes::shinytheme("paper"),
    
    tabPanel(
      title='Explore',
      
      tabsetPanel(
        tabPanel(
          title="Data Input",
          br(),
          p("This is an application for visualising the output of GenoDisc. To start, upload the 'results_package.rds' file output by the GenoDisc pipeline, select a GWAS, and use the tabs to view interactive tables and plots of your results."), 
          p("Click ",a("here", href = "https://github.com/opain/GenoDisc"), " here to learn more about the pipeline. Please cite  ",a("our publication", href = "https://github.com/opain/GenoDisc"), "  and relevent software and datasets included in your analysis."), 
          hr(),
          h5("Choose an .RDS file"),
          fileInput("file", NULL),
          h6('Or'),
          actionButton("loadExample", "Use example data"),
          
          hr(),
          
          h5("Select a GWAS"),
          selectInput("gwas_selector", NULL, ""),
          br(),
          br(),
          br(),
          br(),
          br()
        ),
        tabPanel(
          title="Molecular Associations",
          br(),
          p("This tab shows molecular association results. Select the tabs below to see a summary of results across methods, or method-specific results tables."), 
          hr(),
          tabsetPanel(
            tabPanel(
              title="Summary",
              br(),
              
              fluidPage(
                sidebarPanel(
                  selectInput("selected_methods_mol", "Select methods", "", multiple=T),
                  selectInput("selected_expr_panels_mol", "Select expression panels", "", multiple=T),
                  selectInput("selected_protein_panels_mol", "Select protein panels", "", multiple=T),
                  radioButtons("conf_only_mol", "Show high-confidence only :",
                               choices = c("True" = T,
                                           "False" = F),
                               selected = T),
                  textInput("geneInput_mol", "Enter gene symbols (whitespace- or comma-seperated):")
                ),
                
                mainPanel(
                  uiOutput("message_too_large_mol"),
                  uiOutput("message_no_genes_mol"),
                  uiOutput("mol_assoc_plot.ui")
                )
              )
            )
          )
        )
      )
    )
  ),
  
  br(),
  br(),
  br(),
  br(),
  
  fluidRow(
    class = "footer",
    
    column(6, style = "text-align: left;",
           div("Created by Oliver Pain with colleagues at King's College London."),
    ),
    
    column(6, style = "text-align: right;",
           img(src = "kcl.svg", alt = "Kings College London", style = "height: 60px;"),
           img(src = "Wellcome_logo.png", alt = "The Wellcome Trust", style = "height: 60px;"),
           img(src = "Turing_logo.png", alt = "The Turing Institute", style = "height: 60px;")
    )
  )
)

# Define server logic
server <- function(input, output, session) {

  rds_path <- reactiveVal('')
  
  observeEvent(input$file, {
    rds_path(input$file$datapath)
  })
  
  observeEvent(input$loadExample, {
    rds_path('example.rds')
  })
  
  gwas_data <- reactive({
    req(!is.null(input$file) | input$loadExample > 0)
    readRDS(rds_path())
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
  
  observeEvent(selected_gwas(), {
    req(gwas_data(), selected_gwas())
    
    mol_assoc_summary_data <- reactive({
      print(selected_gwas())
      print(names(gwas_data()))
      print(gwas_data()[[selected_gwas()]]$mol_assoc$finemap$L1)
      
      data.frame(Panel = "SuSie (L=1)",
                 ID=gwas_data()[[selected_gwas()]]$mol_assoc$finemap$L1,
                 Z=1,
                 Sig=F,
                 Coloc=F,
                 Method="SNP\nFine-mapping",
                 Type='')
    })
    
    observeEvent(mol_assoc_summary_data(), {
      all_func_res<-mol_assoc_summary_data()
    })
  })
}

# Run the Shiny app
shinyApp(ui, server)






