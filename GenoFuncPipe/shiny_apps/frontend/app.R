library(shiny)

ui <- fluidPage(
  
  titlePanel("GenoFunc"),
  
  tags$hr(),
  
  h3("Enter email address:"),
  
  textInput(inputId="email", label = NULL, value = ""),
  
  tags$hr(),
  
  h3("Choose GWAS summary statistics file:"),
  h5("Note. file must be gzipped (Max. 600Mb)"),
  
  fileInput("sumstats", 
            label=NULL,
            multiple = F,
            accept = ".gz"),
  
  tags$hr(),
  
  h3("Specify column names:"),
  
  h5("Note. Either rsID or chromosome and position columns are required."),
  h5("Note. Either BETA or OR columns are required."),
  h5("Note. Either P or SE columns are required."),
  
  tags$br(),
  
  textInput('chr_col', "Chromosome (optional):"),
  textInput('bp_col', "Position (optional):"),
  textInput('rsid_col', "rsID (optional):"),
  
  textInput('a1_col', "Allele 1 (A1, Effect allele, required):"),
  textInput('a2_col', "Allele 2 (A2, required):"),
  
  textInput('beta_col', 'BETA (optional):'),
  textInput('or_col', "Odds ratio (OR, optional):"),
  textInput('p_col', "p-value (P, optional):"),
  
  textInput('se_col', "Standard error of BETA or log(OR) (optional):"),
  textInput('n_col', 'Sample Size (optional):'),
  textInput('n_col', 'Frequency of A1 (optional):'),
  textInput('n_col', "Imputation quality (INFO) (optional):"),
  
  tags$hr(),
  
  h3("Select analyses"),

  checkboxInput("ldsc_logical", "LD-score regression intercept and heritability", TRUE),
  
  checkboxInput("clump_logical", "LD-based clumping", TRUE),
  checkboxInput("cojo_logical", "COJO analysis", TRUE),
  checkboxInput("finemap_logical", "SNP-based finemapping", TRUE),
  checkboxInput("magma_gene_logical", "MAGMA gene association", TRUE),
  
  checkboxInput("twas_panel_psychencode_logical", "TWAS using PsychENCODE DLPFC data", TRUE),
  checkboxInput("smr_expression_panel_psychencode_logical", "SMR using PsychENCODE DLPFC data", TRUE),
  
  checkboxInput("smr_expression_panel_metabrain_basalganglia_logical", "SMR using MetaBrain basalganglia data", TRUE),
  checkboxInput("smr_expression_panel_metabrain_cerebellum_logical", "SMR using MetaBrain cerebellum data", TRUE),
  checkboxInput("smr_expression_panel_metabrain_cortex_logical", "SMR using MetaBrain cortex data", TRUE),
  checkboxInput("smr_expression_panel_metabrain_hippocampus_logical", "SMR using MetaBrain hippocampus data", TRUE),
  checkboxInput("smr_expression_panel_metabrain_spinalcord_logical", "SMR using MetaBrain spinalcord data", TRUE),
    
  checkboxInput("pwas_panel_rosmap_logical", "PWAS using ROSMAP DLPFC data", TRUE),
  checkboxInput("pwas_panel_banner_logical", "PWAS using Banner et al. DLPFC data", TRUE),
  checkboxInput("smr_protein_panel_rosmap_logical", "SMR using ROSMAP DLPFC data", TRUE),

  checkboxInput("magma_drugtargetor_logical", "Drug enrichment analysis using MAGMA and Drug Targetor", TRUE),
  checkboxInput("twas_gsea_lincs_logical", "Drug enrichment analysis using TWAS and CMAP and Drug Targetor (TWAS-GSEA method)", TRUE),
  checkboxInput("twas_so_lincs_logical", "Drug enrichment analysis using TWAS and CMAP and Drug Targetor (Average ranks method)", TRUE),
  
  tags$hr(),
  
  actionButton("submit_button","Submit Job"),
  
  tags$hr(),
  
  htmlOutput("text"),
  
  tags$hr(),
  
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  options(shiny.maxRequestSize=600*1024^2) 
  
  output$text <- renderText({
    
    #preparing progress tracker      
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "", value = 0)
    updateProgress <- function(value = NULL, detail = NULL, max=NULL) {
      if(is.null(max))max <- 50
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + 1/max
      }
      progress$set(value = value, detail = detail)
    }

    # Take a dependency on input$submit_button
    if(input$submit_button > 0){
      path <- isolate(input$sumstats[["datapath"]])
      filename <- isolate(input$sumstats[["name"]])
      email <- isolate(input$email)
      
      if(is.null(path)){
        return("No GWAS summary statistics file has been selected - make sure the upload tracker says 'upload completed' before clicking 'Submit Job'")
      }
      if(is.null(email)){
        return("No email address has been specified")
      }
      
      # There will need to be a function initate the analysis
      out <- 'Your job has been submitted. You will receive an email when the analysis is complete.'
      
      return(out)
      
    }
  })
}

# Run the app ----
shinyApp(ui, server)

# To do:
# - Disable submit job until required information has been specified. This can be done using shinyjs

