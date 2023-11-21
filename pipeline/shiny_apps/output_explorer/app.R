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
  ")),
  
  navbarPage(
    title = div(img(src='logo/horizontal.png',
                    style="margin-top: -14px;
                               padding-right:10px;
                               padding-bottom:0px",
                    height = 80)),
    theme = shinythemes::shinytheme("paper"),
    
    tabPanel(
      title='Home',
      
      div(
        style = "display: flex; align-items; padding: 20px;",
        div(
          style = "flex: 4; padding-right: 20px;",
          h3("Welcome to GenoDisc"),
          p("GenoDisc is a comprehensive platform for Genome-Wide Association Study (GWAS) summary statistics analysis. Explore genetic associations, visualize results, and gain insights into your data with our user-friendly tools."),
          p("Get started by uploading your GWAS summary statistics and let GenoDisc help you uncover meaningful patterns in your data."),
          hr(),
          
          # Instructions for submitting data
          h4("Submitting GWAS Summary Statistics"),
          p("You can submit your GWAS summary statistics for analysis via our King's College London server by going to the 'Submit' tab."),
          hr(),
          
          # Instructions for exploring existing results
          h4("Exploring Previous Results"),
          p("If you have already run the GenoDisc pipeline and have results, you can explore them by uploading your results_package.rds file to the 'Explore' tab."),
          hr(),
          
          # Information on citing the platform
          h4("Citing GenoDisc"),
          p("If you use GenoDisc for a publication, please cite our publication describing the platform, as well as the underlying datasets and methods it uses."),
          p("Publication reference: TBD"),
          hr() ,
          
          # Information on citing the platform
          h4("Our Newsletter"),
          p("If you would like to receive updates regarding GenoDisc, provide you email address below:"),
          textInput(inputId="email_newsletter", label = NULL, value = NULL),
          actionButton("email_button","Sign Up"),
          uiOutput("email_submit_text"),
        ),
        
        div(
          style = "flex: 2; padding-left: 20px; text-align: center;",
          img(
            src = "schematic.png",
            alt = "Schematic",
            style = "max-width: 70%; height: auto;"
          )
        )
      )
    ),
    tabPanel(
      title='Submit',
      
      tabsetPanel(
        tabPanel(
          title="Required",
          br(),
          h5("Email Address:"),
          textInput(inputId="email", label = NULL, value = NULL),
          br(),
          h5("GWAS sumstats:"),
          p("Must be gzipped. Max. 600Mb."),
          fileInput("sumstats", 
                    label=NULL,
                    multiple = F,
                    accept = ".gz"),
          br(),
          
          # This will be a table showing the column names in the sumstats, and how they have been interpreted.
          withSpinner(uiOutput("submit_colnames.ui")),
          uiOutput("header_check_messages.ui"),
          uiOutput("specify_n.ui"),
          uiOutput("specify_population.ui"),
          
        ),
        
        tabPanel(
          title="Options",
          br(),
          h5('Heritability Estimation:'),
          checkboxInput("ldsc_logical", "LD-Score Regression", F),
          uiOutput("specify_pop_prev.ui"),
          br(),
          h5('SNP Association'),
          checkboxInput("clump_logical", "LD-based clumping", F),
          checkboxInput("cojo_logical", "COJO analysis", F),
          checkboxInput("finemap_logical", "SNP-based finemapping", F),
          br(),
          h5('Molecular Association'),
          fluidRow(
            column(width=3,
                   h6('Positional Mapping'),
                   checkboxInput("magma_gene_logical", "MAGMA gene association", F)
            ),
            column(width=3,
                   h6('Expression'),
                   selectInput("twas_panels", "Select TWAS panels:", choices=twas_panel_names, multiple=T),
                   selectInput("smr_expr_panels", "Select SMR expression panels:", choices=smr_expr_panel_names, multiple=T)
            ),
            column(width=3,
                   h6('Protein'),
                   selectInput("pwas_panels", "Select PWAS panels:", choices=pwas_panel_names, multiple=T),
                   selectInput("smr_protein_panels", "Select SMR protein panels:", choices=smr_protein_panel_names, multiple=T)
            )
          ),
          br(),
          h5('Enrichment'),
          selectInput("enrich_method", "Select methods:", c('MAGMA','TWAS-GSEA'), multiple=T),
          fluidRow(
            column(width=3,
                   h6('Drug'),
                   selectInput("enrich_data", "Select datasets:", c('DrugTargetor - Drug-Gene Interactions'), multiple=T),
            ),
          ),
        )
      ),
      
      hr(),
      
      actionButton("submit_button", "Submit Job", class = "btn-primary"),
      uiOutput("submit_text"),
      
      # There should be a way for users to cancel their submitted job.
      # Perhaps by emailing me with a specific subject, from the associated email address.

      br(),
      br(),
    ),
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
          title="GWAS QC",
          br(),
          p("This tab shows key quality control statistics for your selected GWAS."), 
          fluidRow(
            column(width=4,
                   dataTableOutput("qc_table"),
            )
          )
        ),
        tabPanel(
          title="SNP Associations",
          br(),
          p("This tab shows SNP association results. Select the Lead variant tab below to view information for independent lead variants identified by either LD-based clumping or COJO. Select the Fine-mapping tab below to view SuSiE Finemaping results."), 
          hr(),
          tabsetPanel(
            tabPanel(
              title="Lead variants",
              br(),            
              fluidPage(
                sidebarPanel(
                  radioButtons("clumping_type", "Select method:",
                               choices = c("COJO" = "cojo_analysis",
                                           "LD-based clumping" = "ld_clumping"),
                               selected = "cojo_analysis"),
                  radioButtons("pvalue_threshold", "Select P-value Threshold:",
                               choices = c("Genome-wide significance (p < 5e-8)" = 5e-8, 
                                           "Suggestive significance (p < 1e-5)" = 1e-5), 
                               selected = 5e-8),
                  width = 3
                ),
                
                mainPanel(
                  dataTableOutput("snp_assoc_lead_table"),
                  width = 9
                )
              )
            ),
            tabPanel(
              title="Fine-mapping",
              br(),
              fluidPage(
                sidebarPanel(
                  radioButtons("l_param", "Select L parameter:",
                               choices = c("L1" = "L1",
                                           "L10" = "L10"),
                               selected = "L1"),
                  width = 3
                ),
                
                mainPanel(
                  dataTableOutput("snp_assoc_finemap_table"),
                  width = 6
                )
              )
            )
          )
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
            ),
            tabPanel(
              title="MAGMA",
              br(),
              p("This tab shows MAGMA gene association results."), 
              hr(),
              br(),
              fluidRow(
                column(width=6,
                       dataTableOutput("mol_assoc_magma_table"),
                )
              ),
              br()
            ),
            tabPanel(
              title="Expression - FUSION",
              br(),
              p("This tab shows differential expression association results from FUSION."), 
              hr(),
              br(),
              fluidRow(
                column(width=9,
                       dataTableOutput("mol_assoc_fusion_expr_table"),
                )
              ),
              br()
            ),
            tabPanel(
              title="Protein - FUSION",
              br(),
              p("This tab shows differential protein level association results from FUSION."), 
              hr(),
              br(),
              fluidRow(
                column(width=9,
                       dataTableOutput("mol_assoc_fusion_protein_table"),
                )
              ),
              br()
            ),
            tabPanel(
              title="Expression - SMR",
              br(),
              p("This tab shows differential expression association results from SMR"), 
              hr(),
              br(),
              fluidRow(
                column(width=9,
                       dataTableOutput("mol_assoc_smr_expr_table"),
                )
              ),
              br()
            ),
            tabPanel(
              title="Protein - SMR",
              br(),
              p("This tab shows differential protein level association results from SMR"), 
              hr(),
              br(),
              fluidRow(
                column(width=9,
                       dataTableOutput("mol_assoc_smr_protein_table"),
                )
              ),
              br()
            )
          )
        ),
        tabPanel(
          title="Enrichment Analysis",
          br(),
          p("This tab shows enrichment analysis results. Select the tabs below to see results for your desired gene annotations and methods"), 
          hr(),
          tabsetPanel(
            tabPanel(
              title="Drug",
              br(),
              tabsetPanel(
                tabPanel(
                  title="Summary",
                  br(),
                  
                  fluidPage(
                    sidebarPanel(
                      selectInput("selected_methods_drug", "Select methods", "", multiple=T),
                      selectInput("selected_expr_panels_drug", "Select expression panels", "", multiple=T),
                      radioButtons("conf_only_drug", "Show FDR significant only :",
                                   choices = c("True" = T,
                                               "False" = F),
                                   selected = T),
                      textInput("drugInput_drug", "Search drug (whitespace- or comma-seperated):"),
                      textInput("atcInput_drug", "Search ATC Code (whitespace- or comma-seperated):"),
                      selectInput("selected_sort_drug", "Sort by:", '', multiple = F)
                      
                    ),
                    
                    mainPanel(
                      uiOutput("message_too_large_drug"),
                      uiOutput("message_no_drugs_drug"),
                      uiOutput("tx_drug_plot.ui")                      
                    )
                  )
                ),
                tabPanel(
                  title="MAGMA",
                  br(),
                  p("This tab shows MAGMA drug enrichment results."), 
                  hr(),
                  br(),
                  fluidRow(
                    column(width=9,
                           dataTableOutput("tx_drug_magma_table"),
                    )
                  ),
                  br()
                ),
                tabPanel(
                  title="GCSC",
                  br(),
                  p("This tab shows GCSC drug enrichment results."), 
                  hr(),
                  br(),
                  fluidRow(
                    column(width=9,
                           dataTableOutput("tx_drug_gcsc_table"),
                    )
                  ),
                  br()
                ),
                tabPanel(
                  title="TWAS-GSEA",
                  br(),
                  p("This tab shows TWAS-GSEA drug enrichment results."), 
                  hr(),
                  br(),
                  fluidRow(
                    column(width=9,
                           dataTableOutput("tx_drug_twas_gsea_table"),
                    )
                  ),
                  br()
                )
              )
            ),
            tabPanel(
              title="ATC",
              br(),
              tabsetPanel(
                tabPanel(
                  title="Summary",
                  br(),
                  
                  fluidPage(
                    sidebarPanel(
                      selectInput("selected_methods_atc", "Select methods", "", multiple=T),
                      selectInput("selected_expr_panels_atc", "Select expression panels", "", multiple=T),
                      radioButtons("conf_only_atc", "Show FDR significant only :",
                                   choices = c("True" = T,
                                               "False" = F),
                                   selected = T),
                      textInput("atcInput_atc", "Search ATC Code (whitespace- or comma-seperated):"),
                      selectInput("selected_sort_atc", "Sort by:", '', multiple = F)
                      
                    ),
                    
                    mainPanel(
                      uiOutput("message_too_large_atc"),
                      uiOutput("message_no_atcs_atc"),
                      uiOutput("tx_atc_plot.ui")                      
                    )
                  )
                ),
                tabPanel(
                  title="MAGMA",
                  br(),
                  p("This tab shows MAGMA ATC enrichment results."), 
                  hr(),
                  br(),
                  fluidRow(
                    column(width=6,
                           dataTableOutput("tx_atc_magma_table"),
                    )
                  ),
                  br()
                ),
                tabPanel(
                  title="GCSC",
                  br(),
                  p("This tab shows GCSC ATC enrichment results."), 
                  hr(),
                  br(),
                  fluidRow(
                    column(width=6,
                           dataTableOutput("tx_atc_gcsc_table"),
                    )
                  ),
                  br()
                ),
                tabPanel(
                  title="TWAS-GSEA",
                  br(),
                  p("This tab shows TWAS-GSEA ATC enrichment results."), 
                  hr(),
                  br(),
                  fluidRow(
                    column(width=8,
                           dataTableOutput("tx_atc_twas_gsea_table"),
                    )
                  ),
                  br()
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
  
  ##############################################################################
  # Home
  ##############################################################################
  
  email_message <- reactiveVal("")
  
  observeEvent(input$email_button, {
    
    email<-input$email_newsletter
    
    if(!is_valid_email_format(email)){
      email_message("Error: Email address format is invalid\n")
    } else {
      write.table(
        email,
        "../uploads/mailing_list.txt",
        col.names = F,
        row.names = F,
        quote = T,
        append = T
      )
      email_message("You will receive a confirmation email.\n")
    }
  })
  
  output$email_submit_text <- renderUI({
    HTML(gsub("\n", "<br>", email_message()))
  })
  
  ##############################################################################
  # Submit
  ##############################################################################
  
  read_ss<-reactive({
    req(input$sumstats)
    
    # Read in the header and interpret column names
    sub_ss<-fread(input$sumstats$datapath, nrows = 1000)
    
    return(sub_ss)
  })
  
  head_interp<-reactive({
    sub_ss<-read_ss()
    
    # Read in the header and interpret column names
    sub_header<-toupper(names(sub_ss))
    
    # Remove columns that are all NA
    sub_ss_comp<-sub_ss[,apply(sub_ss, 2, function(x) !(all(is.na(x)))), with=F]
    sub_header_comp<-toupper(names(sub_ss_comp))
    
    int_header <- sub_header_comp
    for(i in names(ss_head_dict)){
      int_header[int_header %in% ss_head_dict[[i]]] <- i
    }
    int_header[!(int_header %in% unlist(ss_head_dict))]<-NA
    
    # Show original and interpreted header
    header_interp <- data.frame(Original = sub_header_comp,
                                Interpreted = int_header)
    
    # Show columns that are ignored due to be irrelevant, duplicated, or all NA
    header_interp$Keep<-!(
      !(int_header %in% names(ss_head_dict)) | 
        duplicated(int_header)
    )
    
    # Insert reason it was ignored
    header_interp$Reason<-NA
    header_interp$Reason[!(int_header %in% names(ss_head_dict))]<-'Not recognised'
    header_interp$Reason[duplicated(int_header)]<-'Duplicated'
    
    # Show columns ignored due to missingness
    for(i in sub_header[!(sub_header %in% header_interp$Original)]){
      header_interp<-rbind(header_interp, data.frame(Original = i,
                                                     Interpreted = NA,
                                                     Keep=F,
                                                     Reason = 'First 1000 rows NA'))
    }
    
    header_interp[is.na(header_interp)]<-'NA'
    header_interp$Keep<-as.character(header_interp$Keep)
    
    # Insert description of each column after interpretation
    header_labels<-data.frame(Interpreted=c('SNP','CHR','BP','A1','A2','P','OR','BETA','Z','SE','N','N_CAS','N_CON','NEF','FRQ','FRQ_A','FREQ_U','INFO'),
                              Description=c("RSID for variant",
                                            "Chromosome number",
                                            "Base pair position",
                                            "Allele 1 (effect allele)",
                                            "Allele 2",
                                            "P-value of association",
                                            "Odds ratio effect size",
                                            "BETA effect size",
                                            "Z-score",
                                            "Standard error of log(OR) or BETA",
                                            "Total sample size",
                                            "Number of cases",
                                            "Number of controls",
                                            "Effective sample size",
                                            "Allele frequency",
                                            "Allele frequency in cases",
                                            "Allele frequency in controls",
                                            "Imputation quality"))
    
    header_interp<-merge(header_interp, header_labels, by='Interpreted')
    header_interp$Keep<-factor(header_interp$Keep, levels=c('TRUE','FALSE'))
    header_interp<-header_interp[order(header_interp$Keep),]
    header_interp<-header_interp[,c('Original','Interpreted','Keep','Reason','Description')]
    
    return(header_interp)
  })
  
  output$header_check_messages.ui <- renderUI({
    header_interp<-head_interp()
    messages <- NULL
    if(!all(c('A1','A2') %in% header_interp$Interpreted)){
      messages <- c(messages, "Error: Missing A1 and A2 information.\n")
    }
    if(!any(c('OR','BETA','Z') %in% header_interp$Interpreted)){
      messages <- c(messages, "Error: Effect size (BETA, OR, Z) must be present.\n")
    }
    if(!('SNP' %in% header_interp$Interpreted) & !(all(c('CHR','BP') %in% header_interp$Interpreted))){
      messages <- c(messages, "Error: Either SNP, or CHR and BP, must be present.\n")
    }
    if(any(c('FRQ_A','FRQ_U') %in% header_interp$Interpreted) & sum(c('FRQ_A','FRQ_U') %in% header_interp$Interpreted) == 1){
      messages <- c(messages, "Warning: Both FRQ_A and FRQ_U must be present for either to be considered.\n")
    }
    
    # Display the messages
    if (!is.null(messages)) {
      tags$div(
        style = "color: red;",  # Style the message, e.g., red for errors
        messages
      )
    }
  })
  
  output$specify_n.ui <- renderUI({
    header_interp<-head_interp()
    tag_list<-NULL
    if(!(any(c('N','NEF') %in% header_interp$Interpreted)) & !(all(c('N_CAS','N_CON') %in% header_interp$Interpreted))){
      tag_list<-c(tag_list, tagList(
        h6("**Note.** N, NEF, or N_CAS and N_CON, were not present."),
        br(),
        h5("Specify sample size:"),
        numericInput("sample_size", NULL, value = NULL),
      ))
    }
    if((!('NEF' %in% header_interp$Interpreted)) &
       (!all(c('N_CAS','N_CON') %in% header_interp$Interpreted))){
      tag_list<-c(tag_list, tagList(
        br(),
        h5("If binary outcome, specify proportion of cases:"),
        numericInput("samp_prop", NULL, value = NULL)
      ))
    }
    tag_list
  })
  
  output$specify_pop_prev.ui <- renderUI({
    header_interp<-head_interp()
    if(input$ldsc_logical){
      tagList(
        h6("If binary outcome, specify the population prevalence:"),
        p('Note. This is used to convert heritability estimates onto the liability scale.'),
        numericInput("pop_prev", NULL, value = NULL),
        br(),
      )
    }
  })
  
  output$specify_population.ui <- renderUI({
    header_interp<-head_interp()
    
    tagList(
      br(),
      h5("Specify the population the GWAS based on:"),
      p('Note. If mixed, you can select the largest population within the GWAS.'),
      selectInput("population", NULL, selected = NULL, choices=populations),
      br()    
    )
  })
  
  output$sub_ss_head <- renderDataTable({
    tmp<-head_interp()
    
    datatable(tmp,
              options = list(
                columnDefs = list(list(
                  className = 'dt-center', targets = '_all'
                ))
              ),
              selection = 'none', 
              rownames = F)
  })
  
  output$submit_colnames.ui<-renderUI({
    req(input$sumstats)
    tagList(
      h5('Column name interpretation:'),
      fluidRow(
        column(width=6,
               dataTableOutput("sub_ss_head")
        )
      )
    )
  })
  
  # Check whether criteria are met for the button to be clicked
  submit_criteria<-reactive({

    header_interp<-head_interp()
    
    path <- input$sumstats[["datapath"]]
    email <- input$email
    
    error<-NULL
    if(email == ''){
      error<-c(error,"Error: No email address has been specified.\n")
    } else {
      if(!is_valid_email_format(email)){
        error<-c(error,"Error: Email address format is invalid\n")
      }
    }
    if(!all(c('A1','A2') %in% header_interp$Interpreted)){
      error<-c(error,"Error: Missing A1 and A2 information.\n\n")
    }
    if(!any(c('OR','BETA','Z') %in% header_interp$Interpreted)){
      error<-c(error,"Error: Effect size (BETA, OR, Z) must be present.\n")
    }
    if(!('SNP' %in% header_interp$Interpreted) & !(all(c('CHR','BP') %in% header_interp$Interpreted))){
      error<-c(error,"Error: Either SNP, or CHR and BP, must be present.\n")
    }
    if(
      !(any(c('N','NEF') %in% header_interp$Interpreted)) & 
      !(all(c('N_CAS','N_CON') %in% header_interp$Interpreted)) &
      is.na(input$sample_size)){
      error<-c(error,"Error: Sample size must be specified if N, NEFF, or N_CAS and N_CON, are not present.\n")
    }
    if(
      'TWAS-GSEA' %in% input$enrich_method & 
      is.null(input$twas_panels)
    ){
      error<-c(error,"Error: At least one TWAS panel must be selected to use TWAS-GSEA.\n")
    }
    if(
      (is.null(input$enrich_method) &
       !is.null(input$enrich_data)) | 
      (!is.null(input$enrich_method) &
       is.null(input$enrich_data))
    ){
      error<-c(error,"Error: For enrichment analysis, the user must select both an enrichment method and dataset.\n")
    }
    if(
      input$clump_logical == F &
      input$cojo_logical == F &
      input$finemap_logical == F &
      input$magma_gene_logical == F &
      is.null(input$enrich_method) &
      is.null(input$enrich_data) &
      is.null(input$twas_panels) &
      is.null(input$pwas_panels) &
      is.null(input$smr_expr_panels) &
      is.null(input$smr_protein_panels) &
      is.null(input$smr_protein_panels)
    ){
      error<-c(error,"Error: No analyses have been selected. Go to the 'Options' tab.\n")
    }
    if(
      is.null(input$population)
    ){
      error<-c(error,"Error: You must specify which population the GWAS is based on.\n")
    }
    
    return(error)
  })
  
  # Reactive value to track submission state
  fileSubmitted <- reactiveVal(FALSE)
  
  # Reactive value to track messages
  message <- reactiveVal("")
  
  observeEvent(input$submit_button, {
    if(is.null(submit_criteria())){
      message('You will receive a confirmation email shortly.')
    } else {
      message(submit_criteria())
    }
    
    if(is.null(submit_criteria())){
      # Create directory to store output
      dir.create(paste0("../uploads/submissions"), recursive = T)
      dir.create(paste0("../uploads/outputs"), recursive = T)
            
      # Assign job id
      if(!file.exists('../uploads/job_list.txt')){
        job_id<-1
        write.table(job_id,'../uploads/job_list.txt', col.names=F, row.names=F, quote=F)
      } else {
        job_list<-as.numeric(fread('../uploads/job_list.txt', header=F)$V1)
        job_list<-job_list[order(job_list)]
        job_id<-max(job_list)+1
        job_list<-c(job_list, job_id)
        write.table(job_list,'../uploads/job_list.txt', col.names=F, row.names=F, quote=F)
      }
      
      # Create directory to store output
      dir.create(paste0("../uploads/submissions/job_",job_id))
      dir.create(paste0("../uploads/outputs/job_",job_id))
      
      # Copy the file to openstack
      file.copy(input$sumstats$datapath, paste0("../uploads/submissions/job_",job_id,"/job_",job_id,"_ss.gz"))
      
      #############
      # Prepare config file based on user input
      #############
      config_file<-NULL
      config_file<-c(paste0("outdir: ../uploads/outputs/job_",job_id))
      config_file<-c("config_file: ../uploads/submissions/job_1/config.yaml")
      config_file<-c("rosmap_fusion: /scratch/prj/oliverpainfel/Data/ROSMAP-PWAS/ROSMAP.n376.fusion.WEIGHTS.zip")
      config_file<-c("banner_fusion: /scratch/prj/oliverpainfel/Data/Banner-PWAS/Banner.n152.fusion.WEIGHTS.zip")
      config_file<-c("rosmap_smr: /scratch/prj/oliverpainfel/Data/ROSMAP_pQTL/ROSMAP.n376.pQTL.txt")
      
      if(input$ldsc_logical){
        config_file<-c(config_file, "ldsc: T")
      } else {
        config_file<-c(config_file, "ldsc: F")
      }
      
      if(input$clump_logical){
        config_file<-c(config_file, "clump: T")
      } else {
        config_file<-c(config_file, "clump: F")
      }
      
      if(input$cojo_logical){
        config_file<-c(config_file, "cojo: T")
      } else {
        config_file<-c(config_file, "cojo: F")
      }
      
      if(input$finemap_logical){
        config_file<-c(config_file, "finemap: T")
      } else {
        config_file<-c(config_file, "finemap: F")
      }
      
      if(input$magma_gene_logical){
        config_file<-c(config_file, "magma_gene: T")
      } else {
        config_file<-c(config_file, "magma_gene: F")
      }
      
      if('MAGMA' %in% input$enrich_method & 'DrugTargetor - Drug-Gene Interactions' %in% input$enrich_data){
        config_file<-c(config_file, "magma_drugtargetor: T")
        config_file<-gsub('magma_gene: F','magma_gene: T', config_file)
      } else {
        config_file<-c(config_file, "magma_drugtargetor: F")
      }
      
      if('TWAS-GSEA' %in% input$enrich_method & 'DrugTargetor - Drug-Gene Interactions' %in% input$enrich_data){
        config_file<-c(config_file, "twas_gsea_drugtargetor: T")
      } else {
        config_file<-c(config_file, "twas_gsea_drugtargetor: F")
      }
      
      config_file<-c(config_file, "external_weights: F")
      config_file<-c(config_file, "external_weights_pos_path: []")
      config_file<-c(config_file, "twas_conditional: F")
      
      if(any(!(input$twas_panels == 'psychencode'))){
        config_file<-c(config_file, "twas_panel_psychencode: T")
      } else {
        config_file<-c(config_file, "twas_panel_psychencode: F")
      }
      
      if(any(fusion_twas_panel_names$original %in% input$twas_panels)){
        config_file<-c(config_file, "twas_panel_fusion: T")
        config_file<-c(config_file, paste0("gtex_weights: [\"",paste(input$twas_panels[input$twas_panels %in% gtex_fusion_panels$original], collapse="\",\""),"\"]"))
        config_file<-c(config_file, paste0("non_gtex_weights: [\"",paste(input$twas_panels[input$twas_panels %in% nongtex_fusion_panels$original], collapse="\",\""),"\"]"))
      } else {
        config_file<-c(config_file, "twas_panel_fusion: F")
        config_file<-c(config_file, "gtex_weights: []")
        config_file<-c(config_file, "non_gtex_weights: []")
      }
      
      if('smr_expression_panel_psychencode' %in% input$smr_expr_panels){
        config_file<-c(config_file, "smr_expression_panel_psychencode: T")
      } else {
        config_file<-c(config_file, "smr_expression_panel_psychencode: F")
      }
      
      if('smr_expression_panel_metabrain_basalganglia' %in% input$smr_expr_panels){
        config_file<-c(config_file, "smr_expression_panel_metabrain_basalganglia: T")
      } else {
        config_file<-c(config_file, "smr_expression_panel_metabrain_basalganglia: F")
      }
      
      if('smr_expression_panel_metabrain_cerebellum' %in% input$smr_expr_panels){
        config_file<-c(config_file, "smr_expression_panel_metabrain_cerebellum: T")
      } else {
        config_file<-c(config_file, "smr_expression_panel_metabrain_cerebellum: F")
      }
      
      if('smr_expression_panel_metabrain_cortex' %in% input$smr_expr_panels){
        config_file<-c(config_file, "smr_expression_panel_metabrain_cortex: T")
      } else {
        config_file<-c(config_file, "smr_expression_panel_metabrain_cortex: F")
      }
      
      if('smr_expression_panel_metabrain_hippocampus' %in% input$smr_expr_panels){
        config_file<-c(config_file, "smr_expression_panel_metabrain_hippocampus: T")
      } else {
        config_file<-c(config_file, "smr_expression_panel_metabrain_hippocampus: F")
      }
      
      if('smr_expression_panel_metabrain_spinalcord' %in% input$smr_expr_panels){
        config_file<-c(config_file, "smr_expression_panel_metabrain_spinalcord: T")
      } else {
        config_file<-c(config_file, "smr_expression_panel_metabrain_spinalcord: F")
      }
      
      if('smr_expression_panel_eqtlgen' %in% input$smr_expr_panels){
        config_file<-c(config_file, "smr_expression_panel_eqtlgen: T")
      } else {
        config_file<-c(config_file, "smr_expression_panel_eqtlgen: F")
      }
      
      if('rosmap' %in% input$pwas_panels){
        config_file<-c(config_file, "pwas_panel_rosmap: T")
      } else {
        config_file<-c(config_file, "pwas_panel_rosmap: F")
      }
      
      if('banner' %in% input$pwas_panels){
        config_file<-c(config_file, "pwas_panel_banner: T")
      } else {
        config_file<-c(config_file, "pwas_panel_banner: F")
      }
      
      if('smr_protein_panel_rosmap' %in% input$smr_protein_panels){
        config_file<-c(config_file, "smr_protein_panel_rosmap: T")
      } else {
        config_file<-c(config_file, "smr_protein_panel_rosmap: F")
      }
      
      config_file<-c(config_file, "gcsc: F")
      config_file<-c(config_file, "gcsc_tissues: []")
      
      writeLines(config_file, paste0("../uploads/submissions/job_",job_id,"/config.yaml"))
      
      ############
      # Create gwas_list
      ############
      
      gwas_list<-data.frame(
        name='sub',
        path=paste0('../uploads/submissions/job_',job_id,'/job_',job_id,'_ss.gz'),
        population=input$population,
        sampling=ifelse(is.null(input$samp_prop), NA, input$samp_prop),
        prevelance=ifelse(is.null(input$pop_prev), NA, input$pop_prev),
        mean=NA,
        sd=NA,
        label="\"Outcome\""
      )
      write.table(gwas_list, paste0("../uploads/submissions/job_",job_id,"/gwas_list.txt"), col.names=T, row.names=F, quote=F)
      
      ############
      # Save submission data
      ############
      
      write.table(
        input$email,
        paste0("../uploads/submissions/job_",job_id,"/job_", job_id,"_email.txt"),
        col.names = F,
        row.names = F,
        quote = F
      )
      
      fileSubmitted(TRUE)
      
      # Transfer submission to HPC
      # This can be done by sending these output to a mounted HPC folder.
      
      # Trigger analysis on HPC
      # This will depend on your HPC setup; you might run a command like:
      # system("ssh user@hpc 'bash run_analysis_script.sh'")
    }
    
  })
  
  observe({
    if (is.null(input$sumstats)) {
      shinyjs::addClass(selector = "#submit_button", class = "disabled-button")
    } else {
      if(!fileSubmitted()){
        shinyjs::removeClass(selector = "#submit_button", class = "disabled-button")
      } else {
        shinyjs::hide("submit_button")
      }
    }
  })
  
  output$submit_text <- renderUI({
    HTML(gsub("\n", "<br>", message()))
  })
  
  ##############################################################################
  # Explore
  ##############################################################################
  
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
  
  observeEvent(input$gwas_selector, {
    
    ######
    # Read in configuration data
    ######
    
    # Read in config file
    config<-gwas_data()$configuration$config
    
    # Read in GWAS list
    gwas_list<-gwas_data()$configuration$gwas_list
    
    # Determine which analyses were requested
    clump_logical<-config[grepl('clump:',config)] == "clump: T"
    
    cojo_logical<-config[grepl('cojo:',config)] == "cojo: T"
    
    finemap_logical<-config[grepl('finemap:',config)] == "finemap: T"
    
    ldsc_logical<-config[grepl('ldsc:',config)] == "ldsc: T"
    
    magma_gene_logical<-config[grepl('magma_gene:',config)] == "magma_gene: T"
    
    twas_panel_psychencode_logical<-config[grepl('twas_panel_psychencode:',config)] == "twas_panel_psychencode: T"
    twas_panel_fusion_logical<-config[grepl('twas_panel_fusion:',config)] == "twas_panel_fusion: T"
    
    twas_logical<-any(twas_panel_psychencode_logical, twas_panel_fusion_logical)
    
    twas_conditional_logical<-config[grepl('twas_conditional:',config)] == "twas_conditional: T" & twas_logical == T
    
    smr_expression_panel_psychencode_logical<-config[grepl('smr_expression_panel_psychencode:',config)] == "smr_expression_panel_psychencode: T"
    
    smr_expression_panel_metabrain_basalganglia_logical<-config[grepl('smr_expression_panel_metabrain_basalganglia:',config)] == "smr_expression_panel_metabrain_basalganglia: T"
    smr_expression_panel_metabrain_cerebellum_logical<-config[grepl('smr_expression_panel_metabrain_cerebellum:',config)] == "smr_expression_panel_metabrain_cerebellum: T"
    smr_expression_panel_metabrain_cortex_logical<-config[grepl('smr_expression_panel_metabrain_cortex:',config)] == "smr_expression_panel_metabrain_cortex: T"
    smr_expression_panel_metabrain_hippocampus_logical<-config[grepl('smr_expression_panel_metabrain_hippocampus:',config)] == "smr_expression_panel_metabrain_hippocampus: T"
    smr_expression_panel_metabrain_spinalcord_logical<-config[grepl('smr_expression_panel_metabrain_spinalcord:',config)] == "smr_expression_panel_metabrain_spinalcord: T"
    
    smr_expression_panel_eqtlgen_logical<-config[grepl('smr_expression_panel_eqtlgen:',config)] == "smr_expression_panel_eqtlgen: T"
    
    pwas_panel_rosmap_logical<-config[grepl('pwas_panel_rosmap:',config)] == "pwas_panel_rosmap: T"
    pwas_panel_banner_logical<-config[grepl('pwas_panel_banner:',config)] == "pwas_panel_banner: T"
    
    smr_protein_panel_rosmap_logical<-config[grepl('smr_protein_panel_rosmap:',config)] == "smr_protein_panel_rosmap: T"
    
    magma_drugtargetor_logical<-config[grepl('magma_drugtargetor:',config)] == "magma_drugtargetor: T"
    twas_gsea_lincs_logical<-config[grepl('twas_gsea_lincs:',config)] == "twas_gsea_lincs: T"
    twas_so_lincs_logical<-config[grepl('twas_so_lincs:',config)] == "twas_so_lincs: T"
    
    dgi_db_comp_logical<-config[grepl('dgi_db_comp:',config)] == "dgi_db_comp: T"
    
    twas_gsea_drugtargetor_logical<-config[grepl('twas_gsea_drugtargetor:',config)] == "twas_gsea_drugtargetor: T"
    
    # Summarise this information
    mol_assoc_logical<-any(magma_gene_logical,
                           twas_logical,
                           smr_expression_panel_psychencode_logical,
                           smr_expression_panel_metabrain_basalganglia_logical,
                           smr_expression_panel_metabrain_cerebellum_logical,
                           smr_expression_panel_metabrain_cortex_logical,
                           smr_expression_panel_metabrain_hippocampus_logical,
                           smr_expression_panel_metabrain_spinalcord_logical,
                           smr_expression_panel_eqtlgen_logical,
                           pwas_panel_rosmap_logical,
                           pwas_panel_banner_logical,
                           smr_protein_panel_rosmap_logical)
    
    metabrain_logical<-any(smr_expression_panel_metabrain_basalganglia_logical,
                           smr_expression_panel_metabrain_cerebellum_logical,
                           smr_expression_panel_metabrain_cortex_logical,
                           smr_expression_panel_metabrain_hippocampus_logical,
                           smr_expression_panel_metabrain_spinalcord_logical)
    
    smr_expression_logical<-any(smr_expression_panel_psychencode_logical,
                                metabrain_logical,
                                smr_expression_panel_eqtlgen_logical)
    
    expression_logical<-any( twas_logical,
                             smr_expression_logical)
    
    psychencode_logical<-any(twas_panel_psychencode_logical,
                             smr_expression_panel_psychencode_logical)
    
    protein_logical<-any( pwas_panel_rosmap_logical,
                          pwas_panel_banner_logical,
                          smr_protein_panel_rosmap_logical)
    
    drug_logical<-any(dgi_db_comp_logical, 
                      magma_drugtargetor_logical,
                      twas_gsea_lincs_logical,
                      twas_so_lincs_logical,
                      twas_gsea_drugtargetor_logical)
    
    # Assign select GWAS name to variable
    gwas_selected <- input$gwas_selector
    
    ######
    # Prepare table for GWAS QC tab
    ######
    
    # Create a table showing key statistics
    qc_val<-data.table(name=gwas_selected,
                       label=gwas_list$label[gwas_list$name == gwas_selected],
                       n_var_orig=gwas_data()[[gwas_selected]]$gwas_qc$cleaner_dat$val$n_var_orig,
                       build=gwas_data()[[gwas_selected]]$gwas_qc$cleaner_dat$val$build$build,
                       n_snp_final=gwas_data()[[gwas_selected]]$gwas_qc$cleaner_dat$val$n_snp_final,
                       lambda_gc=gwas_data()[[gwas_selected]]$gwas_qc$focus_dat$val$lambda_gc,
                       max_chi2=gwas_data()[[gwas_selected]]$gwas_qc$focus_dat$val$max_chi2,
                       n_sig_snp=gwas_data()[[gwas_selected]]$gwas_qc$focus_dat$val$n_sig_snp,
                       obs_h2=paste0(round(gwas_data()[[gwas_selected]]$gwas_qc$ldsc_dat$val$obs_h2_est,3), " (",round(gwas_data()[[gwas_selected]]$gwas_qc$ldsc_dat$val$obs_h2_se,3),")"),
                       int=paste0(round(gwas_data()[[gwas_selected]]$gwas_qc$ldsc_dat$val$int_est,3), " (",round(gwas_data()[[gwas_selected]]$gwas_qc$ldsc_dat$val$int_se,3),")"))
    
    names(qc_val)<-c('GWAS Name',
                     'GWAS Label',
                     'N variants pre-QC',
                     'Identified genome build',
                     'N variants post-QC',
                     'Lambda GC',
                     'Max. chi^2',
                     'N genome-wide significant variants',
                     "LDSC SNP-heritability (SE; observed scale)",
                     "LDSC intercept (SE)")
    
    qc_val<-t(qc_val)
    qc_val<-data.table(Parameter=dimnames(qc_val)[[1]],
                       Value=qc_val[,1])
    
    snp_assoc_data <- gwas_data()[[gwas_selected]]$snp_assoc$clump
    
    output$qc_table <- renderDataTable({
      datatable(qc_val, options = list(dom = 't', 
                                       ordering=F), 
                selection = 'none', 
                rownames = F, 
                colnames = '') %>%
        formatStyle(columns = c("Parameter"), fontWeight = 'bold', textAlign = "center") %>%
        formatStyle(columns = c("Value"), textAlign = "center")
    })
    
    ########
    # Prepare data for SNP Association tab
    ########
    
    snp_assoc_lead_data <- reactive({
      if (input$clumping_type == "ld_clumping") {
        # Read in LD-clumped SNP-associations
        snp_assoc_lead <- gwas_data()[[gwas_selected]]$snp_assoc$clump
      }
      
      if (input$clumping_type == "cojo_analysis") {
        # Read in the COJO associations
        snp_assoc_lead <- gwas_data()[[gwas_selected]]$snp_assoc$cojo
      }
      
      snp_assoc_lead<-snp_assoc_lead[,names(snp_assoc_lead) %in% c("CHR","BP","SNP","A1","A2","BETA","SE","P","NearestGene"),with=F]
      
      # error on screen non numeric input
      snp_assoc_lead$P <- as.numeric(snp_assoc_lead$P)
      
      # Subset SNP data by p-value and check whether it is not
      snp_assoc_lead <- snp_assoc_lead[snp_assoc_lead$P < as.numeric(input$pvalue_threshold), ]
      
      # Round BETA and SE
      snp_assoc_lead$BETA<-round(snp_assoc_lead$BETA, 2)
      snp_assoc_lead$SE<-round(snp_assoc_lead$SE, 2)
      
      # Sort by CHR and BP
      snp_assoc_lead<-snp_assoc_lead[order(snp_assoc_lead$CHR, snp_assoc_lead$BP),]
      
      return(snp_assoc_lead)
      
    })
    
    output$snp_assoc_lead_table <- renderDataTable({
      # Create java script to force scientific notation for P value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "}"
      )
      
      datatable(snp_assoc_lead_data(), rownames=F, width=7, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:7),
                          list(width = '60px', targets = 7),
                          list(width = '600px', targets = 8))))
    })
    
    snp_assoc_finemap_data <- reactive({
      if (input$l_param == "L1") {
        snp_assoc_finemap <- gwas_data()[[gwas_selected]]$snp_assoc$susie$L1
      }
      
      if (input$l_param == "L10") {
        snp_assoc_finemap <- gwas_data()[[gwas_selected]]$snp_assoc$susie$L10
      }
      
      snp_assoc_finemap$cs_log10bf<-NULL
      snp_assoc_finemap$cs_avg_r2<-NULL
      snp_assoc_finemap$cs_min_r2<-NULL
      snp_assoc_finemap$TopPIP<-round(snp_assoc_finemap$TopPIP,2)
      
      return(snp_assoc_finemap)
      
    })
    
    output$snp_assoc_finemap_table <- renderDataTable({
      datatable(snp_assoc_finemap_data(), rownames=F, width=7, options = list(
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:5))))
    })
    
    #######
    # Prepare data for molecular association tables
    #######
    
    output$mol_assoc_magma_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$mol_assoc$magma
      
      datatable(tmp, rownames=F, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:5))))
    })
    
    output$mol_assoc_fusion_expr_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res
      tmp$TWAS.Z<-round(tmp$TWAS.Z, 3)
      tmp$COLOC_logical<-NULL
      tmp[is.na(tmp)]<-'NA'
      
      names(tmp)[names(tmp) == 'TWAS.Z']<-'Z'
      names(tmp)[names(tmp) == 'TWAS.P']<-'P'
      names(tmp)[names(tmp) == 'TWAS.P.FDR']<-'P.FDR'
      
      tmp<-tmp[, c("PANEL","CHR","P0","P1","Ensembl ID","Gene Symbol","Z","P","P.FDR","COLOC.PP3","COLOC.PP4"), with=F]
      
      datatable(tmp, rownames=F, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:10))))
    })
    
    output$mol_assoc_fusion_protein_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res
      tmp$pwas_all.Z<-round(tmp$pwas_all.Z, 3)
      tmp$COLOC_logical<-NULL
      tmp[is.na(tmp)]<-'NA'
      
      names(tmp)[names(tmp) == 'pwas_all.Z']<-'Z'
      names(tmp)[names(tmp) == 'pwas_all.P']<-'P'
      names(tmp)[names(tmp) == 'pwas_all.P.FDR']<-'P.FDR'
      
      tmp<-tmp[, c("PANEL","CHR","P0","P1","Ensembl ID","Gene Symbol","Z","P","P.FDR","COLOC.PP3","COLOC.PP4"), with=F]
      
      datatable(tmp, rownames=F, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:10))))
    })
    
    output$mol_assoc_smr_expr_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "  var x = data[9];",
        "  $('td:eq(9)', row).html(x.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res
      tmp<-tmp[,c("PANEL","CHR","BP","Ensembl ID","Gene Symbol","b_SMR","se_SMR","p_SMR","p_SMR.FDR","p_HEIDI"), with=F]
      tmp$b_SMR<-round(tmp$b_SMR, 3)
      tmp$se_SMR<-round(tmp$b_SMR, 3)
      
      names(tmp)[names(tmp) == 'b_SMR']<-'BETA'
      names(tmp)[names(tmp) == 'se_SMR']<-'SE'
      names(tmp)[names(tmp) == 'p_SMR']<-'P'
      names(tmp)[names(tmp) == 'p_SMR.FDR']<-'P.FDR'
      names(tmp)[names(tmp) == 'p_HEIDI']<-"P (HEIDI)"
      
      datatable(tmp, rownames=F, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:9))))
    })
    
    output$mol_assoc_smr_protein_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "  var y = data[9];",
        "  $('td:eq(9)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res
      tmp<-tmp[,c("PANEL","CHR","BP","Ensembl ID","Gene Symbol","b_SMR","se_SMR","p_SMR","p_SMR.FDR","p_HEIDI"), with=F]
      tmp$b_SMR<-round(tmp$b_SMR, 3)
      tmp$se_SMR<-round(tmp$b_SMR, 3)
      
      names(tmp)[names(tmp) == 'b_SMR']<-'BETA'
      names(tmp)[names(tmp) == 'se_SMR']<-'SE'
      names(tmp)[names(tmp) == 'p_SMR']<-'P'
      names(tmp)[names(tmp) == 'p_SMR.FDR']<-'P.FDR'
      names(tmp)[names(tmp) == 'p_HEIDI']<-"P (HEIDI)"
      
      datatable(tmp, rownames=F, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:9))))
    })
    
    #######
    # Prepare data for molecular association plot
    #######
    
    mol_assoc_summary_data <- reactive({
      
      all_func_res<-NULL
      
      if(finemap_logical){
        
        finemap_L1_tmp<-data.frame(Panel = "SuSie (L=1)",
                                   ID=gwas_data()[[gwas_selected]]$mol_assoc$finemap$L1,
                                   Z=1,
                                   Sig=F,
                                   Coloc=F,
                                   Method="SNP\nFine-mapping",
                                   Type='')
        
        all_func_res<-rbind(all_func_res, finemap_L1_tmp)
        
      }
      
      if(twas_logical){
        
        twas_tmp<-data.table(Panel=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$PANEL,
                             ID=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$`Gene Symbol`,
                             Z=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$TWAS.Z,
                             Sig=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$TWAS.P.FDR < 0.05,
                             Coloc=gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res$COLOC_logical)
        
        twas_tmp$Method<-'FUSION'
        twas_tmp$Type<-'Expr.'
        twas_tmp$Type[grepl('SPLIC',twas_tmp$Panel)]<-'Splice'
        
        # Retain only the most significant assoc for each gene within PANEL (only relevent for splice panel)
        twas_tmp<-twas_tmp[order(-abs(twas_tmp$Z)),]
        twas_tmp<-twas_tmp[!duplicated(paste0(twas_tmp$Panel, twas_tmp$ID)),]
        
        all_func_res<-rbind(all_func_res, twas_tmp)
      }
      
      if(smr_expression_logical){
        # SMR expression
        smr_expression_res_tmp<-data.table(Panel=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$PANEL,
                                           ID=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$`Gene Symbol`,
                                           Z=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$b_SMR/gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$se_SMR,
                                           Sig=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$p_SMR.FDR < 0.05,
                                           Coloc=gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res$p_HEIDI > 0.05)
        
        smr_expression_res_tmp$Method<-'SMR'
        smr_expression_res_tmp$Type<-'Expr.'
        
        all_func_res<-rbind(all_func_res, smr_expression_res_tmp)
      }
      
      # PWAS
      if(pwas_panel_rosmap_logical & !(pwas_panel_rosmap_logical & pwas_panel_banner_logical)){
        
        pwas_tmp<-data.table( Panel=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$PANEL,
                              ID=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$`Gene Symbol`,
                              Z=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$pwas_all.Z,
                              Sig=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$pwas_all.P.FDR < 0.05,
                              Coloc=gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res$COLOC_logical)
        
        pwas_tmp<-pwas_tmp[order(-abs(pwas_tmp$Z)),]
        pwas_tmp<-pwas_tmp[!duplicated(paste0(pwas_tmp$Panel, pwas_tmp$ID)),]
        pwas_tmp$Method<-'FUSION'
        pwas_tmp$Type<-'Protein'
        
        all_func_res<-rbind(all_func_res, pwas_tmp)
      }
      
      if(smr_protein_panel_rosmap_logical){
        
        # SMR protein
        smr_protein_res_tmp<-data.table(Panel=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$PANEL,
                                        ID=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$`Gene Symbol`,
                                        Z=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$b_SMR/gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$se_SMR,
                                        Sig=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$p_SMR.FDR < 0.05,
                                        Coloc=gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res$p_HEIDI > 0.05)
        
        smr_protein_res_tmp<-smr_protein_res_tmp[order(-abs(smr_protein_res_tmp$Z)),]
        smr_protein_res_tmp<-smr_protein_res_tmp[!duplicated(paste0(smr_protein_res_tmp$Panel, smr_protein_res_tmp$ID)),]
        smr_protein_res_tmp$Method<-'SMR'
        smr_protein_res_tmp$Type<-'Protein'
        
        all_func_res<-rbind(all_func_res, smr_protein_res_tmp)
        
      }
      
      if(magma_gene_logical){
        
        magma_tmp<-data.frame(Panel = 'MAGMA',
                              ID=gwas_data()[[gwas_selected]]$mol_assoc$magma$ID,
                              Z=abs(qnorm(as.numeric(gwas_data()[[gwas_selected]]$mol_assoc$magma$P))),
                              Sig=as.numeric(gwas_data()[[gwas_selected]]$mol_assoc$magma$P.FDR) < 0.05,
                              Coloc=F,
                              Method='MAGMA',
                              Type='')
        
        all_func_res<-rbind(all_func_res, magma_tmp)
        
      }
      
      if(clump_logical){
        
        nearest_tmp<-data.frame(Panel = 'NearestGene',
                                ID=gwas_data()[[gwas_selected]]$mol_assoc$nearest$clump,
                                Z=1,
                                Sig=F,
                                Coloc=F,
                                Method='Nearest\nGene',
                                Type='')
        
        all_func_res<-rbind(all_func_res, nearest_tmp)
        
      }
      
      return(all_func_res)
    })
    
    observe({
      all_func_res<-mol_assoc_summary_data()
      methods<-unique(all_func_res$Method)
      updateSelectInput(session, "selected_methods_mol", choices = methods, selected=methods)
    })
    
    observe({
      all_func_res<-mol_assoc_summary_data()
      expr_panels<-unique(all_func_res$Panel[all_func_res$Type == 'Expr.'])
      updateSelectInput(session, "selected_expr_panels_mol", choices = expr_panels, selected=expr_panels)
    })
    
    observe({
      all_func_res<-mol_assoc_summary_data()
      protein_panels<-unique(all_func_res$Panel[all_func_res$Type == 'Protein'])
      updateSelectInput(session, "selected_protein_panels_mol", choices = protein_panels, selected=protein_panels)
    })
    
    mol_assoc_summary_data_filtered<-reactive({
      all_func_res<-mol_assoc_summary_data()
      
      # Filter results table by user specified methods
      all_func_res<-all_func_res[all_func_res$Method %in% input$selected_methods_mol,]
      
      # Filter results table by user specified expression and protein panels
      if(any(all_func_res$Type == 'Expr.')){
        all_func_res<-all_func_res[!(all_func_res$Type == 'Expr.' & !(all_func_res$Panel %in% input$selected_expr_panels_mol)),]
      }
      if(any(all_func_res$Type == 'Protein')){
        all_func_res<-all_func_res[!(all_func_res$Type == 'Protein' & !(all_func_res$Panel %in% input$selected_protein_panels_mol)),]
      }
      
      # Filter results table if user specifies high confidence genes only
      if(input$conf_only_mol){
        all_func_res<-all_func_res[all_func_res$ID %in% all_func_res$ID[which((all_func_res$Sig == T & all_func_res$Coloc == T) | all_func_res$Panel == "SuSie (L=1)")],]
      }
      
      input_genes <- unlist(strsplit(input$geneInput_mol, "[, ]"))
      selected_genes <- input_genes[input_genes != ""]
      
      # Insert NA rows for all panels and methods so when filtering by gene, all selected panels and methods remain
      na_rows<-all_func_res[!(duplicated(paste0(all_func_res$Panel, all_func_res$Method))),]
      na_rows$ID<-NA
      na_rows$Z<-NA
      na_rows$Sig<-NA
      na_rows$Coloc<-NA
      
      all_func_res<-rbind(na_rows, all_func_res)
      
      if(length(selected_genes) > 0){
        if(sum(selected_genes %in% all_func_res$ID) > 0){
          all_func_res<-all_func_res[all_func_res$ID %in% selected_genes | is.na(all_func_res$ID),]
        } else {
          all_func_res<-data.frame(matrix(nrow=0, ncol=5))
        }
      }
      
      return(all_func_res)
    })
    
    # Identify number of genes
    plot_dim_mol<-reactive({
      all_func_res<-mol_assoc_summary_data_filtered()
      
      if(nrow(all_func_res) > 0){
        num_row <- length(unique(all_func_res$ID))
        plot_height<-(max(nchar(all_func_res$Panel))*3)+(num_row * 20)+100
        num_col <- length(unique(paste0(all_func_res$Panel,'_',all_func_res$Method,'_',all_func_res$Type)))
        plot_width<-(max(nchar(all_func_res$ID))*1.2)+(num_col * 60)
        plot_width<-max(plot_width,(length(unique(all_func_res$Method))*100))
      } else {
        plot_height<-100
        plot_width<-100
      }
      
      return(list(height=plot_height,
                  width=plot_width))
    })
    
    output$mol_assoc_plot<-renderPlot({
      
      all_func_res<-mol_assoc_summary_data_filtered()
      
      if(plot_dim_mol()[['height']] < 10000 & nrow(all_func_res) > 0){
        
        # Insert missing data
        all_func_res_all<-NULL
        for(i in unique(all_func_res$Panel)){
          for(j in unique(all_func_res$Method[all_func_res$Panel == i])){
            
            all_func_res_panel<-all_func_res[all_func_res$Panel == i & all_func_res$Method == j,]
            all_func_res_other<-all_func_res[!(all_func_res$Panel %in% all_func_res_panel$Panel) & !(all_func_res$Method %in% all_func_res_panel$Method),]
            all_func_res_other<-all_func_res_other[!(all_func_res_other$ID %in% all_func_res_panel$ID),]
            all_func_res_other<-unique(all_func_res_other$ID)
            
            if(length(all_func_res_other) > 0){
              all_func_res_panel_missing<-data.frame(ID=all_func_res_other)
              
              all_func_res_panel_missing$Panel=i
              all_func_res_panel_missing$ID=all_func_res_other
              all_func_res_panel_missing$Z=NA
              all_func_res_panel_missing$Sig=0
              all_func_res_panel_missing$Coloc=0
              all_func_res_panel_missing$Method=j
              all_func_res_panel_missing$Type=all_func_res_panel$Type[1]
              all_func_res_panel_missing$Group=all_func_res_panel$Group[1]
              
              all_func_res_panel_missing<-all_func_res_panel_missing[,names(all_func_res_panel)]
              
              all_func_res_all<-rbind(all_func_res_all,all_func_res_panel_missing)
            }
            
            all_func_res_all<-rbind(all_func_res_all,all_func_res_panel)
          }
        }
        
        # Now remove the NA rows
        all_func_res_all<-all_func_res_all[!is.na(all_func_res_all$ID),]
        
        all_func_res_all$Group<-paste0(all_func_res_all$Method,'\n',all_func_res_all$Type )
        all_func_res_all$Group[all_func_res_all$Group == 'SNP\nFine-mapping\n']<-'SuSiE'
        all_func_res_all$Group[all_func_res_all$Group == 'MAGMA\n']<-'MAGMA'
        all_func_res_all$Group[all_func_res_all$Group == 'Nearest\nGene\n']<-'Nearest\nGene'
        all_func_res_all$Group<-factor(all_func_res_all$Group, levels=unique(all_func_res_all$Group))
        
        all_func_res_all<-all_func_res_all[order(as.character(all_func_res_all$ID)),]
        
        all_func_res_all$ID<-factor(all_func_res_all$ID, levels=unique(all_func_res_all$ID))
        
        group_siz<-NULL
        for(i in unique(all_func_res_all$Group)){
          group_siz<-rbind(group_siz, data.frame(Group=i,
                                                 Size=length(unique(all_func_res_all$Panel[all_func_res_all$Group==i]))))
        }
        
        # Set minimum size to 3 to allow space for labels
        group_siz$Size[group_siz$Size < 2]<-2
        group_siz$Prop<-group_siz$Size/sum(group_siz$Size)
        group_siz$Width<-4*group_siz$Prop
        
        all_func_res_all$ID<-factor(all_func_res_all$ID, levels=rev(unique(as.character(all_func_res_all$ID))))
        
        x<-c(-max(abs(all_func_res_all$Z), na.rm=T),0,max(abs(all_func_res_all$Z), na.rm=T))
        x<-(x-min(x))/(max(x)-min(x))
        
        # Create second version of the plot
        heatmap<-ggplot(data = all_func_res_all, aes(x = Panel, y = ID)) +
          theme_bw()	+
          geom_point(data=all_func_res_all[all_func_res_all$Sig == T,], aes(x = Panel, y = ID), colour='black', size=5) +
          geom_point(data=all_func_res_all[all_func_res_all$Coloc ==T & all_func_res_all$Sig == T,], aes(x = Panel, y = ID), colour='black', shape=15, size=6) +
          geom_point(data=all_func_res_all[(all_func_res_all$Method == 'SNP\nFine-mapping' | all_func_res_all$Method == 'Nearest\nGene') & !(is.na(all_func_res_all$Z)),], aes(x = Panel, y = ID), colour='#00FF00', size=5) +
          geom_point(data=all_func_res_all[all_func_res_all$Method != 'SNP\nFine-mapping' & all_func_res_all$Method != 'Nearest\nGene',], aes(colour = Z), size=4) +
          scale_colour_gradientn(colours=c("#0066FF","#0099FF","#FFFFFF","#FF6666","#FF0000"), na.value = NA,name = "Z-score", limits = c(-max(abs(all_func_res_all$Z), na.rm=T), max(abs(all_func_res_all$Z), na.rm=T)), values=x) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
          labs(x='', y='') +
          facet_wrap(~ Group , nrow=1, scales = "free_x") +
          scale_y_discrete(limits= unique(rev(all_func_res_all$ID))) +
          theme(text = element_text(size = 14))
        
        gt = ggplot_gtable(ggplot_build(heatmap))
        
        for(i in 1:nrow(group_siz)){
          gt$widths[gt$layout$l[grep(paste0('panel-',i,'-1'), gt$layout$name)]] = group_siz$Width[i]*gt$widths[gt$layout$l[grep(paste0('panel-',i,'-1'), gt$layout$name)]]
          
        }
        
        grid.draw(gt)
      } else {
        NULL
      }
    })
    
    output$mol_assoc_plot.ui <- renderUI({
      if(plot_dim_mol()[['height']] < 10000 & nrow(mol_assoc_summary_data_filtered()) > 0){
        plotOutput("mol_assoc_plot", height = plot_dim_mol()[['height']], width=plot_dim_mol()[['width']])
      } else {
        NULL
      }
    })
    
    output$message_too_large_mol <- renderUI({
      if(plot_dim_mol()[['height']] > 10000){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "Plot is too large. Restrict to high-confidence genes or specify a list of genes."
        ))
      }
    })
    
    output$message_no_genes_mol <- renderUI({
      if(nrow(mol_assoc_summary_data_filtered()) == 0){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "No genes are present."
        ))        
      }
    })
    
    #######
    # Prepare data for drug-specific association tables
    #######
    
    output$tx_drug_magma_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$tx$drug$magma
      tmp$BETA<-round(tmp$BETA,3)
      tmp$SE<-round(tmp$SE,3)
      
      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 4:5)
          )),
        escape = FALSE
      )
    })
    
    output$tx_drug_gcsc_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$tx$drug$gcsc
      tmp$Enrichment<-round(tmp$Enrichment, 3)
      tmp$SE<-round(tmp$SE, 3)
      tmp$Z<-round(tmp$Z, 3)
      
      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 4:5)
          )),
        escape = FALSE
      )
    })
    
    output$tx_drug_twas_gsea_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$tx$drug$twas_gsea
      tmp$Estimate<-round(tmp$Estimate, 3)
      tmp$SE<-round(tmp$SE, 3)
      
      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 4:5)
          )),
        escape = FALSE
      )
    })
    
    #######
    # Prepare data for drug enrichment plot
    #######
    
    tx_drug_summary_data <- reactive({
      
      ###
      # MAGMA
      ###
      magma_gs<-gwas_data()[[gwas_selected]]$tx$drug$magma
      
      # Convert one-sided p to a Z score
      magma_gs$Z<--qnorm(magma_gs$P)
      magma_gs<-magma_gs[,c('Name','Z','P','P.FDR','ATC Code')]
      magma_gs$Method<-'MAGMA'
      magma_gs$Panel<-'MAGMA'
      
      ###
      # GCSC
      ###
      gcsc_gs<-gwas_data()[[gwas_selected]]$tx$drug$gcsc
      
      gcsc_gs<-gcsc_gs[,c('Name','Z','P','P.FDR','ATC Code')]
      gcsc_gs$Method<-'GCSC'
      gcsc_gs$Panel<-'Brain and Blood'
      
      ###
      # TWAS-GSEA
      ###
      
      gsea_gs<-gwas_data()[[gwas_selected]]$tx$drug$twas_gsea
      gsea_gs$Method<-'TWAS-GSEA'
      gsea_gs$Z<-gsea_gs$Estimate/gsea_gs$SE
      gsea_gs<-gsea_gs[,c('Name','Z','P','P.FDR','Method','Panel','ATC Code')]
      # Flip Z so >0 indicates reversal of GWAS outcome.
      gsea_gs$Z<--gsea_gs$Z
      
      # Insert missing values
      gsea_gs_all<-gsea_gs
      for(i in unique(gsea_gs_all$Panel)){
        gsea_gs_i<-gsea_gs[gsea_gs$Panel == i,]
        gsea_gs_other<-gsea_gs[gsea_gs$Panel != i,]
        gsea_gs_rest<-data.frame(Name=unique(gsea_gs_other$Name[!(gsea_gs_other$Name %in% gsea_gs_i$Name)]),
                                 Z=NA,
                                 P=NA,
                                 P.FDR=NA,
                                 Method='TWAS-GSEA',
                                 Panel=i,
                                 ATC_Code=NA)
        names(gsea_gs_rest)<-gsub('ATC_Code','ATC Code',names(gsea_gs_rest))
        
        gsea_gs_all<-rbind(gsea_gs_all, gsea_gs_rest)
        
      }
      gsea_gs<-gsea_gs_all
      
      ###
      # Combine results
      ###
      
      all_gs<-do.call(rbind, list(magma_gs, gcsc_gs, gsea_gs))
      
      return(all_gs)
    })
    
    observe({
      all_gs<-tx_drug_summary_data()
      methods<-unique(all_gs$Method)
      updateSelectInput(session, "selected_methods_drug", choices = methods, selected=methods)
    })
    
    observe({
      all_gs<-tx_drug_summary_data()
      expr_panels<-unique(all_gs$Panel[all_gs$Method == 'TWAS-GSEA'])
      updateSelectInput(session, "selected_expr_panels_drug", choices = expr_panels, selected=expr_panels)
    })
    
    tx_drug_summary_data_filtered<-reactive({
      all_gs<-tx_drug_summary_data()
      
      # Filter results table by user specified methods
      all_gs<-all_gs[all_gs$Method %in% input$selected_methods_drug,]
      
      # Filter results table by user specified expression
      if(any(all_gs$Method == 'TWAS-GSEA')){
        all_gs<-all_gs[!(all_gs$Method == 'TWAS-GSEA' & !(all_gs$Panel %in% input$selected_expr_panels_drug)),]
      }
      
      # Filter results table if user specifies high confidence genes only
      if(input$conf_only_drug){
        all_gs<-all_gs[all_gs$Name %in% all_gs$Name[which(all_gs$P.FDR < 0.05)],]
      }
      
      input_drugs <- unlist(strsplit(input$drugInput_drug, "[, ]"))
      selected_drugs <- input_drugs[input_drugs != ""]
      
      input_atc <- unlist(strsplit(input$atcInput_drug, "[, ]"))
      selected_atc <- input_atc[input_atc != ""]
      
      # Insert NA rows for all panels and methods so when filtering by drug, all selected panels and methods remain
      na_rows<-all_gs[!(duplicated(paste0(all_gs$Panel, all_gs$Method))),]
      na_rows$Name<-NA
      na_rows$Z<-NA
      na_rows$P<-NA
      na_rows$P.FDR<-NA
      
      all_gs<-rbind(na_rows, all_gs)
      
      if(length(selected_drugs) > 0){
        if(sum(grepl(paste(selected_drugs, collapse='|'), all_gs$Name, ignore.case = T)) > 0){
          all_gs<-all_gs[grepl(paste(selected_drugs, collapse='|'), all_gs$Name, ignore.case = T) & !is.na(all_gs$Name),]
        } else {
          all_gs<-data.frame(matrix(nrow=0, ncol=5))
        }
      }
      
      if(length(selected_atc) > 0){
        if(sum(grepl(paste(selected_atc, collapse='|'), all_gs$`ATC Code`, ignore.case = T)) > 0){
          all_gs<-all_gs[grepl(paste(selected_atc, collapse='|'), all_gs$`ATC Code`, ignore.case = T) & !is.na(all_gs$`ATC Code`),]
        } else {
          all_gs<-data.frame(matrix(nrow=0, ncol=5))
        }
      }
      
      return(all_gs)
    })
    
    # Identify number of drugs
    plot_dim_drug<-reactive({
      all_gs<-tx_drug_summary_data_filtered()
      
      if(nrow(all_gs) > 0){
        num_row <- length(unique(all_gs$Name))
        plot_height<-(max(nchar(all_gs$Panel))*3)+(num_row * 20)+100
        num_col <- length(unique(paste0(all_gs$Panel,'_',all_gs$Method,'_')))
        plot_width<-150+(max(nchar(all_gs$Name), na.rm=T)*2)+(num_col * 50)
        plot_width<-max(plot_width,(length(unique(all_gs$Method))*100))
      } else {
        plot_height<-100
        plot_width<-100
      }
      
      return(list(height=plot_height,
                  width=plot_width))
    })
    
    observe({
      tmp<-tx_drug_summary_data_filtered()
      choices<-'All - Z'
      if(any(tmp$Method == 'MAGMA')){
        choices<-c(choices, 'MAGMA - Z')
      }
      if(any(tmp$Method == 'GCSC')){
        choices<-c(choices, 'GCSC - Z')
      }
      if(any(tmp$Method == 'TWAS.GSEA')){
        choices<-c(choices, 'TWAS.GSEA - Z')
      }
      if(length(unique(tmp$Method)) == 1){
        choices<-choices[choices != 'All - Z']
      }
      choices<-c(choices, 'Alphabetical')
      
      updateSelectInput(session, "selected_sort_drug", choices = choices, selected=choices[1])
    })
    
    output$tx_drug_plot<-renderPlot({
      
      all_gs<-tx_drug_summary_data_filtered()
      all_gs$`ATC Code`<-NULL
      
      if(plot_dim_drug()[['height']] < 10000 & nrow(all_gs) > 0){
        
        # Insert missing data
        all_gs_all<-NULL
        for(i in unique(all_gs$Panel)){
          for(j in unique(all_gs$Method[all_gs$Panel == i])){
            
            all_gs_panel<-all_gs[all_gs$Panel == i & all_gs$Method == j,]
            all_gs_other<-all_gs[!(all_gs$Panel %in% all_gs_panel$Panel) & !(all_gs$Method %in% all_gs_panel$Method),]
            all_gs_other<-all_gs_other[!(all_gs_other$Name %in% all_gs_panel$Name),]
            all_gs_other<-unique(all_gs_other$Name)
            
            if(length(all_gs_other) > 0){
              all_gs_panel_missing<-data.frame(Name=all_gs_other)
              all_gs_panel_missing$Panel=i
              all_gs_panel_missing$Name=all_gs_other
              all_gs_panel_missing$Z=NA
              all_gs_panel_missing$P=NA
              all_gs_panel_missing$P.FDR=NA
              all_gs_panel_missing$Method=j
              
              all_gs_panel_missing<-all_gs_panel_missing[,names(all_gs_panel)]
              
              all_gs_all<-rbind(all_gs_all,all_gs_panel_missing)
            }
            
            all_gs_all<-rbind(all_gs_all,all_gs_panel)
          }
        }
        
        # Now remove the NA rows
        all_gs_all<-all_gs_all[!is.na(all_gs_all$Name),]
        methods<-c('MAGMA','GCSC','TWAS-GSEA')[c('MAGMA','GCSC','TWAS-GSEA') %in% all_gs_all$Method]
        all_gs_all$Method<-factor(all_gs_all$Method, levels=methods)
        
        # Sort according to user input
        if(input$selected_sort_drug == 'Alphabetical'){
          all_gs_all$Name<-factor(all_gs_all$Name, levels=unique(all_gs_all$Name[rev(order(all_gs_all$Name))]))
        }
        if(input$selected_sort_drug == 'All - Z'){
          all_gs_all$Name<-factor(all_gs_all$Name, levels=rev(unique(rev(all_gs_all$Name[order(all_gs_all$Z, na.last=F)]))))
        }
        if(input$selected_sort_drug == 'TWAS-GSEA - Z'){
          all_gs_all$Name <- factor(all_gs_all$Name, levels = rev(unique(rev(all_gs_all$Name[all_gs_all$Method == 'TWAS-GSEA'][order(all_gs_all$Z[all_gs_all$Method == 'TWAS-GSEA'], na.last=F)]))))
        }
        if(input$selected_sort_drug == 'MAGMA - Z'){
          all_gs_all$Name <- factor(all_gs_all$Name, levels = unique(all_gs_all$Name[all_gs_all$Method == 'MAGMA'][order(all_gs_all$Z[all_gs_all$Method == 'MAGMA'], na.last=F)]))
        }
        if(input$selected_sort_drug == 'GCSC - Z'){
          all_gs_all$Name <- factor(all_gs_all$Name, levels = unique(all_gs_all$Name[all_gs_all$Method == 'GCSC'][order(all_gs_all$Z[all_gs_all$Method == 'GCSC'], na.last=F)]))
        }
        
        group_siz<-NULL
        for(i in methods){
          group_siz<-rbind(group_siz, data.frame(Group=i,
                                                 Size=length(unique(all_gs_all$Panel[all_gs_all$Method==i]))))
        }
        
        # Set minimum size to 3 to allow space for labels
        group_siz$Size[group_siz$Size < 2]<-2
        group_siz$Prop<-group_siz$Size/sum(group_siz$Size)
        group_siz$Width<-4*group_siz$Prop
        
        x<-c(-max(abs(all_gs_all$Z), na.rm=T),0,max(abs(all_gs_all$Z), na.rm=T))
        x<-(x-min(x))/(max(x)-min(x))
        
        heatmap<-ggplot(data = all_gs_all, aes(x = Panel, y = Name)) +
          theme_bw()	+
          geom_point(data=all_gs_all, aes(x = Panel, y = Name, colour = Z), size=5) +
          geom_point(data=all_gs_all[which(all_gs_all$P < 0.05),], aes(x = Panel, y = Name), colour='black', fill=NA, size=6) +
          geom_point(data=all_gs_all[which(all_gs_all$P.FDR < 0.05),], aes(x = Panel, y = Name), colour='black', fill=NA, size=7, shape=15) +
          geom_point(data=all_gs_all, aes(x = Panel, y = Name, colour = Z), size=5) +
          # For reason, the factor-based sorted gets messed up with the Z point, but specifying twice fixes it?
          scale_colour_gradientn(colours=c("#0066FF","#0099FF","#FFFFFF","#FF6666","#FF0000"), na.value = NA,name = "Z-score", limits = c(-max(abs(all_gs_all$Z), na.rm=T), max(abs(all_gs_all$Z), na.rm=T)), values=x) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
          labs(x='', y='') +
          facet_wrap(~ Method , nrow=1, scales = "free_x")
        
        library(grid)
        gt = ggplot_gtable(ggplot_build(heatmap))
        
        for(i in 1:nrow(group_siz)){
          gt$widths[gt$layout$l[grep(paste0('panel-',i,'-1'), gt$layout$name)]] = group_siz$Width[i]*gt$widths[gt$layout$l[grep(paste0('panel-',i,'-1'), gt$layout$name)]]
        }
        
        grid.draw(gt)
        
      } else {
        NULL
      }
    })
    
    output$tx_drug_plot.ui <- renderUI({
      if(plot_dim_drug()[['height']] < 10000 & nrow(tx_drug_summary_data_filtered()) > 0){
        plotOutput("tx_drug_plot", height = plot_dim_drug()[['height']], width=plot_dim_drug()[['width']])
      } else {
        NULL
      }
    })
    
    output$message_too_large_drug <- renderUI({
      if(plot_dim_drug()[['height']] > 10000){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "Plot is too large. Restrict to significant drugs or specify a list of drugs."
        ))
      }
    })
    
    output$message_no_drugs_drug <- renderUI({
      if(nrow(tx_drug_summary_data_filtered()) == 0){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "No drugs are present."
        ))        
      }
    })
    
    #######
    # Prepare data for atc-specific association tables
    #######
    
    output$tx_atc_magma_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[2];",
        "  $('td:eq(2)', row).html(x.toExponential(2));",
        "  var y = data[3];",
        "  $('td:eq(3)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$tx$atc$magma
      tmp$Name<-paste0(tmp$`ATC Code`,': ',tmp$`ATC Description`)
      tmp<-tmp[,c('Name','N Drugs','P','P.FDR'), with=F]
      
      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = 0:3),
            list(width = '60px', targets = 2:3)
          )),
        escape = FALSE
      )
    })
    
    output$tx_atc_gcsc_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[2];",
        "  $('td:eq(2)', row).html(x.toExponential(2));",
        "  var y = data[3];",
        "  $('td:eq(3)', row).html(y.toExponential(2));",
        "}"
      )
      
      tmp<-gwas_data()[[gwas_selected]]$tx$atc$gcsc
      tmp$Name<-paste0(tmp$`ATC Code`,': ',tmp$`ATC Description`)
      tmp<-tmp[,c('Name','N Drugs','P','P.FDR'), with=F]
      
      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = 0:3),
            list(width = '60px', targets = 2:3)
          )),
        escape = FALSE
      )
    })
    
    output$tx_atc_twas_gsea_table<-renderDataTable({
      tmp<-gwas_data()[[gwas_selected]]$tx$atc$twas_gsea
      
      tmp$P.FDR_all<-p.adjust(tmp$P, method = 'fdr')
      tmp$P.FDR.onside_all<-p.adjust(tmp$P.oneside, method = 'fdr')
      
      tmp$Z<--qnorm(tmp$P)
      tmp$Z<-tmp$Z*sign(tmp$Estimate)
      tmp$Name<-paste0(tmp$`ATC Code`,': ',tmp$`ATC Description`)
      
      tmp$Estimate<-round(tmp$Estimate,3)
      
      tmp<-tmp[,c("Name","Panel","N Drugs", "Estimate","P","P.FDR_all"), with=F]
      names(tmp)<-c("Name","Panel","N Drugs","Estimate","P","P.FDR")
      
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[4];",
        "  $('td:eq(4)', row).html(x.toExponential(2));",
        "  var y = data[5];",
        "  $('td:eq(5)', row).html(y.toExponential(2));",
        "}"
      )
      
      datatable(
        tmp,
        rownames = F,
        options = list(# Apply javascript for P value column
          rowCallback = JS(js),
          # Centre column contents and fix width of Pvalue column
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '60px', targets = 4:5)
          )),
        escape = FALSE
      )
    })
    
    #######
    # Prepare data for atc enrichment plot
    #######
    
    tx_atc_summary_data <- reactive({
      
      ###
      # MAGMA
      ###
      
      magma_gs_atc<-gwas_data()[[gwas_selected]]$tx$atc$magma
      
      magma_gs_atc$Z<--qnorm(magma_gs_atc$P)
      magma_gs_atc$FDR_Sig<-magma_gs_atc$P.FDR < 0.05
      magma_gs_atc$Nom_Sig<-magma_gs_atc$P < 0.05
      magma_gs_atc$Name<-paste0(magma_gs_atc$`ATC Code`,': ',magma_gs_atc$`ATC Description`)
      magma_gs_atc$Method<-'MAGMA'
      magma_gs_atc$Panel<-'MAGMA'
      
      magma_gs_atc<-magma_gs_atc[,c("Name","Z","FDR_Sig","Nom_Sig","Method","Panel"), with=F]
      
      ###
      # GCSC
      ###
      
      gcsc_gs_atc<-gwas_data()[[gwas_selected]]$tx$atc$gcsc
      
      gcsc_gs_atc$Z<--qnorm(gcsc_gs_atc$P)
      gcsc_gs_atc$FDR_Sig<-gcsc_gs_atc$P.FDR < 0.05
      gcsc_gs_atc$Nom_Sig<-gcsc_gs_atc$P < 0.05
      gcsc_gs_atc$Name<-paste0(gcsc_gs_atc$`ATC Code`,': ',gcsc_gs_atc$`ATC Description`)
      gcsc_gs_atc$Method<-'GCSC'
      gcsc_gs_atc$Panel<-'GCSC'
      
      gcsc_gs_atc<-gcsc_gs_atc[,c("Name","Z","FDR_Sig","Nom_Sig","Method","Panel"), with=F]
      
      ###
      # TWAS-GSEA
      ###
      
      gsea_gs_atc<-gwas_data()[[gwas_selected]]$tx$atc$twas_gsea
      
      gsea_gs_atc$P.FDR_all<-p.adjust(gsea_gs_atc$P, method = 'fdr')
      gsea_gs_atc$P.FDR.onside_all<-p.adjust(gsea_gs_atc$P.oneside, method = 'fdr')
      
      gsea_gs_atc$Z<--qnorm(gsea_gs_atc$P)
      gsea_gs_atc$Z<-gsea_gs_atc$Z*sign(gsea_gs_atc$Estimate)
      gsea_gs_atc$FDR_Sig<-gsea_gs_atc$P.FDR_all < 0.05
      gsea_gs_atc$Nom_Sig<-gsea_gs_atc$P < 0.05
      gsea_gs_atc$Name<-paste0(gsea_gs_atc$`ATC Code`,': ',gsea_gs_atc$`ATC Description`)
      
      gsea_gs_atc$Method<-'TWAS-GSEA'
      
      gsea_gs_atc<-gsea_gs_atc[,c("Name","Z","FDR_Sig","Nom_Sig","Method","Panel"), with=F]
      
      # Insert missing values
      gsea_gs_atc_all<-gsea_gs_atc
      for(i in unique(gsea_gs_atc_all$Panel)){
        gsea_gs_atc_i<-gsea_gs_atc[gsea_gs_atc$Panel == i,]
        gsea_gs_atc_other<-gsea_gs_atc[gsea_gs_atc$Panel != i,]
        gsea_gs_atc_rest<-data.frame(Name=unique(gsea_gs_atc_other$Name[!(gsea_gs_atc_other$Name %in% gsea_gs_atc_i$Name)]),
                                     Z =NA,
                                     FDR_Sig=NA,
                                     Nom_Sig=NA,
                                     Method='TWAS-GSEA',
                                     Panel=i)
        
        gsea_gs_atc_all<-rbind(gsea_gs_atc_all, gsea_gs_atc_rest)
        
      }
      gsea_gs_atc<-gsea_gs_atc_all
      
      all_gs_atc<-do.call(rbind, list(magma_gs_atc, gcsc_gs_atc, gsea_gs_atc))
      
      return(all_gs_atc)
    })
    
    observe({
      all_gs<-tx_atc_summary_data()
      methods<-unique(all_gs$Method)
      updateSelectInput(session, "selected_methods_atc", choices = methods, selected=methods)
    })
    
    observe({
      all_gs<-tx_atc_summary_data()
      expr_panels<-unique(all_gs$Panel[all_gs$Method == 'TWAS-GSEA'])
      updateSelectInput(session, "selected_expr_panels_atc", choices = expr_panels, selected=expr_panels)
    })
    
    tx_atc_summary_data_filtered<-reactive({
      all_gs_atc<-tx_atc_summary_data()
      
      # Filter results table by user specified methods
      all_gs_atc<-all_gs_atc[all_gs_atc$Method %in% input$selected_methods_atc,]
      
      # Filter results table by user specified expression
      if(any(all_gs_atc$Method == 'TWAS-GSEA')){
        all_gs_atc<-all_gs_atc[!(all_gs_atc$Method == 'TWAS-GSEA' & !(all_gs_atc$Panel %in% input$selected_expr_panels_atc)),]
      }
      
      # Filter results table if user specifies high confidence genes only
      if(input$conf_only_atc){
        all_gs_atc<-all_gs_atc[all_gs_atc$Name %in% all_gs_atc$Name[which(all_gs_atc$FDR_Sig == T)],]
      }
      
      input_atcs <- unlist(strsplit(input$atcInput_atc, "[, ]"))
      selected_atcs <- input_atcs[input_atcs != ""]
      
      # Insert NA rows for all panels and methods so when filtering by atc, all selected panels and methods remain
      na_rows<-all_gs_atc[!(duplicated(paste0(all_gs_atc$Panel, all_gs_atc$Method))),]
      na_rows$Name<-NA
      na_rows$Z<-NA
      na_rows$FDR_Sig<-NA
      na_rows$Nom_Sig<-NA
      
      all_gs_atc<-rbind(na_rows, all_gs_atc)
      
      if(length(selected_atcs) > 0){
        if(sum(grepl(paste(selected_atcs, collapse='|'), all_gs_atc$Name, ignore.case = T)) > 0){
          all_gs_atc<-all_gs_atc[grepl(paste(selected_atcs, collapse='|'), all_gs_atc$Name, ignore.case = T) & !is.na(all_gs_atc$Name),]
        } else {
          all_gs_atc<-data.frame(matrix(nrow=0, ncol=5))
        }
      }
      
      return(all_gs_atc)
    })
    
    # Identify number of atcs
    plot_dim_atc<-reactive({
      all_gs_atc<-tx_atc_summary_data_filtered()
      
      if(nrow(all_gs_atc) > 0){
        num_row <- length(unique(all_gs_atc$Name))
        plot_height<-(max(nchar(all_gs_atc$Panel))*2)+(num_row * 20)+100
        num_col <- length(unique(paste0(all_gs_atc$Panel,'_',all_gs_atc$Method,'_')))
        plot_width<-150+(max(nchar(all_gs_atc$Name), na.rm=T)*2)+(num_col * 50)
        plot_width<-max(plot_width,(length(unique(all_gs_atc$Method))*100))
      } else {
        plot_height<-100
        plot_width<-100
      }
      
      return(list(height=plot_height,
                  width=plot_width))
    })
    
    observe({
      tmp<-tx_atc_summary_data_filtered()
      choices<-'All - Z'
      if(any(tmp$Method == 'MAGMA')){
        choices<-c(choices, 'MAGMA - Z')
      }
      if(any(tmp$Method == 'GCSC')){
        choices<-c(choices, 'GCSC - Z')
      }
      if(any(tmp$Method == 'TWAS.GSEA')){
        choices<-c(choices, 'TWAS.GSEA - Z')
      }
      if(length(unique(tmp$Method)) == 1){
        choices<-choices[choices != 'All - Z']
      }
      choices<-c(choices, 'Alphabetical')
      
      updateSelectInput(session, "selected_sort_atc", choices = choices, selected=choices[1])
    })
    
    output$tx_atc_plot<-renderPlot({
      
      all_gs_atc<-tx_atc_summary_data_filtered()
      
      if(plot_dim_atc()[['height']] < 10000 & nrow(all_gs_atc) > 0){
        
        # Insert missing data
        all_gs_atc_all<-NULL
        for(i in unique(all_gs_atc$Panel)){
          for(j in unique(all_gs_atc$Method[all_gs_atc$Panel == i])){
            
            all_gs_atc_panel<-all_gs_atc[all_gs_atc$Panel == i & all_gs_atc$Method == j,]
            all_gs_atc_other<-all_gs_atc[!(all_gs_atc$Panel %in% all_gs_atc_panel$Panel) & !(all_gs_atc$Method %in% all_gs_atc_panel$Method),]
            all_gs_atc_other<-all_gs_atc_other[!(all_gs_atc_other$Name %in% all_gs_atc_panel$Name),]
            all_gs_atc_other<-unique(all_gs_atc_other$Name)
            
            if(length(all_gs_atc_other) > 0){
              all_gs_atc_panel_missing<-data.frame(Name=all_gs_atc_other)
              all_gs_atc_panel_missing$Panel=i
              all_gs_atc_panel_missing$Name=all_gs_atc_other
              all_gs_atc_panel_missing$Z=NA
              all_gs_atc_panel_missing$FDR_Sig=NA
              all_gs_atc_panel_missing$Nom_Sig=NA
              all_gs_atc_panel_missing$Method=j
              
              
              all_gs_atc_panel_missing<-all_gs_atc_panel_missing[,names(all_gs_atc_panel)]
              
              all_gs_atc_all<-rbind(all_gs_atc_all,all_gs_atc_panel_missing)
            }
            
            all_gs_atc_all<-rbind(all_gs_atc_all,all_gs_atc_panel)
          }
        }
        
        # Now remove the NA rows
        all_gs_atc_all<-all_gs_atc_all[!is.na(all_gs_atc_all$Name),]
        methods<-c('MAGMA','GCSC','TWAS-GSEA')[c('MAGMA','GCSC','TWAS-GSEA') %in% all_gs_atc_all$Method]
        all_gs_atc_all$Method<-factor(all_gs_atc_all$Method, levels=methods)
        
        # Shorten long ATC descriptions
        for(i in unique(all_gs_atc_all$Name)){
          if(nchar(i) > 30){
            i_new<-paste0(substr(i, 1, 27),'...')
            i_new<-gsub(' \\.\\.\\.','...',i_new)
            all_gs_atc_all$Name[all_gs_atc_all$Name == i]<-i_new
          }
        }
        
        # Sort according to user input
        if(input$selected_sort_atc == 'Alphabetical'){
          all_gs_atc_all$Name<-factor(all_gs_atc_all$Name, levels=unique(all_gs_atc_all$Name[rev(order(all_gs_atc_all$Name))]))
        }
        if(input$selected_sort_atc == 'All - Z'){
          all_gs_atc_all$Name<-factor(all_gs_atc_all$Name, levels=rev(unique(rev(all_gs_atc_all$Name[order(all_gs_atc_all$Z, na.last=F)]))))
        }
        if(input$selected_sort_atc == 'TWAS-GSEA - Z'){
          all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = rev(unique(rev(all_gs_atc_all$Name[all_gs_atc_all$Method == 'TWAS-GSEA'][order(all_gs_atc_all$Z[all_gs_atc_all$Method == 'TWAS-GSEA'], na.last=F)]))))
        }
        if(input$selected_sort_atc == 'MAGMA - Z'){
          all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = unique(all_gs_atc_all$Name[all_gs_atc_all$Method == 'MAGMA'][order(all_gs_atc_all$Z[all_gs_atc_all$Method == 'MAGMA'], na.last=F)]))
        }
        if(input$selected_sort_atc == 'GCSC - Z'){
          all_gs_atc_all$Name <- factor(all_gs_atc_all$Name, levels = unique(all_gs_atc_all$Name[all_gs_atc_all$Method == 'GCSC'][order(all_gs_atc_all$Z[all_gs_atc_all$Method == 'GCSC'], na.last=F)]))
        }
        
        group_siz<-NULL
        for(i in methods){
          group_siz<-rbind(group_siz, data.frame(Group=i,
                                                 Size=length(unique(all_gs_atc_all$Panel[all_gs_atc_all$Method==i]))))
        }
        
        # Set minimum size to 3 to allow space for labels
        group_siz$Size[group_siz$Size < 2]<-2
        group_siz$Prop<-group_siz$Size/sum(group_siz$Size)
        group_siz$Width<-4*group_siz$Prop
        
        x<-c(-max(abs(all_gs_atc_all$Z), na.rm=T),0,max(abs(all_gs_atc_all$Z), na.rm=T))
        x<-(x-min(x))/(max(x)-min(x))
        
        heatmap<-ggplot(data = all_gs_atc_all, aes(x = Panel, y = Name)) +
          theme_bw() +
          geom_point(data=all_gs_atc_all, aes(colour = Z), size=5) +
          geom_point(data=all_gs_atc_all[which(all_gs_atc_all$Nom_Sig == T),], aes(x = Panel, y = Name), colour='black', fill=NA, size=6) +
          geom_point(data=all_gs_atc_all[which(all_gs_atc_all$FDR_Sig == T),], aes(x = Panel, y = Name), colour='black', fill=NA, size=7, shape=15) +
          geom_point(data=all_gs_atc_all, aes(colour = Z), size=5) +
          scale_colour_gradientn(colours=c("#0066FF","#0099FF","#FFFFFF","#FF6666","#FF0000"), na.value = NA,name = "Z-score", limits = c(min(all_gs_atc_all$Z, na.rm=T), max(all_gs_atc_all$Z, na.rm=T)), values=x) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
          labs(x='', y='') +
          facet_wrap(~ Method , nrow=1, scales = "free_x")
        
        
        library(grid)
        gt = ggplot_gtable(ggplot_build(heatmap))
        
        for(i in 1:nrow(group_siz)){
          gt$widths[gt$layout$l[grep(paste0('panel-',i,'-1'), gt$layout$name)]] = group_siz$Width[i]*gt$widths[gt$layout$l[grep(paste0('panel-',i,'-1'), gt$layout$name)]]
        }
        
        grid.draw(gt)
        
      } else {
        NULL
      }
    })
    
    output$tx_atc_plot.ui <- renderUI({
      if(plot_dim_atc()[['height']] < 10000 & nrow(tx_atc_summary_data_filtered()) > 0){
        plotOutput("tx_atc_plot", height = plot_dim_atc()[['height']], width=plot_dim_atc()[['width']])
      } else {
        NULL
      }
    })
    
    output$message_too_large_atc <- renderUI({
      if(plot_dim_atc()[['height']] > 10000){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "Plot is too large. Restrict to significant atcs or specify a list of atcs."
        ))
      }
    })
    
    output$message_no_atcs_atc <- renderUI({
      if(nrow(tx_atc_summary_data_filtered()) == 0){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "No ATC codes remain"
        ))        
      }
    })
  })
}

# Run the Shiny app
shinyApp(ui, server)






