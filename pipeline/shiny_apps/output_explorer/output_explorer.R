# Load required libraries
library(shiny)
library(data.table)
library(DT)
library(grid)
library(ggplot2)
library(cowplot)
library(shinythemes)

# Define UI for the Shiny app
ui <- fluidPage(

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
           h2("Welcome to GenoDiscover"),
           p("GenoDiscover is your comprehensive platform for Genome-Wide Association Study (GWAS) summary statistics analysis. Explore genetic associations, visualize results, and gain insights into your data with our user-friendly tools."),
           p("Get started by uploading your GWAS summary statistics and let GenoDiscover help you uncover meaningful patterns in your data."),
           hr(),
           
           # Instructions for submitting data
           h3("Submitting GWAS Summary Statistics"),
           p("You can submit your GWAS summary statistics for analysis via our King's College London server by going to the 'Submit' tab."),
           hr(),
           
           # Instructions for exploring existing results
           h3("Exploring Previous Results"),
           p("If you have already run the GenoDiscover pipeline and have results, you can explore them by uploading your results_package.rds file to the 'Explore' tab."),
           hr(),
           
           # Information on citing the platform
           h3("Citing GenoDiscover"),
           p("If you use GenoDiscover for a publication, please cite our publication describing the platform, as well as the underlying datasets and methods it uses."),
           p("Publication reference: TBD"),
           hr(),        
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
      title='Submit'
      
    ),
    tabPanel(
      title='Explore',
  
      tabsetPanel(
        tabPanel(
          title="Data Input",
          br(),
          p("This is an application for visualising the output of GenoDiscover. To start, upload the 'results_package.rds' file output by the GenoDiscover pipeline, select a GWAS, and use the tabs to view interactive tables and plots of your results."), 
          p("Click ",a("here", href = "https://github.com/opain/GenoDiscover"), " here to learn more about the pipeline. Please cite  ",a("our publication", href = "https://github.com/opain/GenoDiscover"), "  and relevent software and datasets included in your analysis."), 
          hr(),
          br(),
          fileInput("file", "Choose an .RDS file"),
          selectInput("gwas_selector", "Select a GWAS", ""),
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
          hr(),
          br(),
          fluidRow(
            column(width=6,
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
                  selectInput("selected_methods", "Select methods", "", multiple=T),
                  selectInput("selected_expr_panels", "Select expression panels", "", multiple=T),
                  selectInput("selected_protein_panels", "Select protein panels", "", multiple=T),
                  radioButtons("conf_only", "Show high-confidence only :",
                               choices = c("True" = T,
                                           "False" = F),
                               selected = T),
                  textInput("geneInput", "Enter gene symbols (whitespace- or comma-seperated):")
                ),
                
                mainPanel(
                  uiOutput("messageDiv1"),
                  uiOutput("messageDiv2"),
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
                column(width=6,
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
                column(width=6,
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
                column(width=6,
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
                column(width=6,
                       dataTableOutput("mol_assoc_smr_protein_table"),
                )
              ),
              br()
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
  
  gwas_data <- reactive({
    req(input$file)
    rds<-readRDS(input$file$datapath)
    rds
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

      datatable(gwas_data()[[gwas_selected]]$mol_assoc$magma, rownames=F, options = list(
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
      
      mol_assoc_fusion_expr_data<-gwas_data()[[gwas_selected]]$mol_assoc$exp$fusion$res
      mol_assoc_fusion_expr_data$TWAS.Z<-round(mol_assoc_fusion_expr_data$TWAS.Z, 3)
        
      datatable(mol_assoc_fusion_expr_data, rownames=F, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:11))))
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
      
      mol_assoc_fusion_protein_data<-gwas_data()[[gwas_selected]]$mol_assoc$protein$fusion$res
      mol_assoc_fusion_protein_data$pwas_all.Z<-round(mol_assoc_fusion_protein_data$pwas_all.Z, 3)
      
      datatable(mol_assoc_fusion_protein_data, rownames=F, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:11))))
    })
    
    output$mol_assoc_smr_expr_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[7];",
        "  $('td:eq(7)', row).html(x.toExponential(2));",
        "  var y = data[8];",
        "  $('td:eq(8)', row).html(y.toExponential(2));",
        "}"
      )
      
      datatable(gwas_data()[[gwas_selected]]$mol_assoc$exp$smr$res, rownames=F, options = list(
        # Apply javascript for P value column
        rowCallback = JS(js), 
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = 0:9))))
    })
    
    output$mol_assoc_smr_protein_table<-renderDataTable({
      # Create java script to force scientific notation for P and P.FDR value column, whilst allowing numeric sorting
      js <- c(
        "function(row, data, displayNum, index){",
        "  var x = data[6];",
        "  $('td:eq(6)', row).html(x.toExponential(2));",
        "  var y = data[7];",
        "  $('td:eq(7)', row).html(y.toExponential(2));",
        "}"
      )
      
      datatable(gwas_data()[[gwas_selected]]$mol_assoc$protein$smr$res, rownames=F, options = list(
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
      updateSelectInput(session, "selected_methods", choices = methods, selected=methods)
    })
    
    observe({
      all_func_res<-mol_assoc_summary_data()
      expr_panels<-unique(all_func_res$Panel[all_func_res$Type == 'Expr.'])
      updateSelectInput(session, "selected_expr_panels", choices = expr_panels, selected=expr_panels)
    })
    
    observe({
      all_func_res<-mol_assoc_summary_data()
      protein_panels<-unique(all_func_res$Panel[all_func_res$Type == 'Protein'])
      updateSelectInput(session, "selected_protein_panels", choices = protein_panels, selected=protein_panels)
    })
    
    mol_assoc_summary_data_filtered<-reactive({
      all_func_res<-mol_assoc_summary_data()
      
      # Filter results table by user specified methods
      all_func_res<-all_func_res[all_func_res$Method %in% input$selected_methods,]
      
      # Filter results table by user specified expression and protein panels
      if(any(all_func_res$Type == 'Expr.')){
        all_func_res<-all_func_res[!(all_func_res$Type == 'Expr.' & !(all_func_res$Panel %in% input$selected_expr_panels)),]
      }
      if(any(all_func_res$Type == 'Protein')){
        all_func_res<-all_func_res[!(all_func_res$Type == 'Protein' & !(all_func_res$Panel %in% input$selected_protein_panels)),]
      }
      
      # Filter results table if user specifies high confidence genes only
      if(input$conf_only){
        all_func_res<-all_func_res[all_func_res$ID %in% all_func_res$ID[which((all_func_res$Sig == T & all_func_res$Coloc == T) | all_func_res$Panel == "SuSie (L=1)")],]
      }
      
      input_genes <- unlist(strsplit(input$geneInput, "[, ]"))
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
    plot_dim<-reactive({
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
      
      if(plot_dim()[['height']] < 10000 & nrow(all_func_res) > 0){
        
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
    }, type = "cairo")
    
    output$mol_assoc_plot.ui <- renderUI({
      if(plot_dim()[['height']] < 10000 & nrow(mol_assoc_summary_data_filtered()) > 0){
        plotOutput("mol_assoc_plot", height = plot_dim()[['height']], width=plot_dim()[['width']])
      } else {
        NULL
      }
    })
    
    output$messageDiv1 <- renderUI({
      if(plot_dim()[['height']] > 10000){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "Plot is too large. Restrict to high-confidence genes or specify a list of genes."
        ))
      }
    })
    
    output$messageDiv2 <- renderUI({
      if(nrow(mol_assoc_summary_data_filtered()) == 0){
        HTML(sprintf(
          "<div style='color: red;'>%s</div>",
          "No genes are present."
        ))        
      }
    })
  })
}

# Run the Shiny app
shinyApp(ui, server)






