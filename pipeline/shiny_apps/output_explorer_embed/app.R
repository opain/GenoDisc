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
source('../../scripts/functions/sumstat_cleaner_functions.R')
source('../output_explorer/functions.R')

options(shiny.maxRequestSize = 600 * 1024 * 1024)

# Define UI for the Shiny app
ui <- fluidPage(

  shinyjs::useShinyjs(),

  tags$style(HTML("
    .custom-panel {
          padding: 20px;
          border: 1px solid #ccc;
    }
  ")),

  theme = shinythemes::shinytheme("paper"),

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
                  textInput("geneInput_mol", "Enter gene symbols (whitespace- or comma-seperated):"),
                  hr(),
                  h5('Select high confidence criteria:'),
                  selectInput("selected_group_hc_mol", "Select methods", "", multiple=T)
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
            ),
            tabPanel(
              title="Panel Info.",
              br(),
              p("This tab shows the number of features and individuals for each panel."),
              hr(),
              br(),
              fluidRow(
                column(width=7,
                       dataTableOutput("panel_info_table"),
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
        ),
        tabPanel(
          title="References",
          br(),
          p("Please be sure to cite the software and datasets used by this pipeline. Relevent citations are shown below:"),
          fluidRow(
            column(width=8,
                   dataTableOutput("reference_table"),
            )
          ),
        ),
        tabPanel(
          title="Configuration",

          br(),
          div(
            class = "custom-panel",
            h6(strong('Repository information:')),
            uiOutput("repo_info")
          ),
          hr(),
          div(
            class = "custom-panel",
            h6(strong('Config file:')),
            DT::dataTableOutput("config_table"),
          ),
          hr(),
          div(
            class = "custom-panel",
            h6(strong('GWAS list:')),
            DT::dataTableOutput("gwas_list")
          ),
          hr()
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

    ######
    # Prepare table for GWAS QC tab
    ######

    # Create a table showing key statistics
    qc_val<-data.table(name=selected_gwas(),
                       label=gwas_list$label[gwas_list$name == selected_gwas()],
                       n_var_orig=gwas_data()[[selected_gwas()]]$gwas_qc$cleaner_dat$val$n_var_orig,
                       build=gwas_data()[[selected_gwas()]]$gwas_qc$cleaner_dat$val$build$build,
                       n_snp_final=gwas_data()[[selected_gwas()]]$gwas_qc$cleaner_dat$val$n_snp_final,
                       lambda_gc=gwas_data()[[selected_gwas()]]$gwas_qc$focus_dat$val$lambda_gc,
                       max_chi2=gwas_data()[[selected_gwas()]]$gwas_qc$focus_dat$val$max_chi2,
                       n_sig_snp=gwas_data()[[selected_gwas()]]$gwas_qc$focus_dat$val$n_sig_snp,
                       obs_h2=paste0(round(gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_dat$val$obs_h2_est,3), " (",round(gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_dat$val$obs_h2_se,3),")"),
                       int=paste0(round(gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_dat$val$int_est,3), " (",round(gwas_data()[[selected_gwas()]]$gwas_qc$ldsc_dat$val$int_se,3),")"))

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

    snp_assoc_data <- gwas_data()[[selected_gwas()]]$snp_assoc$clump

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
        snp_assoc_lead <- gwas_data()[[selected_gwas()]]$snp_assoc$clump
      }

      if (input$clumping_type == "cojo_analysis") {
        # Read in the COJO associations
        snp_assoc_lead <- gwas_data()[[selected_gwas()]]$snp_assoc$cojo
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
        snp_assoc_finemap <- gwas_data()[[selected_gwas()]]$snp_assoc$susie$L1
      }

      if (input$l_param == "L10") {
        snp_assoc_finemap <- gwas_data()[[selected_gwas()]]$snp_assoc$susie$L10
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

      tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$magma

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

      tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$exp$fusion$res
      tmp$TWAS.Z<-round(tmp$TWAS.Z, 3)
      tmp$`High Confidence`  <- tmp$TWAS.P.FDR < 0.05 & tmp$COLOC_logical
      tmp$COLOC_logical<-NULL
      tmp[is.na(tmp)]<-'NA'

      names(tmp)[names(tmp) == 'TWAS.Z']<-'Z'
      names(tmp)[names(tmp) == 'TWAS.P']<-'P'
      names(tmp)[names(tmp) == 'TWAS.P.FDR']<-'P.FDR'

      tmp<-tmp[, c("PANEL","CHR","P0","P1","Ensembl ID","Gene Symbol","Z","P","P.FDR","COLOC.PP3","COLOC.PP4","High Confidence"), with=F]

      datatable(tmp, rownames=F, options = list(
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

      tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$protein$fusion$res
      tmp$pwas_all.Z<-round(tmp$pwas_all.Z, 3)
      tmp$`High Confidence`  <- tmp$pwas_all.P.FDR < 0.05 & tmp$COLOC_logical
      tmp$COLOC_logical<-NULL
      tmp[is.na(tmp)]<-'NA'

      names(tmp)[names(tmp) == 'pwas_all.Z']<-'Z'
      names(tmp)[names(tmp) == 'pwas_all.P']<-'P'
      names(tmp)[names(tmp) == 'pwas_all.P.FDR']<-'P.FDR'

      tmp<-tmp[, c("PANEL","CHR","P0","P1","Ensembl ID","Gene Symbol","Z","P","P.FDR","COLOC.PP3","COLOC.PP4","High Confidence"), with=F]

      datatable(tmp, rownames=F, options = list(
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
        "  var x = data[9];",
        "  $('td:eq(9)', row).html(x.toExponential(2));",
        "}"
      )

      tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$res
      tmp$`High Confidence`<-tmp$p_SMR.FDR < 0.05 & tmp$p_HEIDI > 0.05
      tmp<-tmp[,c("PANEL","CHR","BP","Ensembl ID","Gene Symbol","b_SMR","se_SMR","p_SMR","p_SMR.FDR","p_HEIDI","High Confidence"), with=F]
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
        columnDefs = list(list(className = 'dt-center', targets = 0:10))))
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

      tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$protein$smr$res
      tmp$`High Confidence`<-tmp$p_SMR.FDR < 0.05 & tmp$p_HEIDI > 0.05
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
    # Create table showing panel information
    #######


    output$panel_info_table<-renderDataTable({

      all_panel_info<-NULL

      if(twas_logical){
        tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$exp$fusion$panels
        all_panel_info<-rbind(all_panel_info, tmp)
      }

      if(smr_expression_logical){
        tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$panels
        all_panel_info<-rbind(all_panel_info, tmp)
      }

      if(any(pwas_panel_rosmap_logical, pwas_panel_banner_logical)){
        tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$protein$fusion$panels
        all_panel_info<-rbind(all_panel_info, tmp)
      }

      if(smr_protein_panel_rosmap_logical){
        tmp<-gwas_data()[[selected_gwas()]]$mol_assoc$protein$smr$panels
        all_panel_info<-rbind(all_panel_info, tmp)
      }

      all_panel_info<-all_panel_info[order(all_panel_info$Software, all_panel_info$Type, all_panel_info$Panel),]

      names(all_panel_info)[names(all_panel_info) == 'N_indiv']<-'N Individuals'
      names(all_panel_info)[names(all_panel_info) == 'N_gene']<-'N Genes'

      datatable(all_panel_info, rownames=F, options = list(
        # Centre column contents and fix width of Pvalue column
        columnDefs = list(list(className = 'dt-center', targets = '_all'))))
    })

    #######
    # Prepare data for molecular association plot
    #######

    mol_assoc_summary_data <- reactive({

      all_func_res<-NULL

      if(finemap_logical){

        finemap_L1_tmp<-data.frame(Panel = "SuSie (L=1)",
                                   ID=gwas_data()[[selected_gwas()]]$mol_assoc$finemap$L1,
                                   Z=1,
                                   Sig=F,
                                   Coloc=F,
                                   Method="SNP\nFine-mapping",
                                   Type='')

        all_func_res<-rbind(all_func_res, finemap_L1_tmp)

      }

      if(twas_logical){

        twas_tmp<-data.table(Panel=gwas_data()[[selected_gwas()]]$mol_assoc$exp$fusion$res$PANEL,
                             ID=gwas_data()[[selected_gwas()]]$mol_assoc$exp$fusion$res$`Gene Symbol`,
                             Z=gwas_data()[[selected_gwas()]]$mol_assoc$exp$fusion$res$TWAS.Z,
                             Sig=gwas_data()[[selected_gwas()]]$mol_assoc$exp$fusion$res$TWAS.P.FDR < 0.05,
                             Coloc=gwas_data()[[selected_gwas()]]$mol_assoc$exp$fusion$res$COLOC_logical)

        twas_tmp$Method<-'FUSION'
        twas_tmp$Type<-'Expr.'
        twas_tmp$Type[grepl('SPLIC',twas_tmp$Panel, ignore.case = T)]<-'Splice'

        # Retain only the most significant assoc for each gene within PANEL (only relevent for splice panel)
        twas_tmp<-twas_tmp[order(-abs(twas_tmp$Z)),]
        twas_tmp<-twas_tmp[!duplicated(paste0(twas_tmp$Panel, twas_tmp$ID)),]

        twas_tmp$Type<-factor(twas_tmp$Type, levels=c('Expr.','Splice'))
        twas_tmp<-twas_tmp[order(twas_tmp$Type),]

        all_func_res<-rbind(all_func_res, twas_tmp)
      }

      if(smr_expression_logical){
        # SMR expression
        smr_expr_id<-gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$res$`Gene Symbol`
        smr_expr_id[is.na(smr_expr_id)]<-gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$res$`Ensembl ID`[is.na(smr_expr_id)]
        smr_expression_res_tmp<-data.table(Panel=gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$res$PANEL,
                                           ID=smr_expr_id,
                                           Z=gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$res$b_SMR/gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$res$se_SMR,
                                           Sig=gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$res$p_SMR.FDR < 0.05,
                                           Coloc=gwas_data()[[selected_gwas()]]$mol_assoc$exp$smr$res$p_HEIDI > 0.05)

        smr_expression_res_tmp$Method<-'SMR'
        smr_expression_res_tmp$Type<-'Expr.'

        all_func_res<-rbind(all_func_res, smr_expression_res_tmp)
      }

      # PWAS
      if(any(pwas_panel_rosmap_logical, pwas_panel_banner_logical)){

        pwas_tmp<-data.table( Panel=gwas_data()[[selected_gwas()]]$mol_assoc$protein$fusion$res$PANEL,
                              ID=gwas_data()[[selected_gwas()]]$mol_assoc$protein$fusion$res$`Gene Symbol`,
                              Z=gwas_data()[[selected_gwas()]]$mol_assoc$protein$fusion$res$pwas_all.Z,
                              Sig=gwas_data()[[selected_gwas()]]$mol_assoc$protein$fusion$res$pwas_all.P.FDR < 0.05,
                              Coloc=gwas_data()[[selected_gwas()]]$mol_assoc$protein$fusion$res$COLOC_logical)

        pwas_tmp<-pwas_tmp[order(-abs(pwas_tmp$Z)),]
        pwas_tmp<-pwas_tmp[!duplicated(paste0(pwas_tmp$Panel, pwas_tmp$ID)),]
        pwas_tmp$Method<-'FUSION'
        pwas_tmp$Type<-'Protein'

        all_func_res<-rbind(all_func_res, pwas_tmp)
      }

      if(smr_protein_panel_rosmap_logical){

        # SMR protein
        smr_protein_res_tmp<-data.table(Panel=gwas_data()[[selected_gwas()]]$mol_assoc$protein$smr$res$PANEL,
                                        ID=gwas_data()[[selected_gwas()]]$mol_assoc$protein$smr$res$`Gene Symbol`,
                                        Z=gwas_data()[[selected_gwas()]]$mol_assoc$protein$smr$res$b_SMR/gwas_data()[[selected_gwas()]]$mol_assoc$protein$smr$res$se_SMR,
                                        Sig=gwas_data()[[selected_gwas()]]$mol_assoc$protein$smr$res$p_SMR.FDR < 0.05,
                                        Coloc=gwas_data()[[selected_gwas()]]$mol_assoc$protein$smr$res$p_HEIDI > 0.05)

        smr_protein_res_tmp<-smr_protein_res_tmp[order(-abs(smr_protein_res_tmp$Z)),]
        smr_protein_res_tmp<-smr_protein_res_tmp[!duplicated(paste0(smr_protein_res_tmp$Panel, smr_protein_res_tmp$ID)),]
        smr_protein_res_tmp$Method<-'SMR'
        smr_protein_res_tmp$Type<-'Protein'

        all_func_res<-rbind(all_func_res, smr_protein_res_tmp)

      }

      if(magma_gene_logical){

        magma_tmp<-data.frame(Panel = 'MAGMA',
                              ID=gwas_data()[[selected_gwas()]]$mol_assoc$magma$ID,
                              Z=abs(qnorm(as.numeric(gwas_data()[[selected_gwas()]]$mol_assoc$magma$P))),
                              Sig=as.numeric(gwas_data()[[selected_gwas()]]$mol_assoc$magma$P.FDR) < 0.05,
                              Coloc=F,
                              Method='MAGMA',
                              Type='')

        all_func_res<-rbind(all_func_res, magma_tmp)

      }

      if(clump_logical){

        nearest_tmp<-data.frame(Panel = 'NearestGene',
                                ID=gwas_data()[[selected_gwas()]]$mol_assoc$nearest$clump,
                                Z=1,
                                Sig=F,
                                Coloc=F,
                                Method='Nearest\nGene',
                                Type='')

        all_func_res<-rbind(all_func_res, nearest_tmp)

      }

      return(all_func_res)
    })

    observeEvent(mol_assoc_summary_data(), {
      all_func_res<-mol_assoc_summary_data()

      methods<-unique(all_func_res$Method)
      updateSelectInput(session, "selected_methods_mol", choices = methods, selected=methods)

      expr_panels<-unique(all_func_res$Panel[all_func_res$Type == 'Expr.' | all_func_res$Type == 'Splice'])
      updateSelectInput(session, "selected_expr_panels_mol", choices = expr_panels, selected=expr_panels)

      protein_panels<-unique(all_func_res$Panel[all_func_res$Type == 'Protein'])
      updateSelectInput(session, "selected_protein_panels_mol", choices = protein_panels, selected=protein_panels)

      # Select groups used to define high confidence genes
      res_group<-paste0(all_func_res$Method,'\n',all_func_res$Type )
      res_group[res_group == 'SNP\nFine-mapping\n']<-'SuSiE'
      res_group[res_group == 'MAGMA\n']<-'MAGMA'
      res_group[res_group == 'Nearest\nGene\n']<-'Nearest\nGene'

      hc_groups<-c('SuSiE','FUSION\nExpr.','FUSION\nSplice','SMR\nExpr.','FUSION\nProtein','SMR\nProtein')
      hc_groups<-hc_groups[hc_groups %in% res_group]

      updateSelectInput(session, "selected_group_hc_mol", choices = hc_groups[hc_groups %in% res_group], selected=hc_groups)

    })

    mol_assoc_summary_data_filtered<-reactive({
      all_func_res<-mol_assoc_summary_data()

      # Filter results table by user specified methods
      all_func_res<-all_func_res[all_func_res$Method %in% input$selected_methods_mol,]

      # Filter results table by user specified expression and protein panels
      if(any(all_func_res$Type == 'Expr.' | all_func_res$Type == 'Splice')){
        all_func_res<-all_func_res[!((all_func_res$Type == 'Expr.' | all_func_res$Type == 'Splice') & !(all_func_res$Panel %in% input$selected_expr_panels_mol)),]
      }
      if(any(all_func_res$Type == 'Protein')){
        all_func_res<-all_func_res[!(all_func_res$Type == 'Protein' & !(all_func_res$Panel %in% input$selected_protein_panels_mol)),]
      }

      # Insert NA rows for all panels and methods so when filtering by gene, all selected panels and methods remain
      na_rows<-all_func_res[!(duplicated(paste0(all_func_res$Panel, all_func_res$Method))),]
      na_rows$ID<-'Placeholder'
      na_rows$Z<-NA
      na_rows$Sig<-NA
      na_rows$Coloc<-NA

      all_func_res<-rbind(na_rows, all_func_res)

      # Filter results table if user specifies high confidence genes only
      if(input$conf_only_mol){
        # Create group variable
        res_group<-paste0(all_func_res$Method,'\n',all_func_res$Type )
        res_group[res_group == 'SNP\nFine-mapping\n']<-'SuSiE'
        res_group[res_group == 'MAGMA\n']<-'MAGMA'
        res_group[res_group == 'Nearest\nGene\n']<-'Nearest\nGene'

        hc_genes<-all_func_res$ID[which(((all_func_res$Sig == T & all_func_res$Coloc == T) | all_func_res$Panel == "SuSie (L=1)") & res_group %in% input$selected_group_hc_mol)]
        all_func_res<-all_func_res[all_func_res$ID %in% c(hc_genes,'Placeholder'),]
      }

      input_genes <- unlist(strsplit(input$geneInput_mol, "[, ]"))
      selected_genes <- input_genes[input_genes != ""]

      if(length(selected_genes) > 0){
        if(sum(grepl(paste(selected_genes, collapse='|'), all_func_res$ID, ignore.case = T)) > 0){
          selected_genes<-c(selected_genes, 'Placeholder')
          all_func_res<-all_func_res[grepl(paste(selected_genes, collapse='|'), all_func_res$ID, ignore.case = T) & !is.na(all_func_res$ID),]
        } else {
          all_func_res<-data.frame(matrix(nrow=0, ncol=5))
        }
      }

      return(all_func_res)
    }) %>% debounce(1000)

    # Identify number of genes
    plot_dim_mol<-reactive({
      all_func_res<-mol_assoc_summary_data_filtered()

      if(nrow(all_func_res) > 0){
        num_row <- length(unique(all_func_res$ID))
        plot_height<-(max(nchar(all_func_res$Panel))*3)+(num_row * 20)+100
        num_col <- length(unique(paste0(all_func_res$Panel,'_',all_func_res$Method,'_',all_func_res$Type)))
        num_pan <- length(unique(all_func_res$Method))
        plot_width<-120+(max(nchar(all_func_res$ID), na.rm=T)*4) + (num_col * 27) + (num_pan*15)
        plot_width<-max(plot_width,(length(unique(all_func_res$Method))*140))
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
        all_func_res_all<-all_func_res_all[all_func_res_all$ID != 'Placeholder',]

        all_func_res_all$Group<-paste0(all_func_res_all$Method,'\n',all_func_res_all$Type )
        all_func_res_all$Group[all_func_res_all$Group == 'SNP\nFine-mapping\n']<-'SuSiE'
        all_func_res_all$Group[all_func_res_all$Group == 'MAGMA\n']<-'MAGMA'
        all_func_res_all$Group[all_func_res_all$Group == 'Nearest\nGene\n']<-'Nearest\nGene'

        groups<-c('SuSiE','FUSION\nExpr.','FUSION\nSplice','SMR\nExpr.','FUSION\nProtein','SMR\nProtein','MAGMA','Nearest\nGene')
        groups<-groups[groups %in% all_func_res_all$Group]

        all_func_res_all$Group<-factor(all_func_res_all$Group, levels=groups)

        all_func_res_all<-all_func_res_all[order(as.character(all_func_res_all$ID)),]

        all_func_res_all$ID<-factor(all_func_res_all$ID, levels=unique(all_func_res_all$ID))

        group_siz<-NULL
        for(i in groups){
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

      tmp<-gwas_data()[[selected_gwas()]]$tx$drug$magma
      if(is.null(tmp)) return(NULL)
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

      tmp<-gwas_data()[[selected_gwas()]]$tx$drug$gcsc
      if(is.null(tmp)) return(NULL)
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

      tmp<-gwas_data()[[selected_gwas()]]$tx$drug$twas_gsea
      if(is.null(tmp)) return(NULL)
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
      magma_gs<-gwas_data()[[selected_gwas()]]$tx$drug$magma

      if(!is.null(magma_gs)){
        # Convert one-sided p to a Z score
        magma_gs$Z<--qnorm(magma_gs$P)
        magma_gs<-magma_gs[,c('Name','Z','P','P.FDR','ATC Code')]
        magma_gs$Method<-'MAGMA'
        magma_gs$Panel<-'MAGMA'
      }

      ###
      # GCSC
      ###
      gcsc_gs<-gwas_data()[[selected_gwas()]]$tx$drug$gcsc

      if(!is.null(gcsc_gs)){
        gcsc_gs<-gcsc_gs[,c('Name','Z','P','P.FDR','ATC Code')]
        gcsc_gs$Method<-'GCSC'
        gcsc_gs$Panel<-'Brain and Blood'
      }

      ###
      # TWAS-GSEA
      ###

      gsea_gs<-gwas_data()[[selected_gwas()]]$tx$drug$twas_gsea

      if(!is.null(gsea_gs)){
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
      }

      ###
      # Combine results
      ###

      all_gs<-do.call(rbind, Filter(Negate(is.null), list(magma_gs, gcsc_gs, gsea_gs)))

      return(all_gs)
    })

    observeEvent(tx_drug_summary_data(), {
      all_gs<-tx_drug_summary_data()

      methods<-unique(all_gs$Method)
      updateSelectInput(session, "selected_methods_drug", choices = methods, selected=methods)

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

      # Insert NA rows for all panels and methods so when filtering by drug, all selected panels and methods remain
      na_rows<-all_gs[!(duplicated(paste0(all_gs$Panel, all_gs$Method))),]
      na_rows$Name<-'Placeholder'
      na_rows$`ATC Code`<-'Placeholder'
      na_rows$Z<-NA
      na_rows$P<-NA
      na_rows$P.FDR<-NA

      all_gs<-rbind(na_rows, all_gs)

      # Filter results table if user specifies high confidence genes only
      if(input$conf_only_drug){
        hc_drugs<-all_gs$Name[
          which(
            (all_gs$Method == 'TWAS-GSEA' & all_gs$P.FDR < 0.05) |
              (all_gs$Method == 'MAGMA' & all_gs$P.FDR < 0.05 & all_gs$Z > 0) |
              (all_gs$Method == 'GCSC' & all_gs$P.FDR < 0.05 & all_gs$Z > 0)
          )]

        all_gs<-all_gs[all_gs$Name %in% c(hc_drugs,'Placeholder'),]
      }

      input_drugs <- unlist(strsplit(input$drugInput_drug, "[, ]"))
      selected_drugs <- input_drugs[input_drugs != ""]

      input_atc <- unlist(strsplit(input$atcInput_drug, "[, ]"))
      selected_atc <- input_atc[input_atc != ""]

      if(length(selected_drugs) > 0){
        if(sum(grepl(paste(selected_drugs, collapse='|'), all_gs$Name, ignore.case = T)) > 0){
          selected_drugs<-c(selected_drugs, 'Placeholder')
          all_gs<-all_gs[grepl(paste(selected_drugs, collapse='|'), all_gs$Name, ignore.case = T) & !is.na(all_gs$Name),]
        } else {
          all_gs<-data.frame(matrix(nrow=0, ncol=5))
        }
      }

      if(length(selected_atc) > 0){
        if(sum(grepl(paste(selected_atc, collapse='|'), all_gs$`ATC Code`, ignore.case = T)) > 0){
          selected_atc<-c(selected_atc, 'Placeholder')
          all_gs<-all_gs[grepl(paste(selected_atc, collapse='|'), all_gs$`ATC Code`, ignore.case = T) & !is.na(all_gs$`ATC Code`),]
        } else {
          all_gs<-data.frame(matrix(nrow=0, ncol=5))
        }
      }

      return(all_gs)
    }) %>% debounce(1000)

    # Identify number of drugs
    plot_dim_drug<-reactive({
      all_gs<-tx_drug_summary_data_filtered()

      if(nrow(all_gs) > 0){
        num_row <- length(unique(all_gs$Name))
        plot_height<-(max(nchar(all_gs$Panel))*3)+(num_row * 20)+100
        num_col <- length(unique(paste0(all_gs$Panel,'_',all_gs$Method,'_')))
        num_pan <- length(unique(all_gs$Method))
        plot_width<-120+(max(nchar(all_gs$Name), na.rm=T)*4)+(num_col * 27) + (num_pan*15)
        plot_width<-max(plot_width,(length(unique(all_gs$Method))*140))
      } else {
        plot_height<-100
        plot_width<-100
      }

      return(list(height=plot_height,
                  width=plot_width))
    })

    observeEvent(tx_drug_summary_data_filtered(), {
      tmp<-tx_drug_summary_data_filtered()
      choices<-'All - Z'
      if(any(tmp$Method == 'MAGMA')){
        choices<-c(choices, 'MAGMA - Z')
      }
      if(any(tmp$Method == 'GCSC')){
        choices<-c(choices, 'GCSC - Z')
      }
      if(any(tmp$Method == 'TWAS-GSEA')){
        choices<-c(choices, 'TWAS-GSEA - Z')
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
        all_gs_all<-all_gs_all[all_gs_all$Name != 'Placeholder',]

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
          facet_wrap(~ Method , nrow=1, scales = "free_x") +
          theme(text = element_text(size = 14))

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

      tmp<-gwas_data()[[selected_gwas()]]$tx$atc$magma
      if(is.null(tmp)) return(NULL)
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

      tmp<-gwas_data()[[selected_gwas()]]$tx$atc$gcsc
      if(is.null(tmp)) return(NULL)
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
      tmp<-gwas_data()[[selected_gwas()]]$tx$atc$twas_gsea
      if(is.null(tmp)) return(NULL)

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

      magma_gs_atc<-gwas_data()[[selected_gwas()]]$tx$atc$magma

      if(!is.null(magma_gs_atc)){
        magma_gs_atc$Z<--qnorm(magma_gs_atc$P)
        magma_gs_atc$FDR_Sig<-magma_gs_atc$P.FDR < 0.05
        magma_gs_atc$Nom_Sig<-magma_gs_atc$P < 0.05
        magma_gs_atc$Name<-paste0(magma_gs_atc$`ATC Code`,': ',magma_gs_atc$`ATC Description`)
        magma_gs_atc$Method<-'MAGMA'
        magma_gs_atc$Panel<-'MAGMA'

        magma_gs_atc<-magma_gs_atc[,c("Name","Z","FDR_Sig","Nom_Sig","Method","Panel"), with=F]
      }

      ###
      # GCSC
      ###

      gcsc_gs_atc<-gwas_data()[[selected_gwas()]]$tx$atc$gcsc

      if(!is.null(gcsc_gs_atc)){
        gcsc_gs_atc$Z<--qnorm(gcsc_gs_atc$P)
        gcsc_gs_atc$FDR_Sig<-gcsc_gs_atc$P.FDR < 0.05
        gcsc_gs_atc$Nom_Sig<-gcsc_gs_atc$P < 0.05
        gcsc_gs_atc$Name<-paste0(gcsc_gs_atc$`ATC Code`,': ',gcsc_gs_atc$`ATC Description`)
        gcsc_gs_atc$Method<-'GCSC'
        gcsc_gs_atc$Panel<-'GCSC'

        gcsc_gs_atc<-gcsc_gs_atc[,c("Name","Z","FDR_Sig","Nom_Sig","Method","Panel"), with=F]
      }

      ###
      # TWAS-GSEA
      ###

      gsea_gs_atc<-gwas_data()[[selected_gwas()]]$tx$atc$twas_gsea

      if(!is.null(gsea_gs_atc)){
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
      }

      all_gs_atc<-do.call(rbind, Filter(Negate(is.null), list(magma_gs_atc, gcsc_gs_atc, gsea_gs_atc)))

      return(all_gs_atc)
    })

    observeEvent(tx_atc_summary_data(), {
      all_gs<-tx_atc_summary_data()

      methods<-unique(all_gs$Method)
      updateSelectInput(session, "selected_methods_atc", choices = methods, selected=methods)

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

      # Insert NA rows for all panels and methods so when filtering by atc, all selected panels and methods remain
      na_rows<-all_gs_atc[!(duplicated(paste0(all_gs_atc$Panel, all_gs_atc$Method))),]
      na_rows$Name<-'Placeholder'
      na_rows$Z<-NA
      na_rows$FDR_Sig<-NA
      na_rows$Nom_Sig<-NA

      all_gs_atc<-rbind(na_rows, all_gs_atc)

      # Filter results table if user specifies high confidence genes only
      if(input$conf_only_atc){
        hc_atc<-all_gs_atc$Name[
          which(
            (all_gs_atc$Method == 'TWAS-GSEA' & all_gs_atc$FDR_Sig) |
              (all_gs_atc$Method == 'MAGMA' & all_gs_atc$FDR_Sig & all_gs_atc$Z > 0) |
              (all_gs_atc$Method == 'GCSC' & all_gs_atc$FDR_Sig & all_gs_atc$Z > 0)
          )]

        all_gs_atc<-all_gs_atc[all_gs_atc$Name %in% c(hc_atc,'Placeholder'),]
      }

      input_atcs <- unlist(strsplit(input$atcInput_atc, "[, ]"))
      selected_atcs <- input_atcs[input_atcs != ""]

      if(length(selected_atcs) > 0){
        if(sum(grepl(paste(selected_atcs, collapse='|'), all_gs_atc$Name, ignore.case = T)) > 0){
          selected_atcs<-c(selected_atcs, 'Placeholder')
          all_gs_atc<-all_gs_atc[grepl(paste(selected_atcs, collapse='|'), all_gs_atc$Name, ignore.case = T) & !is.na(all_gs_atc$Name),]
        } else {
          all_gs_atc<-data.frame(matrix(nrow=0, ncol=5))
        }
      }

      return(all_gs_atc)
    }) %>% debounce(1000)

    # Identify number of atcs
    plot_dim_atc<-reactive({
      all_gs_atc<-tx_atc_summary_data_filtered()

      if(nrow(all_gs_atc) > 0){
        num_row <- length(unique(all_gs_atc$Name))
        plot_height<-(max(nchar(all_gs_atc$Panel))*4)+(num_row * 20)+100
        num_col <- length(unique(paste0(all_gs_atc$Panel,'_',all_gs_atc$Method,'_')))
        num_pan <- length(unique(all_gs_atc$Method))
        plot_width<-120+(max(nchar(all_gs_atc$Name), na.rm=T)*4)+(num_col * 27) + (num_pan*15)
        plot_width<-max(plot_width,(length(unique(all_gs_atc$Method))*140))
      } else {
        plot_height<-100
        plot_width<-100
      }

      return(list(height=plot_height,
                  width=plot_width))
    })

    observeEvent(tx_atc_summary_data_filtered(), {
      tmp<-tx_atc_summary_data_filtered()
      choices<-'All - Z'
      if(any(tmp$Method == 'MAGMA')){
        choices<-c(choices, 'MAGMA - Z')
      }
      if(any(tmp$Method == 'GCSC')){
        choices<-c(choices, 'GCSC - Z')
      }
      if(any(tmp$Method == 'TWAS-GSEA')){
        choices<-c(choices, 'TWAS-GSEA - Z')
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
        all_gs_atc_all<-all_gs_atc_all[all_gs_atc_all$Name != 'Placeholder',]

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
          scale_colour_gradientn(colours=c("#0066FF","#0099FF","#FFFFFF","#FF6666","#FF0000"), na.value = NA,name = "Z-score", limits = c(-max(abs(all_gs_atc_all$Z), na.rm=T), max(abs(all_gs_atc_all$Z), na.rm=T)), values=x) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
          labs(x='', y='') +
          facet_wrap(~ Method , nrow=1, scales = "free_x") +
          theme(text = element_text(size = 14))


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

    #######################
    # References
    #######################

    output$reference_table<-renderDataTable({
      ref_table<-NULL

      ref_table<-rbind(ref_table, data.frame(Name = 'LD Score Regression',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/25642630">25642630</a>',
                                             Website = '<a href="https://github.com/bulik/ldsc">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate SNP-h2 and LD-score intercept'))

      ref_table<-rbind(ref_table, data.frame(Name = 'SuSiE',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/35853082">35853082</a>',
                                             Website = '<a href="https://stephenslab.github.io/susie-paper">link</a>',
                                             Type = 'Software',
                                             Use = 'Variant-level finemapping'))

      ref_table<-rbind(ref_table, data.frame(Name = 'FUSION',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/26854917">26854917</a>',
                                             Website = '<a href="http://gusevlab.org/projects/fusion/">link</a>',
                                             Type = 'Software',
                                             Use = 'Perform TWAS/PWAS'))

      ref_table<-rbind(ref_table, data.frame(Name = 'COLOC',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/24830394">24830394</a>',
                                             Website = '<a href="https://chr1swallace.github.io/coloc/">link</a>',
                                             Type = 'Software',
                                             Use = 'Used for colocalisation within FUSION TWAS/PWAS analysis'))

      ref_table<-rbind(ref_table, data.frame(Name = 'SMR',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/29763751">29763751</a>',
                                             Website = '<a href="https://yanglab.westlake.edu.cn/software/smr/#Overview">link</a>',
                                             Type = 'Software',
                                             Use = 'Infer differential expression or protein levels'))

      ref_table<-rbind(ref_table, data.frame(Name = 'MAGMA',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/25885710">25885710</a>',
                                             Website = '<a href="https://ctg.cncr.nl/software/magma">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate gene associations and drug enrichment using DrugTargetor database'))

      ref_table<-rbind(ref_table, data.frame(Name = 'TWAS-GSEA',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/31230729">31230729</a>',
                                             Website = '<a href="https://github.com/opain/TWAS-GSEA">link</a>',
                                             Type = 'Software',
                                             Use = 'Estimate drug enrichment using TWAS and DrugTargetor'))

      ref_table<-rbind(ref_table, data.frame(Name = 'PsychENCODE DLPFC TWAS SNP-weights',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/30545856">30545856</a>',
                                             Website = '<a href="http://resource.psychencode.org/">link</a>',
                                             Type = 'Dataset',
                                             Use = 'TWAS'))

      ref_table<-rbind(ref_table, data.frame(Name = 'PsychENCODE DLPFC eQTL data for SMR',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/30545857">30545857</a>',
                                             Website = '<a href="http://cnsgenomics.com/data/SMR/">link</a>',
                                             Type = 'Dataset',
                                             Use = 'SMR'))

      ref_table<-rbind(ref_table, data.frame(Name = 'MetaBrain eQTL data for SMR',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/36823318">36823318</a>',
                                             Website = '<a href="https://www.metabrain.nl/cis-eqtls.html">link</a>',
                                             Type = 'Dataset',
                                             Use = 'SMR'))

      ref_table<-rbind(ref_table, data.frame(Name = 'ROSMAP/Banner DLPFC PWAS SNP-weights and eQTL in SMR results',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/33510477">33510477</a>',
                                             Website = '<a href="https://www.synapse.org/#!Synapse:syn23627957">link</a>',
                                             Type = 'Dataset',
                                             Use = 'PWAS and SMR'))

      ref_table<-rbind(ref_table, data.frame(Name = 'DrugTargetor',
                                             PMID = '<a href="https://pubmed.ncbi.nlm.nih.gov/30517594">30517594</a>',
                                             Website = '<a href="https://drugtargetor.com/">link</a>',
                                             Type = 'Dataset',
                                             Use = 'Drug repurposing'))

      datatable(ref_table,
                rownames= FALSE,
                escape = -2:-3,
                selection = 'none',
                options = list(
                  paging = FALSE,
                  dom = 'lrt'
                )
                )
    })

    ###############################
    # Configuration
    ###############################

    output$repo_info <- renderUI({
      tagList(
        p(HTML(paste0('<strong>Repo:</strong> ', gwas_data()$configuration$repo$remote))),
        p(HTML(paste0('<strong>Branch:</strong> ', gwas_data()$configuration$repo$branch))),
        p(HTML(paste0('<strong>Commit:</strong> ', gwas_data()$configuration$repo$commit)))
      )
    })

    output$config_table <- DT::renderDataTable({
      # Make table showing config parameters
      config_tab<-config[!grepl('#', config)]
      config_tab<-config_tab[config_tab != '']
      config_tab<-data.frame(do.call(rbind, strsplit(config_tab, ': ')))
      names(config_tab)<-c('Parameter', 'Value')

      datatable(config_tab,
                rownames = FALSE,
                options = list(
                  scrollX = TRUE,
                  ordering = FALSE,
                  columnDefs = list(
                    list(className = "dt-left", targets = "_all"),
                    list(width = '250px', targets = 0)
                  )
                ),
                selection = 'none')
    })

    output$gwas_list <- DT::renderDataTable({

      dat<-gwas_list
      dat[is.na(dat)]<-'NA'

      datatable(dat,
                rownames = FALSE,
                options = list(
                  scrollX = TRUE,
                  dom = 't',
                  ordering = FALSE,
                  columnDefs = list(
                    list(className = "dt-left", targets = "_all")  # Apply the class to all columns
                  )
                ),
                selection = 'none')
    })
  })
}

# Run the Shiny app
shinyApp(ui, server)
