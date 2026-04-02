enrichmentUI <- function(id) {
  ns <- NS(id)

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
                selectInput(ns("selected_methods_drug"), "Select methods", "", multiple=T),
                selectInput(ns("selected_expr_panels_drug"), "Select expression panels", "", multiple=T),
                radioButtons(ns("conf_only_drug"), "Show FDR significant only :",
                             choices = c("True" = T,
                                         "False" = F),
                             selected = T),
                textInput(ns("drugInput_drug"), "Search drug (whitespace- or comma-seperated):"),
                textInput(ns("atcInput_drug"), "Search ATC Code (whitespace- or comma-seperated):"),
                selectInput(ns("selected_sort_drug"), "Sort by:", '', multiple = F)

              ),

              mainPanel(
                uiOutput(ns("message_too_large_drug")),
                uiOutput(ns("message_no_drugs_drug")),
                uiOutput(ns("tx_drug_plot.ui"))
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
                     dataTableOutput(ns("tx_drug_magma_table")),
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
                     dataTableOutput(ns("tx_drug_gcsc_table")),
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
                     dataTableOutput(ns("tx_drug_twas_gsea_table")),
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
                selectInput(ns("selected_methods_atc"), "Select methods", "", multiple=T),
                selectInput(ns("selected_expr_panels_atc"), "Select expression panels", "", multiple=T),
                radioButtons(ns("conf_only_atc"), "Show FDR significant only :",
                             choices = c("True" = T,
                                         "False" = F),
                             selected = T),
                textInput(ns("atcInput_atc"), "Search ATC Code (whitespace- or comma-seperated):"),
                selectInput(ns("selected_sort_atc"), "Sort by:", '', multiple = F)

              ),

              mainPanel(
                uiOutput(ns("message_too_large_atc")),
                uiOutput(ns("message_no_atcs_atc")),
                uiOutput(ns("tx_atc_plot.ui"))
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
                     dataTableOutput(ns("tx_atc_magma_table")),
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
                     dataTableOutput(ns("tx_atc_gcsc_table")),
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
                     dataTableOutput(ns("tx_atc_twas_gsea_table")),
              )
            ),
            br()
          )
        )
      )
    )
  )
}

enrichmentServer <- function(id, gwas_data, selected_gwas) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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

        # Convert to data.table before geom_point subsetting
        all_gs_all <- data.table(all_gs_all)

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
        plotOutput(ns("tx_drug_plot"), height = plot_dim_drug()[['height']], width=plot_dim_drug()[['width']])
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

        # Convert to data.table before geom_point subsetting
        all_gs_atc_all <- data.table(all_gs_atc_all)

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
        plotOutput(ns("tx_atc_plot"), height = plot_dim_atc()[['height']], width=plot_dim_atc()[['width']])
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
