configurationUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    title="Configuration",

    br(),
    div(
      class = "custom-panel",
      h6(strong('Repository information:')),
      uiOutput(ns("repo_info"))
    ),
    hr(),
    div(
      class = "custom-panel",
      h6(strong('Config file:')),
      DT::dataTableOutput(ns("config_table")),
    ),
    hr(),
    div(
      class = "custom-panel",
      h6(strong('GWAS list:')),
      DT::dataTableOutput(ns("gwas_list"))
    ),
    hr()
  )
}

configurationServer <- function(id, gwas_data) {
  moduleServer(id, function(input, output, session) {

    output$repo_info <- renderUI({
      req(gwas_data())
      tagList(
        p(HTML(paste0('<strong>Repo:</strong> ', gwas_data()$configuration$repo$remote))),
        p(HTML(paste0('<strong>Branch:</strong> ', gwas_data()$configuration$repo$branch))),
        p(HTML(paste0('<strong>Commit:</strong> ', gwas_data()$configuration$repo$commit)))
      )
    })

    output$config_table <- DT::renderDataTable({
      req(gwas_data())
      # Make table showing config parameters
      config<-gwas_data()$configuration$config
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
      req(gwas_data())
      gwas_list<-gwas_data()$configuration$gwas_list

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
