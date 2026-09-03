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
      conf <- gd_config(gwas_data())
      # New-format bundles only carry pipeline_version (no repo remote/branch);
      # legacy .rds files carry the old remote/branch/commit trio.
      rows <- list(
        p(HTML(paste0('<strong>Pipeline version:</strong> ',
                      if (!is.na(conf$pipeline_version)) conf$pipeline_version else '(unknown)')))
      )
      if (!is.na(conf$remote) && nzchar(conf$remote)) rows <- c(rows, list(p(HTML(paste0('<strong>Repo:</strong> ',   conf$remote)))))
      if (!is.na(conf$branch) && nzchar(conf$branch)) rows <- c(rows, list(p(HTML(paste0('<strong>Branch:</strong> ', conf$branch)))))
      if (!is.na(conf$commit) && nzchar(conf$commit) && !identical(conf$commit, conf$pipeline_version))
        rows <- c(rows, list(p(HTML(paste0('<strong>Commit:</strong> ', conf$commit)))))
      do.call(tagList, rows)
    })

    output$config_table <- DT::renderDataTable({
      req(gwas_data())
      # Make table showing config parameters
      config<-gd_config(gwas_data())$flags_raw
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
      gwas_list<-gd_config(gwas_data())$gwas_list

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
