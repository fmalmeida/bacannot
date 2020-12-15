#' @import graphics stats shiny methods
#' @export

# Parse sql
parse_sql <- function(sqldb) {
  
  shinyApp(
    ui = fluidPage(
      mainPanel(
        fluidRow(
          #column(3, selectInput("sqltable", "Available SQL Tables:",
          #                      c("Annotation information" = "GFF",
          #                        "Protein gene sequences" = "ProteinFasta",
          #                        "Nucleotide gene sequences" = "NucleotideFasta",
          #                        "Genome summary" = "Genome"))),
          column(2, uiOutput("contigs")),
          column(2, uiOutput("sources")),
          column(2, uiOutput("features")),
          column(2, textInput("start", "Min. start", value = "0")),
          column(2, textInput("end", "Max. end"))
        ),
        fluidRow(
          column(12, dataTableOutput("table"))
        )
      )
    ),
    
    server = function(input, output, session) {
      
      # Create sql string from user selection
      #sqlString <- reactive ({
      #  paste0("SELECT * FROM ", input$sqltable)
      #})
      
      # Open DB in the dataframe
      openTable <- reactive({
        as.data.frame(tbl(sqldb, sql(paste0("SELECT * FROM GFF"))))
      })
      
      # Create selection for contigs
      output$contigs <- renderUI({
        selectInput("contigs", "Select contigs:",
                    unique(sort(openTable()[,1])), 
                    multiple = TRUE)
      })
      
      # Create selection for sources
      output$sources <- renderUI({
        selectInput("sources", "Select sources:",
                    sort(c("CARD", "AMRFinderPlus", "barrnap", "PHAST", "ICEberg", "VFDB", "Victors", "Prodigal")), 
                    multiple = TRUE)
      })
      
      # Create selection for features
      output$features <- renderUI({
        selectInput("features", "Select features:",
                    sort(c("CDS", "tRNA", "rRNA", "ICE", "Resistance", "Virulence", "Prophage")), 
                    multiple = TRUE)
      })
      
      # Grabing users filter selections
      
      ## Contig Filter
      contigFilter <- reactive({
        if (length(as.list(input$contigs)) > 0) {
          as.list(input$contigs)
        } else {
          unique(sort(openTable()[,1]))
        }
      })
      
      ## Source Filter
      sourceFilter <- reactive({
        if (length(as.list(input$sources)) > 0) {
          as.list(input$sources)
        } else {
          sort(c("CARD", "AMRFinderPlus", "barrnap", "PHAST", "ICEberg", "VFDB", "Victors", "Prodigal"))
        }
      })
      
      ## Feature Filter
      featureFilter <- reactive({
        if (length(as.list(input$features)) > 0) {
          as.list(input$features)
        } else {
          sort(c("CDS", "tRNA", "rRNA", "ICE", "Resistance", "Virulence", "Prophage"))
        }
      })
      
      ## Start Filter
      startFilter <- reactive({
        if (isTruthy(input$start)) {
          as.integer(input$start)
        } else {
          as.integer(0)
        }
      })
      
      ## End Filter
      endFilter <- reactive({
        if (isTruthy(input$end)) {
          as.integer(input$end)
        } else {
          Inf
        }
      })
      
      # Render table
      output$table <- 
        renderDataTable({
            datatable(openTable() %>%
                        filter(Contig %in% contigFilter()) %>%
                        filter(str_detect(Source, paste0(sourceFilter(), collapse = "|"))) %>%
                        filter(str_detect(Feature, paste0(featureFilter(), collapse = "|"))) %>%
                        filter(as.integer(Start) >= startFilter()) %>%
                        filter(as.integer(End) <= endFilter()) %>%
                        arrange(ID),
                      rownames = F, selection = 'none',
                      options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15, 20, 50),
                                     scrollX = FALSE,
                                     fixedColumns = FALSE,
                                     autoWidth = TRUE,
                                     searching= FALSE)) %>%
            formatStyle(names(openTable()), lineHeight='80%') 
        })
      
    },
    
    options=list(width = "100%", height = 700)
  )
}