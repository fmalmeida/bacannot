# Shiny GFF/SQLdb parser -- To understand the results of Bacannot
parse_sql <- function(sqldb) {
  
  shinyApp(
    options=list(height = '1500px', width = '100%'),
    ui = fluidPage(
      
      # Custom CSS
      tags$head(
        tags$style(HTML("
      body { 
        /* Normal  */
        font-size: 12pt;
        line-height: 1.5;
      }
        
      blockquote {
          padding: 11px 20px;
          margin: 0 0 20px;
          font-size: 10pt;
          border-left: 5px solid #ec5800;
      }
        
      .table {
        width: 100% !important;
        text-align: center;
        font-size: 10pt;
        margin-left: auto;
        margin-right: auto;
      }

    "))),
      
      # Page theme
      theme = shinytheme("cosmo"),
      
      # Page UIs
      mainPanel(
        h3("Main filters:"),
        fluidRow(
          column(2, uiOutput("contigs")),
          column(2, uiOutput("sources")),
          column(2, uiOutput("features")),
          column(2, textInput("start", "Min. start", value = "0")),
          column(2, textInput("end", "Max. end")),
          column(2, selectInput("strand", "Strand:",
                                sort(c("+", "-")), 
                                multiple = TRUE))
        ),
        h3("Attributes filters:"),
        markdown("> **Obs**: Using a file of patterns users can filter the annotation file based on 
                 the values found in the attributes column. The file must have one pattern per line, 
                 for any of the fields available in the 9th column. Users can download the annotation
                 in excel spreadsheet format to better visualize the values found in the attributes. 
                 The download button is found in the end of the page."),
        fileInput("att_patterns", "Input file of patterns (values)", multiple = FALSE),
        verbatimTextOutput("testeOUT"),
        h3("Results:"),
        markdown("> **Obs**: The rows in the resulting table are selectable for download. 
        By default (with no selection) all rows are used for download. 
        Users can visualize the selected features in the JBrowse by downloading it as 
        a GFF file and inserting it in the Browser by clicking in the button \"Track\", 
        in the Genome Browser."),
        actionButton("reset", "Reset row selection"),
        br(),
        dataTableOutput("table", height = 'auto'), width = "100%",
        selectInput("format", "Download format:",
                    choices = c("gff",  'spreadsheet', "fasta-nucl", "fasta-prot")),
        downloadButton("downloadData", "Download")
      )
    ),
    
    server = function(input, output, session) {
      
      ##############################
      ### Open DBs in dataframes ###
      ##############################
      openTable <- reactive({
        as.data.frame(tbl(sqldb, sql(paste0("SELECT * FROM GFF"))))
      })
      
      openTableNucl <- reactive({
        as.data.frame(tbl(sqldb, sql(paste0("SELECT * FROM NucleotideFasta"))))
      })
      
      openTableProt <- reactive({
        as.data.frame(tbl(sqldb, sql(paste0("SELECT * FROM ProteinFasta"))))
      })
      
      ##############################
      ### Create filter handlers ###
      ##############################
      
      # Create selection for contigs
      output$contigs <- renderUI({
        selectInput("contigs", "Contigs:",
                    unique(sort(openTable()[,1])), 
                    multiple = TRUE)
      })
      
      # Create selection for sources
      output$sources <- renderUI({
        selectInput("sources", "Sources:",
                    sort(c("CARD", "AMRFinderPlus", "barrnap", "PHAST", "ICEberg", "VFDB", "Victors", "Prodigal", "Aragorn")), 
                    multiple = TRUE)
      })
      
      # Create selection for features
      output$features <- renderUI({
        selectInput("features", "Features:",
                    sort(c("CDS", "tRNA", "rRNA", "ICE", "Resistance", "Virulence", "Prophage")), 
                    multiple = TRUE)
      })
      
      #####################################################
      ### Create additional filters based on 9th column ###
      #####################################################
      
      # Parse pattern file
      attributes_patterns <- reactive({
        if (!is.null(input$att_patterns)) {
          file <- input$att_patterns
          paste0(as.list(readLines(file$datapath)), collapse = "|")
        } else {
          ""
        }
      })
      
      #############################
      ### Parse applied filters ###
      #############################
      
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
          sort(c("CARD", "AMRFinderPlus", "barrnap", "PHAST", "ICEberg", "VFDB", "Victors", "Prodigal", "Aragorn"))
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
      
      ## Strand Filter
      strandFilter <- reactive({
        if (length(as.list(input$strand)) > 0) {
          as.list(input$strand)
        } else {
          sort(c("+", "-"))
        }
      })
      
      # Filter data-frame
      filteredGFF <-
        reactive({
          if (!is.null(input$att_patterns)) {
            openTable() %>%
              filter(Contig %in% contigFilter()) %>%
              filter(str_detect(Source, paste0(sourceFilter(), collapse = "|"))) %>%
              filter(str_detect(Feature, paste0(featureFilter(), collapse = "|"))) %>%
              filter(as.integer(Start) >= startFilter()) %>%
              filter(as.integer(End) <= endFilter()) %>%
              filter(Strand %in% strandFilter()) %>%
              arrange(ID) %>%
              filter(str_detect(Attributes, attributes_patterns()))
          } else {
            openTable() %>%
              filter(Contig %in% contigFilter()) %>%
              filter(str_detect(Source, paste0(sourceFilter(), collapse = "|"))) %>%
              filter(str_detect(Feature, paste0(featureFilter(), collapse = "|"))) %>%
              filter(as.integer(Start) >= startFilter()) %>%
              filter(as.integer(End) <= endFilter()) %>%
              filter(Strand %in% strandFilter()) %>%
              arrange(ID)
          }
        })
      
      #############################
      ### Render filtered table ###
      #############################
      output$table <- 
        renderDataTable({
          datatable(filteredGFF() %>%
                      select(ID, Contig, Source, Feature, Start, End, Strand, Attributes),
                    rownames = F, selection = 'multiple', extensions = 'Buttons',
                    options = list(
                      pageLength = 20,
                      lengthMenu = c(5, 10, 15, 20, 50),
                      scrollX = TRUE,
                      scrollY = '600px',
                      fixedColumns = FALSE,
                      autoWidth = TRUE,
                      searching= FALSE,
                      dom = 'lfrtBip',
                      buttons = I('colvis')
                    )
          ) %>% formatStyle(
            columns = colnames(filteredGFF() %>%
                                 select(ID, Contig, Source, Feature, Start, End, Strand, Attributes)), 
            fontSize = '10pt')
        })
      
      ######################################
      ### Understand users row selection ###
      ######################################
      selectedData <- reactive({
        if(length(input$table_rows_selected) > 0) {
          row_count <- input$table_rows_selected
          data <- filteredGFF()[row_count, ]
          data
        } else {
          filteredGFF()
        } 
      })
      
      #########################################
      ### Reset row selection, if necessary ###
      #########################################
      myProxy = DT::dataTableProxy('table')
      observeEvent(input$reset, {
        DT::selectRows(myProxy, NULL)
      })
      
      ##################################################################################
      ### Use selected genes to subset the NUCLEOTIDE and PROTEIN sequence databases ###
      ##################################################################################
      
      # Get NUCL
      filteredNUCL <-
        reactive({
          ids <- getAttributeField(as.character(selectedData()[,10]), "ID", ";")
          nucl <- openTableNucl() %>%
            filter(ID %in% ids) %>%
            arrange(ID)
          
          # Return FASTA
          glue('>{nucl$ID} {nucl$Description}\n{nucl$Sequence}')
        })
      
      # Get PROT
      filteredPROT <-
        reactive({
          ids <- getAttributeField(as.character(selectedData()[,10]), "ID", ";")
          prot <- openTableProt() %>%
            filter(ID %in% ids) %>%
            arrange(ID)
          
          # Return FASTA
          glue('>{prot$ID} {prot$Description}\n{prot$Sequence}')
        })
      
      ##############################################################################
      ### Use selected genes to create the parsed attributes col for spreadsheet ###
      ##############################################################################
      attributes <- reactive({
        
        # attribute fields
        fields <- c("ID", "gene", "product", "inference", "db_xref", "eC_number", "KO",
                    "NDARO:Gene_Name", "NDARO:Gene:Product", "NDARO:Closest_Sequence",
                    "NDARO:Method", "NDARO:Resistance_Category", "NDARO:Resistance:Target",
                    "CARD:Name", "CARD:Product", "CARD:Inference", "CARD:Targeted_drug_class",
                    "VFDB:Product", "VFDB:Target", "victors:Product", "victors:Target",
                    "ICEberg:Product", "ICEberg:Target", "phast:Product", "phast:Target", "note")
        
        # Create empty data.frame for the parsed table
        parsed_attributes <- setNames(data.frame(
          matrix(ncol = 26, nrow = nrow(selectedData()))), fields
          ) 
        
        # Apply the getAttributes to each entry
        for(entry in fields) {
          parsed_attributes[, as.character(entry)] <- 
            getAttributeField(as.character(selectedData()[,10]), as.character(entry), ";")
        }
        
        # Output
        parsed_attributes
      })
      
      
      ########################################
      ### Final step, handle data download ###
      ########################################
      output$downloadData <- downloadHandler(
        
        filename = function() {
          if(input$format == 'gff') {
            paste('data_download', '.gff', sep = '')
          } else if(input$format == 'fasta-nucl') {
            paste('data_download', '.fna', sep = '')
          } else if(input$format == 'fasta-prot') {
            paste('data_download', '.faa', sep = '')
          } else if(input$format == 'spreadsheet') {
            paste('data_download', '.xlsx', sep = '')
          }
        },
        content = function(file) {
          if(input$format == 'gff') {
            write.table(selectedData() %>%
                          select(-ID), file, row.names = F, sep = "\t", quote = F)
          } else if(input$format == 'fasta-nucl') {
            write(filteredNUCL(), file, sep = "\n")
          } else if(input$format == 'fasta-prot') {
            write(filteredPROT(), file, sep = "\n")
          } else if(input$format == 'spreadsheet') {
            write.xlsx(attributes(), file, sheetName="Parsed_Attributes", 
                       col.names=T, row.names=F, append=F)
            # write.xlsx(selectedData() %>%
            #              select(-ID), file, sheetName="Original_GFF", 
            #            col.names=T, row.names=F, append=T)
          }
        }
      )
      
    }
  )
}