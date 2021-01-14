# Shiny app to blast the genome, nucleotide and protein sequences
# used in bacannot pipeline and parser
genome_blast <- function(genome, nucl, prot, sqldb) {
  
  shinyApp(
    options=list(height = '2000px', width = '100%'),
    ui <- fluidPage(
      
      # Custom CSS
      tags$head(
        tags$style(HTML("
      body { 
        /* Normal  */
        font-size: 12pt;
        line-height: 1.5;
      }
        
      blockquote {
          padding: 10px 20px;
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
      
      # Choose shiny theme
      theme = shinytheme("cosmo"),
      
      # This block gives us all the inputs:
      mainPanel(
        width = "100%",
        h2('Blast the genome!'),
        br(),
        textAreaInput('query', 'Input sequence:', value = "", placeholder = "", height="200px", width = '850px'),
        fluidRow(
          column(2, selectInput("db", "Database:", choices=c("genome","genes-nt", "genes-aa"))),
          column(2, selectInput("program", "Program:", choices=c("blastn", "tblastn", "blastx", "blastp"))),
          column(2, selectInput("eval", "e-value:", choices=c(1,0.001,1e-4,1e-5,1e-10)))
        ),
        actionButton("blast", "BLAST!"),
        h3("Blast results"),
        br(),
        actionButton("reset", "Reset row selection"),
        br(),
        dataTableOutput("blastResults", width = '100%', height = 'auto'),
        selectInput("format", "Download format:",
                    choices = c('tsv', 'csv', 'spreadsheet')),
        downloadButton("downloadData", "Download"),
        h3("Annotation intersection!"),
        br(),
        actionButton("reset2", "Reset row selection"),
        br(),
        dataTableOutput("intersection", width = '100%', height = 'auto'),
        selectInput("format2", "Download format:",
                    choices = c('gff', 'fasta-nucl', 'fasta-prot','spreadsheet')),
        downloadButton("downloadData2", "Download"),
      )
    ),
    
    server = function(input, output, session) {
      
      #######################################################
      ### Step 0 - create a watcher for the blast results ###
      #######################################################
      # open dbs in data.frames
      openTable <- reactive({
        as.data.frame(tbl(sqldb, sql(paste0("SELECT * FROM GFF")))) %>%
          select(-ID)
      })
      
      openTableNucl <- reactive({
        as.data.frame(tbl(sqldb, sql(paste0("SELECT * FROM NucleotideFasta"))))
      })
      
      openTableProt <- reactive({
        as.data.frame(tbl(sqldb, sql(paste0("SELECT * FROM ProteinFasta"))))
      })
      
      # Reactive values
      rv <- reactiveValues()
      rv$data <- NULL
      
      #####################################################
      ### Step 1- Grab user input sequence for BLASTing ###
      ### and the required BLAST parameters             ###
      #####################################################
      observeEvent(input$blast, {
        
        # Save query and tmp file
        query <- input$query
        tmp <- tempfile(fileext = ".fa")
        
        # Load query seqs in tmp file
        if (startsWith(query, ">")){
          writeLines(query, tmp)
        } else {
          writeLines(paste0(">Query\n", query), tmp)
        }
        
        # Check user desired database
        if (input$db == "genome") {
          subject <- genome
        } else if (input$db == "genes-nt") {
          subject <- nucl
        } else if (input$db == "genes-aa") {
          subject <- prot
        }
        
        ##########################################
        ### Step 2 - BLAST the query sequence! ###
        ##########################################
        showModal(modalDialog("Executing your BLAST, please wait!", footer=NULL))
        rv$data <- read.table(
          text = system(
            paste0(input$program, " -query ", tmp, " -subject ", subject, " -evalue ", input$eval, " -outfmt 6 ", collapse = ""),
            intern = TRUE
          ), col.names = c("query_id", "subject_id", "pct_identity", "aln_length", "n_of_mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"), comment.char = "#")
        removeModal()
        
        # Call event
        rv$data
        
      }, ignoreNULL = TRUE)
      
      ###################################################
      ### Step 3 - Print the BLAST results to a table ###
      ###################################################
      output$blastResults <- renderDataTable({
        
        data <- rv$data
        data[,1] <- as.factor(data[,1])
        datatable(data,
                  rownames = F,
                  selection = 'multiple',
                  filter = 'none',
                  options = list(pageLength = 5,
                                 lengthMenu = c(5, 10, 15, 20, 50),
                                 scrollX = FALSE,
                                 scrollY = FALSE,
                                 fixedColumns = FALSE,
                                 autoWidth = TRUE,
                                 searching= FALSE))
        
      })
      
      ###############################################
      ### Step 4 - Understand users row selection ###
      ###############################################
      selectedData <- reactive({
        
        if(length(input$blastResults_rows_selected) > 0) {
          row_count <- input$blastResults_rows_selected
          data <- rv$data
          data <- data[row_count, ]
          data
        } else {
          data <- rv$data
          data
        }
        
      })
      
      ##################################################
      ### Step 5 - Reset row selection, if necessary ###
      ##################################################
      myProxy = DT::dataTableProxy('blastResults')
      observeEvent(input$reset, {
        DT::selectRows(myProxy, NULL)
      })
      
      ##############################################
      ### Step 6 - Handle BLAST results download ###
      ##############################################
      output$downloadData <- downloadHandler(
        
        filename = function() {
          if(input$format == 'tsv') {
            paste('data_download', '.tsv', sep = '')
          } else if(input$format == 'csv') {
            paste('data_download', '.csv', sep = '')
          } else if(input$format == 'spreadsheet') {
            paste('data_download', '.xlsx', sep = '')
          }
        },
        content = function(file) {
          if(input$format == 'tsv') {
            write.table(selectedData(), file, row.names = F, sep = "\t", quote = F)
          } else if(input$format == 'csv') {
            write.csv(selectedData(), file, row.names = F, sep = ",", quote = F)
          } else if(input$format == 'spreadsheet') {
            write.xlsx(selectedData(), file, sheetName="Sheet1", 
                       col.names=T, row.names=F, append=F)
          }
        }
      )
      
      ##########################################################################
      ### Step 7 - Detect intersection between blast hits and the annotation ###
      ##########################################################################
      annotation_intersection <- eventReactive(input$blast, {
        
        # Default
        out <- NULL
        
        # Requirement
        req(rv$data)
        
        # write gff in tmp file
        file.create("temp.gff")
        write.table(openTable(), "temp.gff", row.names = F, sep = "\t", quote = F, col.names = F)
        
        # BLAST BED
        if (input$db == 'genome') {
          
          # Get bed
          in_bed <- selectedData() %>%
            select(2, 9, 10) %>%
            as.data.frame()
          
          # Parse bed to fix start and end
          parsed_bed <- read.table(
            text = apply(in_bed, 1, function(x) {
            
            if (as.integer(x[2]) >= as.integer(x[3])) {
              st_bp = x[3]
              ed_bp = x[2]
            } else {
              st_bp = x[2]
              ed_bp = x[3]
            }
            
            paste(x[1], st_bp, ed_bp, sep = '\t')
            
          }), sep = '\t')
          
          # save blast as bed
          file.create("temp.bed")
          write.table(parsed_bed, "temp.bed", row.names = F, sep = "\t", quote = F, col.names = F)
          
          # Call intersection
          out <- read.table(
            text = system(
            "intersectBed -wa -a temp.gff -b temp.bed",
            intern = TRUE
          ), sep = '\t',
          col.names = c("Contig", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
          
        } else {
          
          # Get subject ids
          ids <- selectedData() %>%
            select(2) %>%
            as.data.frame()
          
          # Write this IDs in files
          file.create("temp_ids.txt")
          write.table(ids, "temp_ids.txt", row.names = F, sep = "\t", quote = F, col.names = F)
          
          # Use grep to filter the GFF based on the IDs
          out <- read.table(
            text = system(
              "grep -w -f temp_ids.txt temp.gff",
              intern = TRUE
            ), sep = '\t',
            col.names = c("Contig", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes")
          )
        }
        
        # Call result
        out
        
      })

      output$intersection <- renderDataTable({
        
        # Requirement
        req(annotation_intersection())

        # Call datatable
        datatable(annotation_intersection() %>%
                    unique(),
                  rownames = F,
                  selection = 'multiple',
                  filter = 'none',
                  options = list(pageLength = 5,
                                 lengthMenu = c(5, 10, 15, 20, 50),
                                 scrollX = FALSE,
                                 scrollY = FALSE,
                                 fixedColumns = FALSE,
                                 autoWidth = TRUE,
                                 searching= FALSE))

      })
      
      ################################################################################
      ### Step 8 - Understand users row selection in bedtools intersection results ###
      ################################################################################
      selectedData2 <- reactive({
        if(length(input$intersection_rows_selected) > 0) {
          row_count <- input$intersection_rows_selected
          data <- annotation_intersection()[row_count, ]
          data
        } else {
          annotation_intersection()
        } 
      })
      
      #################################################################
      ### Reset row selection, if necessary (Bedtools intersection) ###
      #################################################################
      myProxy = DT::dataTableProxy('intersection')
      observeEvent(input$reset2, {
        DT::selectRows(myProxy, NULL)
      })
      
      ###########################################################################################
      ### Step 9 - Use selected genes to subset the NUCLEOTIDE and PROTEIN sequence databases ###
      ###########################################################################################
      
      # Get NUCL
      filteredNUCL <-
        reactive({
          ids <- getAttributeField(as.character(selectedData2()[,9]), "ID", ";")
          nucl <- openTableNucl() %>%
            filter(ID %in% ids) %>%
            arrange(ID)
          
          # Return FASTA
          glue('>{nucl$ID} {nucl$Description}\n{nucl$Sequence}')
        })
      
      # Get PROT
      filteredPROT <-
        reactive({
          ids <- getAttributeField(as.character(selectedData2()[,9]), "ID", ";")
          prot <- openTableProt() %>%
            filter(ID %in% ids) %>%
            arrange(ID)
          
          # Return FASTA
          glue('>{prot$ID} {prot$Description}\n{prot$Sequence}')
        })
      
      ########################################################################################
      ### Step 10 - Use selected genes to create the parsed attributes col for spreadsheet ###
      ########################################################################################
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
          matrix(ncol = 26, nrow = nrow(selectedData2()))), fields
        ) 
        
        # Apply the getAttributes to each entry
        for(entry in fields) {
          parsed_attributes[, as.character(entry)] <- 
            getAttributeField(as.character(selectedData2()[,9]), as.character(entry), ";")
        }
        
        # Output
        parsed_attributes
      })
      
      
      #################################################################
      ### Final step, handle data download of bedtools intersection ###
      #################################################################
      output$downloadData2 <- downloadHandler(
        
        filename = function() {
          if(input$format2 == 'gff') {
            paste('data_download', '.gff', sep = '')
          } else if(input$format2 == 'fasta-nucl') {
            paste('data_download', '.fna', sep = '')
          } else if(input$format2 == 'fasta-prot') {
            paste('data_download', '.faa', sep = '')
          } else if(input$format2 == 'spreadsheet') {
            paste('data_download', '.xlsx', sep = '')
          }
        },
        content = function(file) {
          if(input$format2 == 'gff') {
            write.table(selectedData2(), file, row.names = F, sep = "\t", quote = F, col.names = F)
          } else if(input$format2 == 'fasta-nucl') {
            write(filteredNUCL(), file, sep = "\n")
          } else if(input$format2 == 'fasta-prot') {
            write(filteredPROT(), file, sep = "\n")
          } else if(input$format2 == 'spreadsheet') {
            write.xlsx(attributes(), file, sheetName="Parsed_Attributes", 
                       col.names=T, row.names=F, append=F)
          }
        }
      )
      
    }
  )
  
}