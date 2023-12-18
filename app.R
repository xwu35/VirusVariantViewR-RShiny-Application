#----- Source Definitions -----#

source("./global.R")

options(repos = BiocManager::repositories()) # for deployment
options(shiny.sanitize.errors = FALSE) # need to see the error
options(ucscChromosomeNames = FALSE) # for Gvis

#----- Format, Run App -----#

ui <- tagList(
  
  useShinyjs(),
  
  navbarPage(
    
    #theme = shinytheme("flatly"),
    title = "VariantViewR", id = "navbarpage",
    
    tabPanel(title = "Data Set Selection",
           selectInput(inputId = "dataSetSelect",
                       label = "Available Data Sets:",
                       choices = c(dataSet),
                       width = "50%"),

           actionButton(inputId = "go", label = "Go"),
           verbatimTextOutput("buttonValue")),
    
    tabPanel(title = "Sample Table", value = "sampleTab",
             # checkbox for selecting/deselecting all
             #checkboxInput(inputId = "dt_sel", "sel/desel all"),
             #verbatimTextOutput(outputId = "selected_cells", TRUE),
             h4(textOutput("current_dataset")),
             h3("Select from the Sample column to view information:"),
             br(),
             br(),
           DT::dataTableOutput(outputId = "sampleDataTable"),
           downloadButton("downloadData", "Download Sample Table"),
           actionButton(inputId = "inspVariants", 
                        label = "Inspect Variants",
                        icon("fas fa-search")),
           span(actionButton(inputId = "stackCoverage", label = "Stack Coverage",
                             style = "margin-left: 10px;",
                             icon("fas fa-stream"))),
           span(actionButton(inputId = "commonVariants", label = "Find Common Variants",
                             style = "margin-left:  10px",
                             icon("far fa-list-alt"))),
           h6("GetSampleVec()"),
           p(verbatimTextOutput("verbatimOutput1")),
           h6("_cells_selected()"),
           p(verbatimTextOutput("verbatimOutput2")),
           h6("SelectedSampleIndices() (filtered _cells_selected)"),
           p(verbatimTextOutput("verbatimOutput3"))),
    
    tabPanel(title = "Coverage", value = "coverageTab",
           textInput(inputId = "coverageTabInput", label = "Current Sample:"),
           h4(strong("Display Options:")),
           checkboxInput(inputId = "showVarTables", 
                         label = "Show Variant Tables",
                         value = FALSE),
           checkboxInput(inputId = "showCovPlots",
                         label = "Show Coverage Plots",
                         value = TRUE),
           uiOutput(outputId = "dynamicVarTables"),
           #hr(style = "border-top: 3px solid #2C3E50; margin-top: 100px; margin-bottom: 100px;"),
           #hr(style = "border: 0; height: 1px; background: #2C3E50; 
              #background-image: linear-gradient(to right, #ccc, #2C3E50, #ccc);
              #margin-top: 100px; margin-bottom: 100px;"),
           uiOutput(outputId = "dynamicCovPlots")),
    
    tabPanel(title = "Variants", value = "variantTab",
             textInput(inputId = "variantTabTextBox", label = "Current Sample:"),
             DT::DTOutput(outputId = "variantTabDT"),
             plotOutput(outputId = "variantTabCoveragePlot")),
    
    tabPanel(title = "Common Variants", value = "commonVariantsTab",
             textInput(inputId = "commonVariantsTabTextBox",
                       label = "Current Samples:"),
             DT::DTOutput(outputId = "commonVarTable"))
  )
)

server <- function(input, output, session) {
  
  hideTab(inputId = "navbarpage", target = "sampleTab")
  hideTab(inputId = "navbarpage", target = "coverageTab")
  hideTab(inputId = "navbarpage", target = "variantTab")
  hideTab(inputId = "navbarpage", target = "commonVariantsTab")
  
  # define the vector that will hold user-selected samples for coverage display:
  sampleVec <- vector(mode = "character") 
  
  #----- Go -----#
  
  # Get data set, hold in reactive variable so other functions may access
  GetSampleData <- eventReactive(input$go, {
    sampleData <- GenerateSampleData(dataSet = input$dataSetSelect)
    return(sampleData)
  })
  
  observeEvent(input$go, {
    showTab(inputId = "navbarpage", target = "sampleTab", select = TRUE)
    hideTab(inputId = "navbarpage", target = "coverageTab")
    hideTab(inputId = "navbarpage", target = "variantTab")
    hideTab(inputId = "navbarpage", target = "commonVariantsTab")
    
    # clear the vector that will hold user-selected samples for coverage disaplay:
    sampleVec <<- vector(mode = "character")
    
    #----- Data set name -----#
    output$current_dataset <- renderText({
      paste("Current data set: ", input$dataSetSelect)
    })

    #----- Display Sample Data Table -----#
  
    # columnDefs - hide the 3rd column (index 2) (number primary alignments)
    # formatStyle - add CSS style 'cursor: pointer' to the 1st column (i.e. sample)
    #   which is index 1 in in the DT object
    # Table set up to allow multiple cell selections:
    output$sampleDataTable <- DT::renderDataTable({
      data = datatable(GetSampleData(), 
                       selection = list(mode = "multiple", target = "cell"), 
                       rownames = FALSE,
                       options = list(pageLength = 10)) %>%
                         #bLengthChange = 0)) %>%
        formatStyle(columns = 1, cursor = "pointer") 
    })
  })
  
  #----- Download Sample Data Table -----#
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$dataSetSelect, "_sample_table.csv")
    },
    content = function(file) {
      write.csv(GetSampleData(), file, row.names = FALSE)
    }
  )

  #----- Sample Selection -----#
  
  #----- Select/Deselect All -----#
  dt_proxy <- DT::dataTableProxy(outputId = "sampleDataTable")
  
  observeEvent(input$dt_sel, {
    if (input$dt_sel == TRUE) {
      # select all cells in column 0
      DT::selectCells(dt_proxy, 
                      selected = as.matrix(cbind(1:nrow(GetSampleData()), 0)))
    } else {
      DT::selectCells(dt_proxy, selected = NULL)
    }
    output$selected_cells <- renderPrint(print(input$sampleDataTable_cells_selected))
  })

  # need to handle passing the sample vec differently if that checkbox is TRUE
  
  #----- Individual Sample Selection -----#
  
  # this needs to respond to the changing matrix instead of input$sampleDataTable_cells_selected
  GetSampleVec <- eventReactive(input$sampleDataTable_cells_selected, {
    # Returns character vector of Sample IDs.
    currentCell <- input$sampleDataTable_cell_clicked # row, col, value info
    currentCellCol <- currentCell$col
    if (currentCellCol == 0) { # value from Sample column must be selected
      currentSample <- currentCell$value
      if (!(currentSample %in% sampleVec)) {
        sampleVec <<- c(sampleVec, currentSample)
      } else { # remove sample
        sampleVec <<- sampleVec[sampleVec != currentSample]
      }
    }
    return(sampleVec)
  })
  
  # show value of sample vector in the sample table tab
  observeEvent(input$sampleDataTable_cells_selected, {
    output$verbatimOutput1 <- renderPrint(GetSampleVec())
    })
  
  # Calculate a matrix of correctly-selected cell indices.
  #   Used to determine when to enable/disable the action buttons by seeing if
  #   anything is selected on the sample data table.  I have to do this because
  #   I can't figure out how to prevent crashing from GetSampleVec()'s initial 
  #   empty/NULL value.  Solving that would simplify the code a little bit.
  SelectedSampleIndices <- eventReactive(input$sampleDataTable_cells_selected, {
    # User is free to select "incorrect" cells so they need to be filtered first.
    sampleMtxOriginal <- input$sampleDataTable_cells_selected
    sampleMtxFiltered <- NULL
    if (!(all(is.na(sampleMtxOriginal)))) {
      sampleMtxFiltered <- sampleMtxOriginal[sampleMtxOriginal[, 2] == 0, ,
                                             drop = FALSE]
    }
    return(sampleMtxFiltered)
  })
  
  # Enable/Disable Buttons & Hide/Show Tabs:
  observeEvent(input$sampleDataTable_cells_selected, {
    # Do nothing if nothing has been clicked, or the clicked cell isn't in the
    #   first column (which is index 0 for DT objects)
    if (is.null(SelectedSampleIndices()) | all(is.na(SelectedSampleIndices()))) {
      disable("stackCoverage")
      disable("inspVariants")
      disable("commonVariants")
      hideTab(inputId = "navbarpage", target = "coverageTab")
      hideTab(inputId = "navbarpage", target = "variantTab")
      hideTab(inputId = "navbarpage", target = "commonVariantsTab")
      return() 
    } else {
      numIndices <- nrow(SelectedSampleIndices())
    }
    
    if (numIndices == 1) { # single sample selected
      # enable inspVariants, disable stackCoverage and commonVariants
      # hide covTab if showing
      enable("inspVariants")
      disable("stackCoverage")
      hideTab(inputId = "navbarpage", target = "coverageTab")
      disable("commonVariants")
      hideTab(inputId = "navbarpage", target = "commonVariantsTab")
    } else if (numIndices >= 2) { # 2 + samples selected
      # enable stackCoverage and commonVariants, disable inspVariants
      # hide variantTab if showing
      enable("stackCoverage")
      enable("commonVariants")
      disable("inspVariants")
      hideTab(inputId = "navbarpage", target = "variantTab")
    }
  })
  
  observe({
    # what's selected?
    # indices of ALL selected cells
    output$verbatimOutput2 <- renderPrint(input$sampleDataTable_cells_selected)
    # indices of CORRECT selected cells
    output$verbatimOutput3 <- renderPrint(SelectedSampleIndices())
  })
  
  
  #----- Stack Coverage -----#  
  
  observeEvent(input$stackCoverage, {
    # make the variant tab visible & switch to it
    showTab(inputId = "navbarpage", target = "coverageTab", select = TRUE)
    # reset the checkboxes - don't show variant tables, do show cov plots
    updateCheckboxInput(session, inputId = "showVarTables", value = FALSE)
    updateCheckboxInput(session, inputId = "showCovPlots", value = TRUE)

    updateTextInput(session, inputId = "coverageTabInput", value = GetSampleVec())
    
    # for showing/hiding variant tables
    observeEvent(input$showVarTables, {
      if(input$showVarTables == TRUE) {
        shinyjs::show(id = "dynamicVarTables")
      } else {
        shinyjs::hide(id = "dynamicVarTables")
      }
    })
    
    # for showing/hiding coverage plots
    observeEvent(input$showCovPlots, {
      if(input$showCovPlots == TRUE) {
        shinyjs::show(id = "dynamicCovPlots")
      } else {
        shinyjs::hide(id = "dynamicCovPlots")
      }
    })
    
    # render the variant tables and coverage plots dynamically:
    output$dynamicVarTables <- renderUI({
      # Dynamic Variant Tables:
      allVariantData <- map(isolate(GetSampleVec()), GetVCF, 
                            dataSet = input$dataSetSelect)
      
      allVariantTables <- map(allVariantData, function(x) {
        p(h4(strong(comment(x)), style = "text-align: center;"),
          DT::renderDataTable(expr =
                              (data = datatable(x,
                                                #caption = htmltools::tags$caption(
                                                  #style = 'caption-side: top; text-align: center; color:black;',
                                                  #htmltools::strong(comment(x))),
                                                selection = list(mode = "multiple",
                                                                 target = "cell"),
                                                options = list(pageLength = 5),
                                                               #bLengthChange = 0,
                                                               #paging = FALSE,
                                                               #bFilter = 0),
                                                rownames = FALSE)) %>%
                              formatStyle(columns = 2, cursor = "pointer")))
        })
    })
    
    # Dynamic Coverage Plots:
    output$dynamicCovPlots <- renderUI({
        allCoveragData <- map(GetSampleVec(), function(y) {
          p(h4(strong(y), style = "text-align: center;"),
            renderPlot({PlotCoverage(dataSet = input$dataSetSelect, sample = y)},
                       height = 250), style = "margin-top: 30 px;")
        })
      })
    })
  
  #----- Find Common Variants -----#
  
  observeEvent(input$commonVariants, {
    showTab(inputId = "navbarpage", target = "commonVariantsTab", select = TRUE)
    updateTextInput(session, inputId = "commonVariantsTabTextBox",
                    value = GetSampleVec())
    
    # Put the selected samples and their variants into Sample objects
    currentSampleObjects <- BuildSampleObjects(dataSet = input$dataSetSelect,
                                               sampleVector = GetSampleVec())
    
    output$commonVarTable <- DT::renderDataTable({
      data = datatable(data = FindExactMatches(sampleObjectList = currentSampleObjects))
    })
    
  })
  
  #----- Inspect Variants -----#  
  
  observeEvent(input$inspVariants, {
    showTab(inputId = "navbarpage", target = "variantTab", select = TRUE)
    updateTextInput(session, inputId = "variantTabTextBox", 
                    value = GetSampleVec())
    
    # Interactive variant table
    output$variantTabDT <- DT::renderDataTable({
      data = datatable(GetVCF(dataSet = input$dataSetSelect,
                              sample = GetSampleVec()),
                       selection = list(mode = "multiple", target = "cell"),
                       options = list(pageLength = 5),
                       rownames = FALSE) %>%
        formatStyle(columns = 2, cursor = "pointer")
    })
    
    # First renderPlot (no highlight tracks)
    output$variantTabCoveragePlot <- renderPlot({
      PlotCoverage(dataSet = input$dataSetSelect, sample = GetSampleVec())
    })
    
    # Re-set global variable previousMtxSize (used to determine when to
    #   re-draw variantTabCoveragePlot)
    previousMtxSize <<- 0
  })
  
  #----- Variant Selection -----#
  
  observeEvent(input$variantTabDT_cells_selected, {
    # Logical to determine whether to re-paint a blank plot
    redrawBlank <- FALSE
    # What's selected on the variant table?
    originalMatrix <- input$variantTabDT_cells_selected
    if (!(all(is.na(originalMatrix)))) {
      # Select only values that are from the correct column in the variant 
      #   table ("positions" - col 1)
      filteredMatrix <- originalMatrix[originalMatrix[ , 2] == 1, , 
                                       drop = FALSE]
      # Is anything left in the filteredMatrix?
      if (!(all(is.na(filteredMatrix)))) {
        # If something is left in the filtered matrix, continue...
        # First adjust the index values in the returned matrix,
        # because the DT object is 0-indexed, but the data frame is 1-indexed
        # (add 1 to the column values in the filtered matrix)
        currentMtxSize <- nrow(filteredMatrix)
        # if the current size of filtered matrix is different than before, 
        #   we need to re-draw the plot with new variants
        if (previousMtxSize != currentMtxSize) {
          # Adjust the indices
          newMatrix <- cbind((filteredMatrix[ , 1]), filteredMatrix[ , 2] + 1)
          currentSampleValue <- input$sampleDataTable_cell_clicked$value
          # Pull the data for selected variants from the variant call file
          vcf <- as.data.frame(GetVCF(dataSet = input$dataSetSelect,
                                      sample = GetSampleVec()))[newMatrix]
          # Make the plot with variants identified
          output$variantTabCoveragePlot <- renderPlot({
            PlotCoverage(dataSet = input$dataSetSelect,
                         sample = GetSampleVec(), 
                         positions = as.numeric(vcf))
          })
        }
        previousMtxSize <<- as.numeric(currentMtxSize)
        
      } else if (previousMtxSize != 0) {
        redrawBlank <- TRUE
        previousMtxSize <<- 0
      }
    } else {
      redrawBlank <- TRUE
      previousMtxSize <<- 0
    }
    if (redrawBlank == TRUE) {
      output$variantTabCoveragePlot <- renderPlot({
        PlotCoverage(dataSet = input$dataSetSelect,
                     sample = GetSampleVec())
      })
    }
  }) # end of variant selection observeEvent
    
}

shinyApp(ui = ui, server = server)
