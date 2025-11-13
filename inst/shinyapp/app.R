library(magrittr)
library(SummarizedExperiment)
library(shiny)
library(shinyBS)
library(shinyWidgets)
library(plotly)
library(DT)
library(shinycssloaders)
library(MetAlyzer)
library(limma)
library(tidyverse)
library(htmlwidgets)
library(svglite)
library(writexl)
library(bslib)
library(viridis)
library(viridisLite)
library(gridExtra)

ui <- fluidPage(
  # Define notification styles for 'Process' and 'Revert/Default' buttons
  tags$head(
    tags$style(HTML("
      .shiny-notification {
        width: 300px;
        position: fixed;
        top: 15px;
        right: 15px;
        opacity: 1;
      }
    "))
  ),
  titlePanel('Biocrates Metabolomics Analysis'),
  tabsetPanel(
    type = 'tabs',
    tabPanel(
      'Overview',
      sidebarLayout(
        sidebarPanel(
          tags$h4('Data Uploading', style = 'color:steelblue;font-weight:bold'),
          uiOutput('updateFileInput'),
          div(textOutput('textFileInput'), style = 'color:IndianRed;font-weight:bold;font-size:110%'),
          checkboxInput('exampleFile',
                        HTML('Explore app with example dataset from <b>Biocrates</b>'), 
                        value = FALSE),
          # Show data processing options only after file is uploaded
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           tags$h4('Data Processing', style = 'color:steelblue;font-weight:bold'),
                           shinyBS::bsCollapse(
                             open = c('Sample filtering', 'Filtering Log'), multiple = T,
                             shinyBS::bsCollapsePanel('Sample filtering', style = 'info',
                                                      selectInput('smpChoicesFiltering', 'Select sample(s) to remove:',
                                                                  choices = character(0), multiple = T)),
                             shinyBS::bsCollapsePanel(
                               'Metabolite filtering', style = 'info',
                               selectInput('featChoicesFiltering', 'Select metabolite(s) to remove:',
                                           choices = character(0), multiple = T),
                               fluidRow(
                                 column(width = 8, sliderInput('featCompleteCutoffFiltering',
                                                               'Select % of quantified values per metabolite:',
                                                               min = 0, max = 100, value = 80))
                               ),
                               fluidRow(
                                 column(width = 8, sliderInput('featValidCutoffFiltering',
                                                               'Select % of valid quantifications per metabolite:',
                                                               min = 0, max = 100, value = 0)),
                                 column(width = 3, offset = 1,
                                        checkboxGroupInput('featValidStatusFiltering', 'Validity',
                                                           choices = c('Valid', 'LOQ', 'LOD', 'Invalid'),
                                                           selected = c('Valid', 'LOQ'))),
                                 shinyBS::bsTooltip('featCompleteCutoffFiltering',
                                                    'Metabolites with quantification rates below this cutoff are removed.'),
                                 shinyBS::bsTooltip('featValidCutoffFiltering',
                                                    'Metabolites with valid quantifications below this cutoff are removed.'),
                                 shinyBS::bsTooltip('featValidStatusFiltering',
                                                    'The selected is considered valid quantification for filtering.')
                               )
                             ),
                             shinyBS::bsCollapsePanel(
                               'Imputation and Normalization', style = 'info',
                               shinyWidgets::materialSwitch('imputation', 'Half-minimum (HM) imputation',
                                                            value = T, status = 'primary', right = T),
                               #### Needs further discussion on median normalization
                               selectInput('normalization', 'Select normalization method to use:',
                                           choices = c('None', 'Log2 transformation',
                                                       'Median normalization',
                                                       'Total ion count (TIC) normalization'),
                                           selected = 'Log2 transformation', multiple = F)
                               # shinyBS::bsTooltip('imputation', paste('Missing values are replaced with half of the minimum of,
                               #                                        observed values in each metabolite.')),
                               # shinyBS::bsTooltip('normalization', '')
                             )
                           ),
                           fluidRow(
                             column(width = 6, actionButton('updateProcessing', 'Process', width = '100%')),
                             column(width = 6, actionButton('revertProcessing', 'Revert/Default', width = '100%')),
                             shinyBS::bsTooltip('revertProcessing', 'The processed data and specified parameters revert to the origins.')
                           ),
                           #### Hide following functionalities until work is done
                           # tags$br(),
                           # tags$h4('Abundance Matrix Download', style = 'color:steelblue;font-weight:bold'),
                           # fluidRow(
                           #   style = "display: flex; align-items: flex-end;",
                           #   column(width = 7, selectInput('metaChoicesRawAbunExport', 'Select sample metadata to include:',
                           #                                 choices = character(0), multiple = T)),
                           #   column(width = 5, downloadButton('downloadRawAbunExport', 'Download',
                           #                                    style = 'width:100%; margin-bottom: 15px')),
                           #   shinyBS::bsTooltip('downloadRawAbunExport', 'This output can directly be used for MetaboAnalyst.')
                           # ),
                           tags$br(),
                           tags$br(),
                           tags$h4('Metabolite Identifier Table Download', style = 'color:steelblue;font-weight:bold'),
                           checkboxInput('ifUploadedFilePrior2023', 'Check box if uploaded file obtained prior 2023', value = F),
                           fluidRow(
                             style = "display: flex; align-items: flex-end;",
                             column(width = 7, uiOutput('loadFeatIdChoicesExport')),
                             column(width = 5, downloadButton('downloadFeatIdsExport', 'Download',
                                                              style = 'width:100%; margin-bottom: 15px'))
                             # shinyBS::bsTooltip('featIdChoicesExport', 'The identifiers were generated by MetaboAnalyst using HMDB IDs.',
                             #                    placement = 'top')
                           )
          )
        ),
        mainPanel(
          # JavaScript: Double equals (==) is loose equality with type coercion and
          # triple equals (===) is strict equality without type coercion
          conditionalPanel(condition = "output.ifValidUploadedFile === undefined",
                           HTML('<br>
                                <h3>This tool is here to help you analyze and explore your Biocrates data.</h3>
                                <br>
                                <p>Please ensure your Excel file contains the cell <b>"Class"</b> and the column 
                                <b>"Sample Type"</b> with correct indentation, as these are used as anchor cells!</p>
                                <br>
                                <p>Additional information is available by hovering over all action elements.</p>')
          ),
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           shinyBS::bsCollapse(
                             id = 'panelDatOverviewViz',
                             # Note that collapsed panel does not render output until it is expanded
                             open = c('Data distribution', 'Data completeness', 'Quantification status', 'Sample metadata (All)'),
                             multiple = T,
                             shinyBS::bsCollapsePanel('Sample metadata (All)', style = 'primary',
                                                      DT::dataTableOutput('tblSmpMetadat') %>%
                                                      #### Collapsed panel: Workaround (SHOULD WORK BUT NOT WORK)
                                                      # Render table outside collapsed panel and display it
                                                      # uiOutput('uiTblSmpMetadat') %>%
                                                        shinycssloaders::withSpinner(color="#56070C")),
                             shinyBS::bsCollapsePanel('Data distribution', style = 'primary',
                                                      fluidRow(
                                                        style = 'display:flex; align-items: center;',
                                                        column(width = 7, shinyBS::bsCollapse(
                                                          id = 'hintDatDist',
                                                          shinyBS::bsCollapsePanel('💡Hint', style = 'success',
                                                                                   textOutput('summDatDist', container = strong))
                                                        )),
                                                        column(width = 4, offset = 1, uiOutput('updateGpColsDatDist'))
                                                      ),
                                                      plotly::plotlyOutput('plotDatDist') %>%
                                                        shinycssloaders::withSpinner(color="#56070C"),
                                                      fluidRow(style="display:flex; justify-content:right; margin-top:1rem;",
                                                               column(width = 2,
                                                                      selectInput("formatDatDist", label = NULL,
                                                                                  choices = c("html", "png", "pdf", "svg"),
                                                                                  selected = "html")),
                                                               column(width = 2, downloadButton("downloadDatDist",
                                                                                                "Download")))
                                                      ),
                             shinyBS::bsCollapsePanel('Data completeness', style = 'primary',
                                                      fluidRow(
                                                        column(width = 7,
                                                               shinyBS::bsCollapse(
                                                                 id = 'hintDatComplete',
                                                                 shinyBS::bsCollapsePanel('💡Hint', style = 'success',
                                                                                          textOutput('summDatComplete', container = strong))
                                                               )
                                                        )
                                                      ),
                                                      plotly::plotlyOutput('plotDatComplete') %>%
                                                        shinycssloaders::withSpinner(color="#56070C"),
                                                      fluidRow(style="display:flex; justify-content:right; margin-top:1rem;",
                                                               column(width = 2,
                                                                      selectInput("formatDatComplete", label = NULL,
                                                                                  choices = c("html", "png", "pdf", "svg"),
                                                                                  selected = "html")),
                                                               column(width = 2, downloadButton("downloadDatComplete",
                                                                                                "Download")))
                                                      ),
                             shinyBS::bsCollapsePanel('Quantification status', style = 'primary',
                                                      #### Adjust tab size
                                                      fluidRow(
                                                        column(width = 7,
                                                               shinyBS::bsCollapse(
                                                                 id = 'hintQuanStatus',
                                                                 shinyBS::bsCollapsePanel('💡Hint', style = 'success',
                                                                                          textOutput('summQuanStatus', container = strong))
                                                               )
                                                        )
                                                      ),
                                                      plotly::plotlyOutput('plotQuanStatus') %>%
                                                        shinycssloaders::withSpinner(color="#56070C"),
                                                      fluidRow(style="display:flex; justify-content:right; margin-top:1rem;",
                                                               column(width = 2,
                                                                      selectInput("formatQuanStatus", label = NULL,
                                                                                  choices = c("html", "png", "pdf", "svg"),
                                                                                  selected = "html")),
                                                               column(width = 2, downloadButton("downloadQuanStatus",
                                                                                                "Download")))
                                                      )
                           )
          )
        )
      )
    ), # TabPanel 1 End
    tabPanel(
      'Analysis',
      sidebarLayout(
        conditionalPanel(condition = "output.ifValidUploadedFile",
                         sidebarPanel(
                           tags$h4('Differential Analysis', style = 'color:steelblue;font-weight:bold'),
                           fluidRow(
                             column(width = 5, selectInput('smpChoiceGpsLog2FC', 'Compare between:',
                                                           choices = 'Not available', multiple = F)),
                             column(width = 3, selectInput('smpChoicesLog2FC_1', 'Group1',
                                                           choices = character(0), multiple = F),
                                    offset = 1),
                             column(width = 3, selectInput('smpChoicesLog2FC_2', 'Group2',
                                                           choices = character(0), multiple = F))
                           ),
                           fluidRow(
                             style = "display: flex; align-items: center;",
                             column(width = 5, actionButton('computeLog2FC', 'Perform', width = '100%')),
                             column(width = 6, offset = 1,
                                    conditionalPanel(condition = "input.computeLog2FC > 0", # Check if the compute button has been clicked
                                                     downloadButton('downloadLog2FC', 'Download analysis result table', width = '100%')
                                    )
                             ),
                           ),
                           tags$br(),
                           tags$h4('Vulcano Plot', style = 'color:steelblue;font-weight:bold'), #Log₂(FC) Visualization
                           #### Highlighting in scatter plot is to be fixed 
                           fluidRow(
                             style = "display: flex; align-items: flex-end;",
                             column(width = 6, selectInput('metabChoicesVulcano',
                                                           'Select metabolite(s) to highlight:',
                                                           choices = character(0), multiple = T)),
                             column(width = 6, shinyWidgets::materialSwitch('highlightVulcano', 'Highlight',
                                                                            value = F, status = 'primary'))
                           ),
                           tags$h4('Select cutoffs:', style = 'font-weight:bold;font-size:14px'),
                           fluidRow(
                             column(width = 6, sliderInput('plotVolcanoLog2FCCutoff', 'Log₂(FC)',
                                                           min = 0, max = 10, value = 1, step = 0.1, ticks = F)),
                             column(width = 6, selectInput('plotVolcanoPValCutoff', 'q-value',
                                                           choices = c('0.0001', '0.001', '0.01', '0.05', '0.1'), #to avoid scientific notation
                                                           multiple = F, selected = 0.05))
                           )
                         )
        ),
        mainPanel(
          conditionalPanel(condition = "output.ifValidUploadedFile & input.computeLog2FC == 0",
                           div(textOutput('textLog2FC'), style = 'color:IndianRed;font-weight:bold;font-size:110%;margin-top:1rem;')),
          conditionalPanel(condition = "input.computeLog2FC",
                           tags$h3(strong('Vulcano Plot'), style = "margin-top:1rem; color: steelblue;"),
                           plotly::plotlyOutput('plotVolcano') %>%
                             shinycssloaders::withSpinner(color="#56070C"),
                           fluidRow(style="display:flex; justify-content:right; margin-top:1rem;",
                                    column(width = 2, 
                                           selectInput("formatVulcano", label = NULL, choices = c("html", "png", "pdf", "svg"), selected = "html")),
                                    column(width = 2, downloadButton("downloadVulcanoPlot", "Download vulcano plot"))
                           ),
                           tags$br(),
                           tags$h3(strong('Scatter Plot'), style = "margin-top:1rem; color: steelblue;"),
                           fluidRow(
                             column(width = 9, style = "z-index:2;", plotly::plotlyOutput('plotScatter') %>%
                                      shinycssloaders::withSpinner(color="#56070C")),
                             column(width = 3, style = "margin-left: -175px; z-index:1;",
                                    imageOutput('plotScatterLegend'))
                           ),
                           #### Downloaded scatter plot is to be fixed
                           fluidRow(style="display:flex; justify-content:right; margin-top:1rem; margin-bottom:1rem;",
                                    column(width = 2, 
                                           selectInput("formatScatter", label = NULL, choices = c("html", "png", "pdf", "svg"), selected = "html")),
                                    column(width = 2, downloadButton("downloadScatterPlot", "Download scatter plot"))
                           )
          )
        )
      )
    ), # TabPanel 2 End
    tabPanel(
      'Network',
      # Lower network plot a bit
      tags$br(),
      conditionalPanel(condition = "input.computeLog2FC",
                       shinyBS::bsCollapse(open = "",
                        shinyBS::bsCollapsePanel("Advanced Options",
                                                div(style = "display: grid;
                                                              grid-template-columns: repeat(5, 1fr);
                                                              gap: 15px;
                                                              align-items: center;
                                                              width: 80%;
                                                              margin: 0 auto 20px;",
                                                    # Plot Height Slider
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                        sliderInput("networkPlotHeight", 
                                                                    "Plot Height [100px]", 
                                                                    min = 4, max = 20, 
                                                                    value = 10, step = 1)
                                                    ),
                                                    # Metabolite Node Size Slider
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                        sliderInput("networkMetaboliteNodeSize", 
                                                                    "Metabolite Node Size", 
                                                                    min = 5, max = 50, 
                                                                    value = 11, step = 1)
                                                    ),
                                                    # Connection Width Slider
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                        sliderInput("networkConnectionWidth", 
                                                                    "Connection Width", 
                                                                    min = 0.5, max = 5, 
                                                                    value = 1.25, step = 0.25)
                                                    ),
                                                    # Pathway Text Size Slider
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                        sliderInput("networkPathwayTextSize", 
                                                                    "Pathway Text Size", 
                                                                    min = 10, max = 50, 
                                                                    value = 20, step = 1)
                                                    ),
                                                    # Pathway Width Slider
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                        sliderInput("networkPathwayWidth", 
                                                                    "Pathway Width", 
                                                                    min = 5, max = 30, 
                                                                    value = 10, step = 1)
                                                    ),
                                                    # Color Scale Selector --- ADDED
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                        selectInput("networkColorScale", "Color Scale",
                                                                    choices = c("Viridis", "Plasma", "Magma", "Inferno", 
                                                                                "Cividis", "Rocket", "Mako", "Turbo"),
                                                                    selected = "Viridis")
                                                    ),
                                                    # Exclude Pathways Selector --- ADDED
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                        selectInput("networkExcludePathways", "Exclude Pathways",
                                                                    choices = c("Bile Acids", "Eicosanoid Synthesis", "Hormones", 
                                                                                "Beta Oxidation", "Choline/Betaine metabolism", 
                                                                                "Lysine Metabolism", "Poly Amines", "Urea Cycle", 
                                                                                "TCA Cycle", "Glutamate Metabolism", 
                                                                                "Monoamine Metabolism", "Indole/Tryptophane Metabolism"),
                                                                    multiple = TRUE)
                                                    ),
                                                    # Column Name Selector --- ADDED
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                        selectInput("networkValueColumn", "Plotted Value",
                                                                    choices = c("log2FC", "pval", "qval", "tval"),
                                                                    selected = "log2FC")
                                                    ),
                                                    #Download controls --- MOVED
                                                    div(style = "min-width: 150px;",
                                                        selectInput("formatNetwork", label = "Format", choices = c("html", "png", "pdf", "svg"), selected = "html", width = '100%')
                                                    ),
                                                    div(style = "min-width: 150px;",
                                                        downloadButton("downloadNetworkPlot", "Download", style = "width: 100%; margin-top: 5px;")
                                                    ),
                                                    #Default Settings --- MOVED
                                                    div(style = "min-width: 150px; margin: 5px;"),
                                                    div(style = "min-width: 150px; margin: 5px;"),
                                                    div(style = "min-width: 150px; margin: 5px;",
                                                      actionButton('defaultNetworkPlotStyles', 'Default', width = '100%'),
                                                      shinyBS::bsTooltip('defaultNetworkPlotStyles', 'The changed plot style parameters revert to default.')
                                                    ),
                                                )
                        )
                       )
      ),
      conditionalPanel(condition = "output.ifValidUploadedFile & input.computeLog2FC == 0",
                       div(textOutput('textLog2FC_2'), style = 'color:IndianRed;font-weight:bold;font-size:110%;text-align:center')),
      conditionalPanel(condition = "input.computeLog2FC",
                       div(style = "width: 100%;",
                           # Use the height value from the slider to control the plot's height
                           plotly::plotlyOutput('plotNetwork', height = "auto") %>%
                             shinycssloaders::withSpinner(color="#56070C"),
                       ),
                       div(style = "width: 80%; margin: auto; margin-bottom: 200px;",
                           tags$h3(strong('Node Stats'), style = "margin-top:1rem;"),
                           selectInput("selectedNodesVulcano", "Select node(s) to view:",
                                       choices = character(0), multiple = TRUE),
                           plotly::plotlyOutput('plotVolcanoNodes') %>%
                             shinycssloaders::withSpinner(color="#56070C"),
                           downloadButton("downloadNodesExcel", "Download node stats as Excel"))
                           
      )
    ),
    tabPanel(
      'History',
      conditionalPanel(condition = "output.ifValidUploadedFile",
                       style = 'margin-left: 5mm;',
                       fluidRow(
                         column(width = 7, HTML('<br>
                                                <h3 style = "font-weight: bold; color: steelblue;">Processing History</h3>
                                                <h5>This section logs processing
                                                parameters used in different trials (every time you hit Process button).</h5>
                                                <h5><strong>Instructions:</strong> Click on a row to
                                                view the processed data at the specific point.</h5>
                                                <br>')),
                         column(width = 2, div(style = "width: 60%; margin-top: 30%;",
                                               actionButton('clearParamLog', 'Clear', width = '100%'))),
                         shinyBS::bsTooltip('clearParamLog', 'The parameter log table will be cleared.')
                       ),
                       fluidRow(
                         column(width = 7, DT::dataTableOutput("tblParamLog")),
                         column(width = 5,
                                conditionalPanel(condition = "typeof input.tblParamLog_rows_selected === 'undefined' || input.tblParamLog_rows_selected.length == 0",
                                                 HTML('<br>
                                                      <h4 style="color: IndianRed; text-align: left;">Please select a row.</h4>')
                                ),
                                conditionalPanel(condition = "typeof input.tblParamLog_rows_selected !== 'undefined' && input.tblParamLog_rows_selected.length > 0",
                                                 shinyBS::bsCollapse(
                                                   id = 'panelTrialLog',
                                                   open = c('Features Removed', 'Data distribution', 'Data completeness', 'Quantification status'),
                                                   multiple = T,
                                                   shinyBS::bsCollapsePanel('Features Removed', style = 'info',
                                                                            textOutput('textRmFeatsLog') %>%
                                                                              shinycssloaders::withSpinner(color="#56070C")),
                                                   shinyBS::bsCollapsePanel('Data distribution', style = 'info',
                                                                            plotly::plotlyOutput('plotDatDistLog') %>%
                                                                              shinycssloaders::withSpinner(color="#56070C")),
                                                   shinyBS::bsCollapsePanel('Data completeness', style = 'info',
                                                                            plotly::plotlyOutput('plotDatCompleteLog') %>%
                                                                              shinycssloaders::withSpinner(color="#56070C")),
                                                   shinyBS::bsCollapsePanel('Quantification status', style = 'info',
                                                                            plotly::plotlyOutput('plotQuanStatusLog') %>%
                                                                              shinycssloaders::withSpinner(color="#56070C"))
                                                 )
                                )
                         )
                       ),
                       tags$br(),
                       tags$br(),
                       fluidRow(
                         column(width = 7,
                                tags$h3(strong('Command History (Last Operation)'), style = 'color: steelblue;'),
                                verbatimTextOutput('textCommands'),
                                downloadButton('downloadSessionInfo', 'Download session info'),
                                shinyBS::bsTooltip('downloadSessionInfo',
                                                   'Critical for scientific reporting (e.g., package versions).'))
                       )
      )
    )
  )
)

server <- function(input, output, session) {
  # Create reactive objects for storing up-to-date data
  reactMetabObj <- reactiveValues(metabObj = NULL, oriMetabObj = NULL, tmpMetabObj = NULL)
  reactOriSmpMetadatTbl <- reactiveVal()
  reactLog2FCTbl <- reactiveVal()
  reactVulcanoHighlight <- reactiveVal()
  # Create reactive parameter list for detecting any parameter change by users,
  # which is to make data always follow built processing procedure
  reactParamList <- reactiveValues(smpFiltering = c(), featFiltering = c(), featCompleteCutoff = 0,
                                   featValidCutoff = 0, featValidStatus = c(), imputation = F,
                                   normalization = 'None')
  # Create reactive container to store processing parameters and processed data
  # of different trials for logging history and comparing results
  reactAnalysisLog <- reactiveValues(
    paramLogTbl = data.frame(idx = integer(),
                             smpFiltering = character(),
                             featClassFiltering = character(),
                             featFilteringNum = numeric(),
                             featCompleteCutoff = numeric(),
                             featValidCutoff = numeric(),
                             featValidStatus = character(),
                             imputation = character(),
                             normalization = character(),
                             stringsAsFactors = FALSE),
    rmFeatList = list(),
    # metabObjList = list(),
    # smpMetadatTblList = list(),
    metabAggreTblList = list()
  )
  # Create reactive object for recording R commands executed to show users
  reactCodeHistory <- reactiveVal()
  
  
  # Show fileInput only when example data is not used
  output$updateFileInput <- renderUI({
    if (!input$exampleFile) {
      fileInput('uploadedFile', NULL, multiple = F, accept = '.xlsx',
                placeholder = 'No .xlsx file selected')
    }
  })
  # Show hint for uploading data when example data is being used
  output$textFileInput <- renderText({
    if (input$exampleFile) {
      'Uncheck example data box to upload your own data.'
    }
  })
  # Bring UI back to origin when example data box is unchecked
  observe({
    req(!input$exampleFile)
    reactMetabObj$metabObj <- NULL
    
    # Uncheck box indicating if file was generated before 2023 when example file is dropped
    updateCheckboxInput(session, 'ifUploadedFilePrior2023', value = F)
  })
  # Initialize MetAlyzer SE object with example data
  observeEvent(input$exampleFile, {
    req(input$exampleFile)
    #### Collapsed panel: Workaround
    # Expand all panels to render outputs
    # shinyBS::updateCollapse(session, 'panelDatOverviewViz',
    #                         open = c('Data distribution', 'Data completeness', 'Quantification status', 'Sample metadata (All)'))
    
    metabObj <- MetAlyzer::read_webidq(file_path = MetAlyzer::load_demodata_biocrates(), silent = T)
    # Exclude 'Metabolism Indicators' from subsequent processing and analysis
    metabObj <- MetAlyzer::filter_metabolites(metabObj,
                                              drop_metabolites = 'Metabolism Indicators',
                                              drop_NA_concentration = F)
    reactMetabObj$metabObj <- metabObj
    # Make copy of original data for filtering reset
    reactMetabObj$oriMetabObj <- metabObj
    # Make copy of original data for sanity check of filtering
    reactMetabObj$tmpMetabObj <- metabObj
    # Monitor whether example data is used to show further operations on client side
    output$ifValidUploadedFile <- reactive({
      !is.null(reactMetabObj$metabObj)
    })
    outputOptions(output, 'ifValidUploadedFile', suspendWhenHidden = F)
    
    # Prepare static sample metadata table 
    smpMetadatTbl <- colData(reactMetabObj$oriMetabObj) %>%
      tibble::as_tibble(rownames = 'ID') %>%
      dplyr::mutate(ID = paste0('Smp', ID))
    # Use original column names whose spaces are not replaced with '.'
    colnames(smpMetadatTbl) <- c('ID', colnames(colData(reactMetabObj$oriMetabObj)))
    reactOriSmpMetadatTbl(smpMetadatTbl)
    
    # Set parameters back to default
    updateSliderInput(session, 'featCompleteCutoffFiltering', value = 80)
    updateSliderInput(session, 'featValidCutoffFiltering', value = 0)
    updateCheckboxGroupInput(session, 'featValidStatusFiltering', selected = c('Valid', 'LOQ'))
    updateMaterialSwitch(session, 'imputation', value = T)
    updateSelectInput(session, 'normalization', selected = 'Log2 transformation')
    doneImputation(0)
    doneNormalization(0)
    ifParamChange(0)
    doneSmpFiltering(0)
    doneFeatFiltering(0)
    
    updateSelectInput(session, 'featIdChoicesExport', selected = character(0))
    
    # Empty parameter log
    reactParamList$smpFiltering <- c()
    reactParamList$featFiltering <- c()
    reactParamList$featCompleteCutoff <- 0
    reactParamList$featValidCutoff <- 0
    reactParamList$featValidStatus <- c()
    reactParamList$imputation <- F
    reactParamList$normalization <- 'None'
    
    # Check box indicating if file was generated before 2023, because example file is
    updateCheckboxInput(session, 'ifUploadedFilePrior2023', value = T)
    
    # Record commands executed
    reactCodeHistory('######## Data Preparation ########')
    reactCodeHistory(c(reactCodeHistory(),
                       '# Initialize MetAlyzer SE object with example data'))
    reactCodeHistory(c(reactCodeHistory(),
                       'MetAlyzer::read_webidq(file_path = MetAlyzer::load_demodata_biocrates())'))
    reactCodeHistory(c(reactCodeHistory(),
                       '# Exclude "Metabolism Indicators"'))
    reactCodeHistory(c(reactCodeHistory(),
                       'metabObj <- MetAlyzer::filter_metabolites(metalyzer_se = metabObj, drop_metabolites = "Metabolism Indicators")'))
    reactCodeHistory(c(reactCodeHistory(), '\n'))
    reactCodeHistory(c(reactCodeHistory(),
                       '# Visualize data distribution, missing pattern, and quantification quality using plotly::ggplotly()'))
    reactCodeHistory(c(reactCodeHistory(),
                       '# Check source code about how data missingness and quantification status are summarized'))
  })
  # Initialize MetAlyzer SE object with uploaded data
  observeEvent(input$uploadedFile, {
    #### Collapsed panel: Workaround
    # Expand all panels to render outputs
    # shinyBS::updateCollapse(session, 'panelDatOverviewViz',
    #                         open = c('Data distribution', 'Data completeness', 'Quantification status', 'Sample metadata (All)'))
    
    validUploadedFile <- try(
      metabObj <- MetAlyzer::read_webidq(file_path = input$uploadedFile$datapath,
                                         sheet = 1, silent = T),
      silent = T)
    if (!is(validUploadedFile, 'try-error')) {
      # Exclude 'Metabolism Indicators' from subsequent processing and analysis
      if ('Metabolism Indicators' %in% unique(rowData(metabObj)$metabolic_classes)) {
        metabObj <- MetAlyzer::filter_metabolites(metabObj,
                                                  drop_metabolites = 'Metabolism Indicators',
                                                  drop_NA_concentration = F)
      }
      reactMetabObj$metabObj <- metabObj
      # Make copy of original data for filtering reset
      reactMetabObj$oriMetabObj <- metabObj
      # Make copy of original data for sanity check of filtering
      reactMetabObj$tmpMetabObj <- metabObj
      # Monitor whether uploaded file is valid to show further operations on client side
      output$ifValidUploadedFile <- reactive({
        !is.null(reactMetabObj$metabObj)
      })
      outputOptions(output, 'ifValidUploadedFile', suspendWhenHidden = F)
      
      # Prepare static sample metadata table 
      smpMetadatTbl <- colData(reactMetabObj$oriMetabObj) %>%
        tibble::as_tibble(rownames = 'ID') %>%
        dplyr::mutate(ID = paste0('Smp', ID))
      # Use original column names whose spaces are not replaced with '.'
      colnames(smpMetadatTbl) <- c('ID', colnames(colData(reactMetabObj$oriMetabObj)))
      reactOriSmpMetadatTbl(smpMetadatTbl)
      
      # Set parameters back to default
      updateSliderInput(session, 'featCompleteCutoffFiltering', value = 80)
      updateSliderInput(session, 'featValidCutoffFiltering', value = 0)
      updateCheckboxGroupInput(session, 'featValidStatusFiltering', selected = c('Valid', 'LOQ'))
      updateMaterialSwitch(session, 'imputation', value = T)
      updateSelectInput(session, 'normalization', selected = 'Log2 transformation')
      doneImputation(0)
      doneNormalization(0)
      ifParamChange(0)
      doneSmpFiltering(0)
      doneFeatFiltering(0)
      
      updateCheckboxInput(session, 'ifUploadedFilePrior2023', value = F)
      updateSelectInput(session, 'featIdChoicesExport', selected = character(0))
      
      # Empty parameter log
      reactParamList$smpFiltering <- c()
      reactParamList$featFiltering <- c()
      reactParamList$featCompleteCutoff <- 0
      reactParamList$featValidCutoff <- 0
      reactParamList$featValidStatus <- c()
      reactParamList$imputation <- F
      reactParamList$normalization <- 'None'
      
      # Record commands executed
      reactCodeHistory('######## Data Preparation ########')
      reactCodeHistory(c(reactCodeHistory(),
                         '# Initialize MetAlyzer SE object'))
      reactCodeHistory(c(reactCodeHistory(),
                         'metabObj <- MetAlyzer::read_webidq(file_path = "path_to_your_file")'))
      reactCodeHistory(c(reactCodeHistory(),
                         '# Exclude "Metabolism Indicators" if they exist'))
      reactCodeHistory(c(reactCodeHistory(),
                         'metabObj <- MetAlyzer::filter_metabolites(metalyzer_se = metabObj, drop_metabolites = "Metabolism Indicators")'))
      reactCodeHistory(c(reactCodeHistory(), '\n'))
      reactCodeHistory(c(reactCodeHistory(),
                         '# Visualize data distribution, missing pattern, and quantification quality using plotly::ggplotly()'))
      reactCodeHistory(c(reactCodeHistory(),
                         '# Check source code about how data missingness and quantification status are summarized'))
    } else {
      showModal(modalDialog(
        title = 'Uploaded file reading failed...',
        tags$strong('Please check if the uploaded file is exported from webidq Software.'),
        checkboxInput("check1", "Is the file exported from WebIDQ?", value = FALSE),
        checkboxInput("check2", "Does the file contain the 'Class' cell?", value = FALSE),
        checkboxInput("check3", "Does the file contain the 'Sample Type' column?", value = FALSE),
        checkboxInput("check4", "The Sample Type column contains the value 'Sample'? (only for rows with samples)", value = FALSE),
        checkboxInput("check5", "Are there NO duplicate Metabolite columns?", value = FALSE),
        easyClose = TRUE,
        footer = NULL
      ))
      reactMetabObj$metabObj <- NULL
    }
  })
  #### Collapsed panel: Workaround
  # Close expanded panels (due to rendering) right after output is rendered
  observeEvent(reactOriSmpMetadatTbl(), {
    shinyBS::updateCollapse(session, 'panelDatOverviewViz', close = 'Sample metadata (All)')
  })
  
  # Retrieve abundance data and sample metadata and compute feature completeness
  # level and quantification status validity level for showing data overviews
  datOverviewPack <- reactive({
    req(reactMetabObj$metabObj)
    metabAggreTbl <- metadata(reactMetabObj$metabObj)$aggregated_data %>%
      #### Do not know why Metabolite column is grouped by MetAlyzer
      dplyr::ungroup() %>%
      dplyr::mutate(ID = paste0('Smp', ID))
    smpMetadatTbl <- colData(reactMetabObj$metabObj) %>%
      tibble::as_tibble(rownames = 'ID') %>%
      dplyr::mutate(ID = paste0('Smp', ID))
    # Use original column names whose spaces are not replaced with '.'
    colnames(smpMetadatTbl) <- c('ID', colnames(colData(reactMetabObj$metabObj)))
    # Prepare ID levels for displaying samples in order
    idLevels <- rownames(colData(reactMetabObj$metabObj))
    metabAggreTbl <- dplyr::left_join(metabAggreTbl, smpMetadatTbl, by = 'ID') %>%
      dplyr::mutate(ID = factor(ID, levels = paste0('Smp', idLevels)))
    
    # Compute completeness level of each feature for doing filtering
    featCompleteLvTbl <- dplyr::select(metabAggreTbl, Metabolite, Concentration) %>%
      #### 0 is true zero in biocrates data
      dplyr::mutate(Obs = dplyr::case_when(!Concentration %in% NA ~ 1), #c(0, NA)
                    Miss = dplyr::case_when(Concentration %in% NA ~ 1)) %>% #c(0, NA)
      dplyr::group_by(Metabolite) %>%
      dplyr::summarise(ObsCount = sum(Obs, na.rm = T), MissCount = sum(Miss, na.rm = T)) %>%
      dplyr::mutate(TotalCount = ncol(reactMetabObj$metabObj),
                    CompleteRatio = ObsCount / TotalCount)
    # Compute status counts of each feature for doing filtering
    featStatusCountTbl <- dplyr::select(metabAggreTbl, Metabolite, Status) %>%
      dplyr::mutate(Valid = dplyr::case_when(Status %in% 'Valid' ~ 1),
                    LOQ = dplyr::case_when(Status %in% 'LOQ' ~ 1),
                    LOD = dplyr::case_when(Status %in% 'LOD' ~ 1),
                    Invalid = dplyr::case_when(Status %in% 'Invalid' ~ 1)) %>%
      dplyr::group_by(Metabolite) %>%
      dplyr::summarise(ValidCount = sum(Valid, na.rm = T), LOQCount = sum(LOQ, na.rm = T),
                       LODCount = sum(LOD, na.rm = T), InvalidCount = sum(Invalid, na.rm = T)) %>%
      dplyr::mutate(TotalCount = ncol(reactMetabObj$metabObj))
    
    return(list(metabAggreTbl = metabAggreTbl, smpMetadatTbl = smpMetadatTbl,
                featCompleteLvTbl = featCompleteLvTbl, featStatusCountTbl = featStatusCountTbl))
  })
  
  
  # Prepare choices for feature filtering and log2(FC) vulcano highlighting
  featChoices <- reactive({
    req(reactMetabObj$metabObj)
    list(Class = as.list(unique(rowData(reactMetabObj$tmpMetabObj)$metabolic_classes)),
         Metabolite = as.list(rownames(reactMetabObj$tmpMetabObj)))
  })
  # Prepare needed information for sample filtering (i.e., choices, dictionary for
  # choice group searches, and list of choice groups containing duplicated choices)
  # and log2(FC) calculation
  smpChoicePack <- reactive({
    req(reactMetabObj$metabObj)
    smpChoiceList <- as.list(colData(reactMetabObj$tmpMetabObj))
    # Uniquify choices in every choice group and convert NA into character so
    # that it can be shown on client side
    for (i in seq_along(smpChoiceList)) {
      choices <- unique(smpChoiceList[[i]])
      if (any(is.na(choices))) {
        choices <- choices[-which(is.na(choices))]
        smpChoiceList[[i]] <- c(choices, 'NA')
      } else {
        smpChoiceList[[i]] <- choices
      }
    }
    # Remove choice groups with only one level
    # rmChoiceGps <- which(sapply(smpChoiceList, length) == 1)
    # if (length(rmChoiceGps) != 0) {
    #   smpChoiceList <- smpChoiceList[-rmChoiceGps]
    # }
    # Make identical choices in different choice groups unique by adding suffixes
    allSmpChoices <- unlist(smpChoiceList)
    dupSmpChoices <- allSmpChoices[duplicated(allSmpChoices)]
    # Collect choice groups with duplicated choices for later choice group searches
    # (due to design of MetAlyzer::filter_meta_data)
    dupSmpChoiceGps <- c()
    if (length(dupSmpChoices) != 0) {
      for (i in seq_along(smpChoiceList)) {
        choices <- smpChoiceList[[i]]
        if (any(dupSmpChoices %in% choices)) {
          dupIdx <- which(choices %in% dupSmpChoices)
          choices[dupIdx] <- paste0(choices[dupIdx], ' (', names(smpChoiceList)[i], ')')
          smpChoiceList[[i]] <- choices
          dupSmpChoiceGps <- c(dupSmpChoiceGps, names(smpChoiceList)[i])
        }
      }
    }
    
    # Create dictionary for later choice group searches of certain choices (due
    # to design of MetAlyzer::filter_meta_data)
    choiceGpSizes <- sapply(smpChoiceList, length)
    smpChoices2Gps <- lapply(seq_along(choiceGpSizes), function(i) {
      rep(names(choiceGpSizes)[i], choiceGpSizes[i])
    })
    smpChoices2Gps <- unlist(smpChoices2Gps)
    names(smpChoices2Gps) <- unlist(smpChoiceList)
    
    return(list(smpChoiceList = smpChoiceList, smpChoices2Gps = smpChoices2Gps,
                dupSmpChoiceGps = dupSmpChoiceGps))
  })
  # Update choices for feature filtering
  observe({
    req(featChoices())
    if (nrow(reactMetabObj$metabObj) == nrow(reactMetabObj$oriMetabObj)) {
      #### 'Metabolism Indicators' are excluded from subsequent processing and analysis for now
      # if ('Metabolism Indicators' %in% unlist(featChoices())) {
      #   updateSelectInput(session, 'featChoicesFiltering', choices = featChoices(),
      #                     selected = 'Metabolism Indicators')
      # } else {
      updateSelectInput(session, 'featChoicesFiltering', choices = featChoices())
      # }
    }
  })
  # Update choices for sample filtering
  observe({
    req(smpChoicePack()$smpChoiceList)
    # Make choices lists so that sole choice in certain choice group can be shown
    smpChoiceList <- lapply(smpChoicePack()$smpChoiceList, function(choices) {
      as.list(choices)
    })
    if (ncol(reactMetabObj$metabObj) == ncol(reactMetabObj$oriMetabObj)) {
      updateSelectInput(session, 'smpChoicesFiltering', choices = smpChoiceList)
    }
  })
  
  # Create reactive values to monitor if data imputation and normalization are done
  # to avoid data normalized more than one time
  doneImputation <- reactiveVal(0)
  doneNormalization <- reactiveVal(0)
  # Create reactive values to monitor if any parameter is changed
  ifParamChange <- reactiveVal(0)
  # Create reactive values to monitor if filtering is conducted for recording commands executed
  doneSmpFiltering <- reactiveVal(0)
  doneFeatFiltering <- reactiveVal(0)
  
  # Rerun data processing if any parameter is changed, so that processing can follow
  # order from sample filtering, feature filtering, imputation, to normalization
  # Note that every first trial must be logged because output of selectInput will
  # never be c(), instead, NULL when nothing is selected
  observeEvent(input$updateProcessing, {
    if (any(!identical(sort(reactParamList$smpFiltering), sort(input$smpChoicesFiltering)),
            !identical(sort(reactParamList$featFiltering), sort(input$featChoicesFiltering)),
            !identical(reactParamList$featCompleteCutoff, input$featCompleteCutoffFiltering),
            !identical(reactParamList$featValidCutoff, input$featValidCutoffFiltering),
            !identical(reactParamList$featValidStatus, input$featValidStatusFiltering),
            !identical(reactParamList$imputation, input$imputation),
            !identical(reactParamList$normalization, input$normalization))) {
      # Revert processed temporary MetAlyzer object and parameter log back to origins
      reactMetabObj$tmpMetabObj <- reactMetabObj$oriMetabObj
      doneImputation(0)
      doneNormalization(0)
      doneSmpFiltering(0)
      doneFeatFiltering(0)
      # Remove record of 'Data Preprocessing' part if preprocessing is re-conducted
      reactCodeHistory(reactCodeHistory()[seq(8)])
      
      reactParamList$smpFiltering <- c()
      reactParamList$featFiltering <- c()
      reactParamList$featCompleteCutoff <- 0
      reactParamList$featValidCutoff <- 0
      reactParamList$featValidStatus <- c()
      reactParamList$imputation <- F
      reactParamList$normalization <- 'None'
      
      # Turn reactive ifParamChange to 1, so used parameters and processed data
      # will be logged
      ifParamChange(1)
    }
  })
  # Do sample filtering (place this chunk before 'Do feature filtering', so that
  # feature filtering is executed based on sample filtered data)
  observeEvent(input$updateProcessing, {
    req(input$smpChoicesFiltering)
    for (selectedChoice in input$smpChoicesFiltering) {
      selectedChoiceGp <- smpChoicePack()$smpChoices2Gps[selectedChoice]
      # Prepare up-to-date choices of selected choice group to avoid error that
      # metadata of same sample is selected multiple times at once
      choices <- unique(colData(reactMetabObj$tmpMetabObj)[[selectedChoiceGp]])
      # Convert NA into character so that it can be searched
      if (any(is.na(choices))) {
        choices <- choices[-which(is.na(choices))]
        choices <- c(choices, 'NA')
      }
      if (!selectedChoiceGp %in% smpChoicePack()$dupSmpChoiceGps) {
        if (selectedChoice %in% choices) {
          if (selectedChoice == 'NA') {
            reactMetabObj$tmpMetabObj <- MetAlyzer::filter_meta_data(reactMetabObj$tmpMetabObj,
                                                                     !is.na(.data[[selectedChoiceGp]]))
          } else {
            reactMetabObj$tmpMetabObj <- MetAlyzer::filter_meta_data(reactMetabObj$tmpMetabObj,
                                                                     !.data[[selectedChoiceGp]] %in% selectedChoice)
          }
        }
      } else {
        # Prepare manually added suffix if selected choice group is in collected
        # choice groups with duplicates (see smpChoicePack() section)
        regexSuffix <- paste0(' \\(', selectedChoiceGp, '\\)')
        modSelectedChoice <- stringr::str_remove(selectedChoice, regexSuffix)
        if (modSelectedChoice %in% choices) {
          if (modSelectedChoice == 'NA') {
            reactMetabObj$tmpMetabObj <- MetAlyzer::filter_meta_data(reactMetabObj$tmpMetabObj,
                                                                     !is.na(.data[[selectedChoiceGp]]))
          } else {
            reactMetabObj$tmpMetabObj <- MetAlyzer::filter_meta_data(reactMetabObj$tmpMetabObj,
                                                                     !.data[[selectedChoiceGp]] %in% modSelectedChoice)
          }
        }
      }
      doneSmpFiltering(1) #for command history
    }
    # Avoid app crash when no sample is left
    if (ncol(reactMetabObj$tmpMetabObj) == 0) {
      showModal(modalDialog(
        title = 'No sample is left after sample filtering...',
        'Please respecify parameters for sample filtering.',
        easyClose = T,
        footer = NULL
      ))
    }
  })
  # Do feature filtering
  #### Add advanced feature filtering: Modified 80% rule
  observeEvent(input$updateProcessing, {
    req(datOverviewPack())
    # Collect features to remove based on user's selection
    rmSelectedFeats <- input$featChoicesFiltering
    # Collect features to remove based on missingness
    featCompleteLvTbl <- datOverviewPack()$featCompleteLvTbl
    featCompleteLevels <- featCompleteLvTbl$CompleteRatio
    featCompleteCutoff <- input$featCompleteCutoffFiltering / 100
    rmMissFeats <- featCompleteLvTbl$Metabolite[featCompleteLevels < featCompleteCutoff]
    # Collect features to remove based on validity
    # featStatusCountTbl <- datOverviewPack()$featStatusCountTbl
    featValidCutoff <- input$featValidCutoffFiltering / 100
    featValidStatus <- input$featValidStatusFiltering
    # featValidCounts <- rep(0, nrow(featStatusCountTbl))
    # for (status in featValidStatus) {
    #   featValidCounts <- featValidCounts + featStatusCountTbl[[paste0(status, 'Count')]]
    # }
    # featValidLevels <- featValidCounts / featStatusCountTbl[['TotalCount']]
    # rmInvalidFeats <- featStatusCountTbl$Metabolite[featValidLevels < featValidCutoff]
    # # Summarize features to remove
    # rmFeats <- unique(c(rmSelectedFeats, rmMissFeats)) # rmInvalidFeats
    # Skip feature filtering if no sample was left in sample filtering part since
    # feature filtering is based on sample filtered data
    if (ncol(reactMetabObj$tmpMetabObj) != 0) {
      #### Filtering results depend on whether rmSelectedFeats is Null when specifying
      #### rmFeats to 'drop_metabolites'????
      reactMetabObj$tmpMetabObj <- MetAlyzer::filter_metabolites(reactMetabObj$tmpMetabObj,
                                                                 drop_metabolites = rmSelectedFeats,
                                                                 drop_NA_concentration = FALSE)
      reactMetabObj$tmpMetabObj <- MetAlyzer::filter_metabolites(reactMetabObj$tmpMetabObj,
                                                                 drop_metabolites = rmMissFeats,
                                                                 drop_NA_concentration = FALSE,
                                                                 min_percent_valid = featValidCutoff,
                                                                 valid_status = featValidStatus,
                                                                 per_group = NULL)
      if (any(!is.null(rmSelectedFeats), featCompleteCutoff != 0,
              featValidCutoff != 0 & !all(c('Valid', 'LOQ', 'LOD', 'Invalid') %in% featValidStatus))) {
        doneFeatFiltering(1) #for command history
      }
      
      # Avoid app crash when no feature is left, e.g., min_percent_valid > 0 and valid_status == c()
      if (nrow(reactMetabObj$tmpMetabObj) == 0) {
        showModal(modalDialog(
          title = 'No metabolite is left after feature filtering...',
          'Please respecify parameters for metabolite filtering.',
          easyClose = T,
          footer = NULL
        ))
      }
    }
  })
  # Do data imputation and normalization
  observeEvent(input$updateProcessing, {
    req(reactMetabObj$tmpMetabObj)
    # Skip imputation and normalization if no sample or feature was left after filtering
    if (all(ncol(reactMetabObj$tmpMetabObj) != 0, nrow(reactMetabObj$tmpMetabObj) != 0)) {
      if (all(input$imputation, doneImputation() == 0)) {
        reactMetabObj$tmpMetabObj <- MetAlyzer:::data_imputation(reactMetabObj$tmpMetabObj)
        doneImputation(1)
      }
      if (doneNormalization() == 0) {
        if (input$normalization == 'Total ion count (TIC) normalization') {
          reactMetabObj$tmpMetabObj <- MetAlyzer:::data_normalization(reactMetabObj$tmpMetabObj,
                                                                      norm_method = 'TIC')
          doneNormalization(1)
        } else if (input$normalization == 'Median normalization') {
          reactMetabObj$tmpMetabObj <- MetAlyzer:::data_normalization(reactMetabObj$tmpMetabObj,
                                                                      norm_method = 'median')
          doneNormalization(1)
        } else if (input$normalization == 'Log2 transformation') {
          reactMetabObj$tmpMetabObj <- MetAlyzer:::data_normalization(reactMetabObj$tmpMetabObj,
                                                                      norm_method = 'log2')
          doneNormalization(1)
        }
      }
      # Update main MetAlyzer object for further analysis
      reactMetabObj$metabObj <- reactMetabObj$tmpMetabObj
      
      
      # Log parameters and processed data
      if (ifParamChange() == 1) {
        # Log parameters
        reactParamList$smpFiltering <- input$smpChoicesFiltering
        reactParamList$featFiltering <- input$featChoicesFiltering
        reactParamList$featCompleteCutoff <- input$featCompleteCutoffFiltering
        reactParamList$featValidCutoff <- input$featValidCutoffFiltering
        reactParamList$featValidStatus <- input$featValidStatusFiltering
        reactParamList$imputation <- input$imputation
        reactParamList$normalization <- input$normalization
        
        
        # Prepare removed feature list
        selectedFeats <- input$featChoicesFiltering
        featClasses <- unique(rowData(reactMetabObj$oriMetabObj)$metabolic_classes)
        rmFeatClasses <- selectedFeats[selectedFeats %in% featClasses]
        # Prepare removed feature class list
        oriFeatSpace <- rownames(reactMetabObj$oriMetabObj)
        filtFeatSpace <- rownames(reactMetabObj$metabObj)
        rmFeats <- oriFeatSpace[!oriFeatSpace %in% filtFeatSpace]
        
        # Prevent empty cells for better readability
        selectedSmps <- input$smpChoicesFiltering
        if (length(selectedSmps) == 0) {
          selectedSmps <- 'None'
        }
        
        if (length(rmFeatClasses) == 0) {
          rmFeatClasses <- 'None'
        }
        
        selectedValidStatus <- input$featValidStatusFiltering
        if (length(selectedValidStatus) == 0) {
          selectedValidStatus <- 'None'
        }
        
        newParamLog <- data.frame(
          idx = nrow(reactAnalysisLog$paramLogTbl) + 1,
          #### What is effect of I(list())?
          smpFiltering = I(list(paste(selectedSmps, collapse = ' | '))),
          featClassFiltering = paste(rmFeatClasses, collapse = ' | '),
          featFilteringNum = length(rmFeats),
          featCompleteCutoff = reactParamList$featCompleteCutoff,
          featValidCutoff = reactParamList$featValidCutoff,
          featValidStatus = I(list(paste(selectedValidStatus, collapse = ' | '))),
          imputation = as.character(reactParamList$imputation),
          normalization = reactParamList$normalization,
          stringsAsFactors = FALSE
        )
        reactAnalysisLog$paramLogTbl <- rbind(reactAnalysisLog$paramLogTbl, newParamLog)
        
        # Log processed data
        logIdx <- nrow(reactAnalysisLog$paramLogTbl)
        reactAnalysisLog$rmFeatList[[logIdx]] <- rmFeats
        # reactAnalysisLog$metabObjList[[logIdx]] <- reactMetabObj$metabObj
        # reactAnalysisLog$smpMetadatTblList[[logIdx]] <- datOverviewPack()$smpMetadatTbl
        reactAnalysisLog$metabAggreTblList[[logIdx]] <- datOverviewPack()$metabAggreTbl
        
        # Return reactive ifParamChange to 0, so this chunk will run only if any
        # parameter is changed
        ifParamChange(0)
        
        # Record commands executed
        if (any(doneSmpFiltering() != 0, doneFeatFiltering() != 0,
                doneImputation() != 0, doneNormalization() != 0)) {
          # Remove 'Differential Analysis' part if it is already recorded
          if ('######## Differential Analysis ########' %in% reactCodeHistory()) {
            reactCodeHistory(head(reactCodeHistory(), -10))
          }
          reactCodeHistory(c(reactCodeHistory(), '\n'))
          reactCodeHistory(c(reactCodeHistory(), '######## Data Preprocessing ########'))
          if (doneSmpFiltering() != 0) {
            reactCodeHistory(c(reactCodeHistory(),
                               '# Filter samples based on selected samples/sample groups'))
            reactCodeHistory(c(reactCodeHistory(),
                               'metabObj <- MetAlyzer::filter_meta_data(metalyzer_se = metabObj, ...)'))
          }
          if (doneFeatFiltering() != 0) {
            reactCodeHistory(c(reactCodeHistory(),
                               '# Filter metabolites based on selected metabolites/metabolic classes or specified parameters'))
            reactCodeHistory(c(reactCodeHistory(),
                               'metabObj <- MetAlyzer::filter_metabolites(metalyzer_se = metabObj, ...)'))
          }
          if (doneImputation() != 0) {
            reactCodeHistory(c(reactCodeHistory(),
                               '# Impute data using half-minimum'))
            reactCodeHistory(c(reactCodeHistory(),
                               'metabObj <- MetAlyzer:::data_imputation(metalyzer_se = metabObj)'))
          }
          if (doneNormalization() != 0) {
            reactCodeHistory(c(reactCodeHistory(),
                               '# Normalize data using selected method ("log2", "median", "TIC")'))
            reactCodeHistory(c(reactCodeHistory(),
                               'metabObj <- MetAlyzer:::data_normalization(metalyzer_se = metabObj, norm_method = "selected_method")'))
          }
        }
      } 
    } else {
      # Revert temporary MetAlyzer object to origin for redoing filtering
      reactMetabObj$tmpMetabObj <- reactMetabObj$oriMetabObj
    }
  })
  # Revert processed data and specified parameters back to origins
  observeEvent(input$revertProcessing, {
    req(reactMetabObj$metabObj)
    reactMetabObj$metabObj <- reactMetabObj$oriMetabObj
    reactMetabObj$tmpMetabObj <- reactMetabObj$oriMetabObj
    doneImputation(0)
    doneNormalization(0)
    ifParamChange(0)
    doneSmpFiltering(0)
    doneFeatFiltering(0)
    # Remove record of 'Data Preprocessing' part if preprocessing is re-conducted
    reactCodeHistory(reactCodeHistory()[seq(8)])
    
    # Set parameters for feature filtering back to default
    # if ('Metabolism Indicators' %in% unlist(featChoices())) {
    #   updateSelectInput(session, 'featChoicesFiltering', choices = featChoices(),
    #                     selected = 'Metabolism Indicators')
    # } else {
    updateSelectInput(session, 'featChoicesFiltering', choices = featChoices())
    # }
    updateSliderInput(session, 'featCompleteCutoffFiltering', value = 80)
    updateSliderInput(session, 'featValidCutoffFiltering', value = 0)
    updateCheckboxGroupInput(session, 'featValidStatusFiltering', selected = c('Valid', 'LOQ'))
    # Set parameters for sample filtering back to default
    # Make choices lists so that sole choice in certain choice group can be shown
    smpChoiceList <- lapply(smpChoicePack()$smpChoiceList, function(choices) {
      as.list(choices)
    })
    updateSelectInput(session, 'smpChoicesFiltering', choices = smpChoiceList)
    # Set parameters for imputation and normalization back to default
    updateMaterialSwitch(session, 'imputation', value = T)
    updateSelectInput(session, 'normalization', selected = 'Log2 transformation')
    
    # Empty parameter log
    reactParamList$smpFiltering <- c()
    reactParamList$featFiltering <- c()
    reactParamList$featCompleteCutoff <- 0
    reactParamList$featValidCutoff <- 0
    reactParamList$featValidStatus <- c()
    reactParamList$imputation <- F
    reactParamList$normalization <- 'None'
  })
  # Clear parameter log table
  observeEvent(input$clearParamLog, {
    reactAnalysisLog$paramLogTbl <- data.frame(idx = integer(),
                                               smpFiltering = character(),
                                               featClassFiltering = character(),
                                               featFilteringNum = numeric(),
                                               featCompleteCutoff = numeric(),
                                               featValidCutoff = numeric(),
                                               featValidStatus = character(),
                                               imputation = character(),
                                               normalization = character(),
                                               stringsAsFactors = FALSE)
    reactAnalysisLog$rmFeatList <- list()
    # reactAnalysisLog$metabObjList = list()
    # reactAnalysisLog$smpMetadatTblList = list()
    reactAnalysisLog$metabAggreTblList = list()
    # Turn reactive ifParamChange to 1, so current unchanged parameters can be logged
    # if Process button is hit
    ifParamChange(1)
  })
  
  
  # Update sample choice groups for log2(FC) calculation
  observe({
    req(smpChoicePack()$smpChoiceList)
    smpChoiceList <- smpChoicePack()$smpChoiceList
    # Remove choice groups whose level sizes are one or equal to sample size
    choiceGpSizes <- sapply(smpChoiceList, length)
    rmChoiceGps <- which(choiceGpSizes == 1 | choiceGpSizes == ncol(reactMetabObj$metabObj))
    if (length(rmChoiceGps) != 0) {
      smpChoiceList <- smpChoiceList[-rmChoiceGps]
    }
    smpChoiceGps <- names(smpChoiceList)
    if (length(smpChoiceGps) != 0) {
      updateSelectInput(session, 'smpChoiceGpsLog2FC', choices = smpChoiceGps)
    } else {
      updateSelectInput(session, 'smpChoiceGpsLog2FC', choices = 'Not available')
    }
  })
  # Update sample choices according to selected choice group for log2(FC) calculation
  observe({
    req(smpChoicePack()$smpChoiceList)
    selectedChoiceGp <- input$smpChoiceGpsLog2FC
    smpChoices <- smpChoicePack()$smpChoiceList[[selectedChoiceGp]]
    if (selectedChoiceGp %in% smpChoicePack()$dupSmpChoiceGps) {
      regexSuffix <- paste0(' \\(', selectedChoiceGp, '\\)')
      smpChoices <- stringr::str_remove(smpChoices, regexSuffix)
    }
    # Turn NULL (when selectedChoiceGp is 'Not available') into empty vector, or
    # choices will not be updated, i.e., previous choices are shown
    if (length(smpChoices) == 0) {
      smpChoices <- character(0)
    }
    updateSelectInput(session, 'smpChoicesLog2FC_1', choices = smpChoices)
    updateSelectInput(session, 'smpChoicesLog2FC_2', choices = smpChoices, selected = smpChoices[2])
  })
  # Compute log2(FC)
  observeEvent(input$computeLog2FC, {
    if (input$smpChoicesLog2FC_1 != input$smpChoicesLog2FC_2) {
      metabObj <- reactMetabObj$metabObj
      # Do log2 transformation on data for 'calc_log2FC' if normalization was not performed
      if (doneNormalization() == 0) {
        oriConc <- metadata(metabObj)$aggregated_data$Concentration
        metadata(metabObj)$aggregated_data$Concentration <- MetAlyzer:::glog2(oriConc)
      }
      # Extract samples of interest
      selectedChoiceGp <- input$smpChoiceGpsLog2FC
      selectedChoices <- c(input$smpChoicesLog2FC_1, input$smpChoicesLog2FC_2)
      if (!'NA' %in% selectedChoices) {
        metabObj <- MetAlyzer::filter_meta_data(metabObj, .data[[selectedChoiceGp]] %in% selectedChoices)
      } else {
        metabObj <- MetAlyzer::filter_meta_data(metabObj, is.na(.data[[selectedChoiceGp]]) |
                                                .data[[selectedChoiceGp]] %in% selectedChoices)
      }
      metabObj <- MetAlyzer:::calc_log2FC(metalyzer_se = metabObj,
                                          group = selectedChoiceGp,
                                          group_level = selectedChoices)
      reactLog2FCTbl(MetAlyzer:::log2FC(metabObj))
      
      # Record commands executed
      # Avoid recording if differential analysis is re-conducted
      if (!'######## Differential Analysis ########' %in% reactCodeHistory()) {
        reactCodeHistory(c(reactCodeHistory(), '\n'))
        reactCodeHistory(c(reactCodeHistory(), '######## Differential Analysis ########'))
        reactCodeHistory(c(reactCodeHistory(),
                           '# Compare selected Group1 with Group2 stored in sample metadata variable'))
        reactCodeHistory(c(reactCodeHistory(),
                           'metabObj <- MetAlyzer:::calc_log2FC(metalyzer_se = metabObj, group = "sample_metadata_variable")'))
        reactCodeHistory(c(reactCodeHistory(),
                           '# Make interactive vulcano plot with selected cutoffs'))
        reactCodeHistory(c(reactCodeHistory(),
                           'MetAlyzer:::plotly_vulcano(MetAlyzer::log2FC(metabObj), ...)'))
        reactCodeHistory(c(reactCodeHistory(),
                           '# Make interactive scatter plot'))
        reactCodeHistory(c(reactCodeHistory(),
                           'MetAlyzer:::plotly_scatter(MetAlyzer::log2FC(metabObj))$Plot'))
        reactCodeHistory(c(reactCodeHistory(),
                           '# Make interactive network diagram with specified parameters'))
        reactCodeHistory(c(reactCodeHistory(),
                           'MetAlyzer:::plotly_network(MetAlyzer::log2FC(metabObj), ...)'))
      }
      
      # Update the slider input, for custom inputs
      updateSliderInput(session, "plotVolcanoLog2FCCutoff",
                        max = floor(max(na.omit(reactLog2FCTbl()$log2FC))))
    } else {
      showModal(modalDialog(
        title = 'Differential analysis failed...',
        'Please select two different sample groups.',
        easyClose = T,
        footer = NULL
      ))
    }
  })
  
  
  # Update choices for sample metadata to include in exported abundance csv file
  observe({
    req(reactMetabObj$metabObj)
    updateSelectInput(session, 'metaChoicesRawAbunExport', choices = colnames(colData(reactMetabObj$metabObj)))
  })
  # Update choices for feature identifiers to include in exported annotation csv file
  featIdTbl <- reactive({
    app_data_dir <- system.file("shinyapp", "data", package = "MetAlyzer")
    readr::read_csv(file.path(app_data_dir, 'BiocratesFeatureTable.csv'))
  })
  output$loadFeatIdChoicesExport <- renderUI({
    annoChoices <- colnames(featIdTbl())
    # Remove columns that will be included in exported file
    annoChoices <- annoChoices[-which(annoChoices %in% c('TrivialName', 'TrivialName_Prior2023', 'Class'))]
    selectInput('featIdChoicesExport', 'Select metabolite annotation(s) to include:',
                choices = annoChoices, multiple = T)
  })
  
  
  # Set plot theme
  # th <- theme_bw(base_size = 15) +
  #   theme(axis.title = element_text(face = 'bold'),
  #         axis.text = element_text(face = 'bold'),
  #         axis.ticks = element_line(linewidth = 0.8),
  #         legend.text = element_text(size = 15))
  
  # Display Log Elements
  # Log of parameters
  output$tblParamLog <- DT::renderDataTable({
    paramLogTbl <- reactAnalysisLog$paramLogTbl
    # Reverse parameter log table to display latest log on top
    paramLogTbl <- paramLogTbl[order(-paramLogTbl$idx),] %>%
      dplyr::mutate(imputation = dplyr::case_when(imputation %in% 'TRUE' ~ 'True',
                                                  imputation %in% 'FALSE' ~ 'False'),
                    featCompleteCutoff = paste0(featCompleteCutoff, '%'),
                    featValidCutoff = paste0(featValidCutoff, '%'))
    colnames(paramLogTbl) <- c('Process', 'Samples Removed', 'Feature Classes Removed',
                               'Count of Features Removed', 'Completeness Cutoff',
                               'Validity Cutoff', 'Validity Statuses', 'Imputation', 'Normalization')
    DT::datatable(paramLogTbl, rownames = F, options = list(dom = 't', ordering = F),
                  selection = list(mode = 'single', target = 'row'), style = 'bootstrap')
  })
  
  # Prepare row index clicked
  clickedRowIdx <- reactive({
    req(input$tblParamLog_rows_selected)
    clickedRow <- input$tblParamLog_rows_selected
    # Because of reversed parameter log table
    sort(reactAnalysisLog$paramLogTbl$idx, decreasing = T)[clickedRow]
  })
  # Expand all panels each time trial (row) is clicked, so that output can be rendered
  observeEvent(clickedRowIdx(), {
    shinyBS::updateCollapse(session, id = 'panelTrialLog',
                            open = c('Features Removed', 'Data distribution', 'Data completeness', 'Quantification status'))
  })
  
  # Log of all removed features
  output$textRmFeatsLog <- renderText({
    req(clickedRowIdx())
    rmFeats <- reactAnalysisLog$rmFeatList[[clickedRowIdx()]]
    # Place it here so that empty removed feature list got length of 0 in parameter log table
    if (length(rmFeats) == 0) {
      rmFeats <- 'None'
    }
    paste(sort(rmFeats), collapse = ' | ')
  })
  
  # Log of processed data distribution
  output$plotDatDistLog <- plotly::renderPlotly({
    req(clickedRowIdx())
    actMetabAggreTbl <- reactAnalysisLog$metabAggreTblList[[clickedRowIdx()]]
    g <- ggplot(actMetabAggreTbl, aes(x=ID, y=Concentration)) +
      geom_boxplot() +
      labs(x = 'Sample', y = 'Metabolite abundance') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # Put abundance on log10 scale for better visualization if data is not normalized 
    if (reactAnalysisLog$paramLogTbl$normalization[clickedRowIdx()] == 'None') {
      g <- g + scale_y_log10()
    }
    plotly::ggplotly(g)
  })
  
  # Log of processed data quantification status
  output$plotQuanStatusLog <- plotly::renderPlotly({
    req(clickedRowIdx())
    actMetabAggreTbl <- reactAnalysisLog$metabAggreTblList[[clickedRowIdx()]]
    smpStatusCountTbl <- actMetabAggreTbl %>%
      dplyr::group_by(ID, Status) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::ungroup()
    # Prepare colors for quantification statuses
    #### Add color for NA quantification status
    status2Color <- c(Valid = '#33a02c', LOQ = '#1f78b4', LOD = '#ff7f00', Invalid = '#e31a1c')
    statusCols <- status2Color[names(status2Color) %in% unique(smpStatusCountTbl$Status)]
    g <- ggplot(smpStatusCountTbl, aes(x = ID, y = Count, fill = Status)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = statusCols) +
        labs(x = "Sample") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    plotly::ggplotly(g)
  })
  
  # Log of processed data completeness
  output$plotDatCompleteLog <- plotly::renderPlotly({
    req(clickedRowIdx())
    actMetabAggreTbl <- reactAnalysisLog$metabAggreTblList[[clickedRowIdx()]]
    smpCompleteCountTbl <- actMetabAggreTbl %>%
      dplyr::mutate(Completeness = dplyr::case_when(!Concentration %in% NA ~ 'Observed', #c(0, NA)
                                                    Concentration %in% NA ~ 'Missing')) %>% #c(0, NA)
      dplyr::group_by(ID, Completeness) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::ungroup()
    ggplotly(
      ggplot(smpCompleteCountTbl, aes(x=ID, y=Count, fill=Completeness)) +
        geom_col(position = 'stack') +
        scale_fill_manual(values = c(Missing = 'grey', Observed = 'black')) +
        labs(x = 'Sample') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    )
  })
  
  # Log of commands executed
  output$textCommands <- renderText({
    req(reactMetabObj$metabObj)
    paste0(reactCodeHistory(), collapse = '\n')
  })
  # Downloading of R session info
  output$downloadSessionInfo <- downloadHandler(
    filename = function() {
      paste0("session_info_", Sys.Date(), '.txt')
    },
    content = function(file) {
      rSessionInfo <- utils::capture.output(utils::sessionInfo())
      writeLines(rSessionInfo, con = file)
    }
  )
  
  
  # Create reactive object for storing static overview plots for downloading
  reactOverviewPlots <- reactiveValues(datDist = NULL, datComplete = NULL, quanStatus = NULL)
  # Show data overviews
  # Data distribution
  output$summDatDist <- renderText({
    req(datOverviewPack()$metabAggreTbl)
    'Is the dataset already log-transformed? If so, normalization should not be applied again.'
  })
  #### Collapsed panel: Workaround for light rendering
  # Need it to render output when panel is collapsed
  outputOptions(output, 'summDatDist', suspendWhenHidden = F)
  output$updateGpColsDatDist <- renderUI({
    req(smpChoicePack()$smpChoiceList)
    smpChoiceList <- smpChoicePack()$smpChoiceList
    # Remove choice groups whose level sizes are one or equal to sample size
    choiceGpSizes <- sapply(smpChoiceList, length)
    rmChoiceGps <- which(choiceGpSizes == 1 | choiceGpSizes == ncol(reactMetabObj$metabObj))
    if (length(rmChoiceGps) != 0) {
      smpChoiceGps <- names(smpChoiceList)[-rmChoiceGps]
    } else {
      smpChoiceGps <- names(smpChoiceList)
    }
    # Make label on left side, rather than top
    div(
      style = "display: flex; align-items: center;",
      tags$label("Color by:", style = "margin-right: 10px; margin-bottom: 15px;"),
      selectInput('gpColsDatDist', NULL,
                  choices = c('None', smpChoiceGps),
                  selected = 'None', multiple = F)
    )
  })
  output$plotDatDist <- plotly::renderPlotly({
    req(datOverviewPack()$metabAggreTbl, input$gpColsDatDist)
    metabAggreTbl <- datOverviewPack()$metabAggreTbl
    if (input$gpColsDatDist == 'None') {
      g <- ggplot(metabAggreTbl, aes(x=ID, y=Concentration)) +
        geom_boxplot()
    } else {
      g <- ggplot(metabAggreTbl, aes(x=ID, y=Concentration, fill=.data[[input$gpColsDatDist]])) +
        geom_boxplot(alpha = 1) +
        scale_fill_brewer(palette = 'Set1')
    }
    g <- g +
      labs(x = 'Sample', y = 'Metabolite abundance') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # Put abundance on log10 scale for better visualization if data is not normalized 
    if (doneNormalization() == 0) {
      g <- g + scale_y_log10()
    }
    # Update reactive object
    reactOverviewPlots$datDist <- g
    
    ggplotly(g)
  })
  
  # Sample metadata
  output$tblSmpMetadat <- DT::renderDataTable({
    # Static
    req(reactOriSmpMetadatTbl())
    DT::datatable(reactOriSmpMetadatTbl(), rownames = F, filter = list(position = 'top', clear = T, plain = F),
                  selection = list(mode = 'single', target = 'row'), style = 'bootstrap',
                  options = list(pageLength = 5))
    
    # Reactive
    # req(datOverviewPack()$smpMetadatTbl)
    # smpMetadatTbl <- datOverviewPack()$smpMetadatTbl
    # DT::datatable(smpMetadatTbl, rownames = F, filter = list(position = 'top', clear = T, plain = F),
    #               selection = list(mode = 'single', target = 'row'), style = 'bootstrap',
    #               options = list(pageLength = 5))
  })
  outputOptions(output, 'tblSmpMetadat', suspendWhenHidden = F)
  #### Collapsed panel: Workaround (SHOULD WORK BUT NOT WORK)
  # Render table outside collapsed panel and display it
  # output$uiTblSmpMetadat <- renderUI({
  #   DT::dataTableOutput('tblSmpMetadat')
  # })
  
  # Data completeness
  output$summDatComplete <- renderText({
    req(datOverviewPack()$featCompleteLvTbl)
    featCompleteLevels <- datOverviewPack()$featCompleteLvTbl$CompleteRatio
    paste(sum(featCompleteLevels < 0.8), 'out of', nrow(reactMetabObj$metabObj),
          'metabolites fail to fulfil the 80% rule, the recommended filtering threshold.',
          'Besides, is there any sample with high level of missingness?')
  })
  #### Collapsed panel: Workaround for light rendering
  # Need it to render output when panel is collapsed
  outputOptions(output, 'summDatComplete', suspendWhenHidden = F)
  output$plotDatComplete <- plotly::renderPlotly({
    req(datOverviewPack()$metabAggreTbl)
    smpCompleteCountTbl <- datOverviewPack()$metabAggreTbl %>%
      dplyr::mutate(Completeness = dplyr::case_when(!Concentration %in% NA ~ 'Observed', #c(0, NA)
                                                    Concentration %in% NA ~ 'Missing')) %>% #c(0, NA)
      dplyr::group_by(ID, Completeness) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::ungroup()
    g <- ggplot(smpCompleteCountTbl, aes(x=ID, y=Count, fill=Completeness)) +
      geom_col(position = 'stack') +
      scale_fill_manual(values = c(Missing = 'grey', Observed = 'black')) +
      labs(x = 'Sample') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # Update reactive object
    reactOverviewPlots$datComplete <- g
    
    ggplotly(g)
  })
  
  # Quantification status
  output$summQuanStatus <- renderText({
    req(datOverviewPack()$featStatusCountTbl)
    featStatusValid <- datOverviewPack()$featStatusCountTbl %>%
      dplyr::mutate(TotalValidCount = ValidCount + LOQCount,
                    ValidRatio = TotalValidCount / TotalCount)
    paste(sum(featStatusValid$ValidRatio < 0.5), 'out of', nrow(reactMetabObj$metabObj),
          'metabolites have fewer than 50% valid quantifications (Valid, LOQ) and',
          'could be filtered out. Besides, is there any sample containing just few valid quantifications?')
  })
  #### Collapsed panel: Workaround for light rendering
  # Need it to render output when panel is collapsed
  outputOptions(output, 'summQuanStatus', suspendWhenHidden = F)
  output$plotQuanStatus <- plotly::renderPlotly({
    req(datOverviewPack()$metabAggreTbl)
    smpStatusCountTbl <- datOverviewPack()$metabAggreTbl %>%
      dplyr::group_by(ID, Status) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::ungroup()
    status2Color <- c('Valid'='#33a02c', 'LOQ'='#1f78b4', 'LOD'='#ff7f00', 'Invalid'='#e31a1c')
    g <- ggplot(smpStatusCountTbl, aes(x = ID, y = Count, fill = Status)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = status2Color) +
        labs(x = "Sample") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # Update reactive object
    reactOverviewPlots$quanStatus <- g    
    
    plotly::ggplotly(g)
})
  
  # Export abundance matrix
  output$downloadRawAbunExport <- downloadHandler(
    filename = paste("abundance_matrix_", format(Sys.Date(), "%Y-%m-%d"), ".csv", sep = ""),
    content = function(file) {
      #### Add option of exporting also processed data?
      #### What inputs MetaboAnalyst needs?
      MetAlyzer::exportAbunValues(reactMetabObj$metabObj,
                                  input$metaChoicesRawAbunExport,
                                  file_path = file)
    }
  )
  # Export metabolite identifiers
  output$downloadFeatIdsExport <- downloadHandler(
    filename = paste("metabolite_identifiers_", format(Sys.Date(), "%Y-%m-%d"), ".csv", sep = ""),
    content = function(file) {
      if (!input$ifUploadedFilePrior2023) {
        exportfeatIdTbl <- featIdTbl()[featIdTbl()$TrivialName %in% rownames(reactMetabObj$metabObj),] %>%
          dplyr::rename(BiocratesName = TrivialName) %>%
          dplyr::select(-TrivialName_Prior2023)
      } else {
        exportfeatIdTbl <- featIdTbl()[featIdTbl()$TrivialName_Prior2023 %in% rownames(reactMetabObj$metabObj),] %>%
          dplyr::rename(BiocratesName = TrivialName_Prior2023) %>%
          dplyr::select(-TrivialName)
      }
      exportfeatIdTbl <- dplyr::select(exportfeatIdTbl, BiocratesName, Class, input$featIdChoicesExport)
      write.csv(exportfeatIdTbl, file)
    }
  )
  # Export differential analysis result table
  output$downloadLog2FC <- downloadHandler(
    filename = function() {
      paste0("Log2FC_results_", Sys.Date(), ".xlsx") 
    },
    content = function(file) {
      data_to_export <- reactLog2FCTbl()
      writexl::write_xlsx(data_to_export, path = file)
    }
  )
  
  # Visualize log2(FC)
  # Give sign before log2(FC) calculation
  output$textLog2FC <- renderText({
    'Please perform differential analysis first.'
  })
  output$textLog2FC_2 <- renderText({
    'Please perform differential analysis first.'
  })
  
  # Update choices for metabolite highlighting
  observe({
    req(featChoices())
    updateSelectInput(session, 'metabChoicesVulcano', choices = featChoices())
  })
  # Create new column for highlighting metabolites in vulcano and scatter plots
  observe({
    req(reactLog2FCTbl())
    selectedChoices <- input$metabChoicesVulcano
    highlightTbl <- reactLog2FCTbl()
    if (!is.null(selectedChoices)) {
      highlightTbl$highlight <- FALSE
      highlightTbl$highlight[highlightTbl$Metabolite %in% selectedChoices] <- TRUE
      highlightTbl$highlight[highlightTbl$Class %in% selectedChoices] <- TRUE
    } else if (is.null(selectedChoices) | input$highlightVulcano) {
      highlightTbl$highlight <- FALSE
    }
    reactVulcanoHighlight(highlightTbl)
  })
  
  # Vulcano plot
  output$plotVolcano <- plotly::renderPlotly({
    req(reactLog2FCTbl())
    if (!input$highlightVulcano) {
      MetAlyzer:::plotly_vulcano(reactLog2FCTbl(),
                                 x_cutoff = input$plotVolcanoLog2FCCutoff,
                                 y_cutoff = as.numeric(input$plotVolcanoPValCutoff))
    } else {
      MetAlyzer:::plotly_vulcano(reactVulcanoHighlight(),
                                 x_cutoff = input$plotVolcanoLog2FCCutoff,
                                 y_cutoff = as.numeric(input$plotVolcanoPValCutoff))
    }
  })
  
  # Scatter plot
  output$plotScatter <- plotly::renderPlotly({
    req(reactLog2FCTbl())
    # if (!input$highlightVulcano) {
    #   plot <- MetAlyzer:::plotly_scatter(reactLog2FCTbl())
    # } else {
    #   plot <- MetAlyzer:::plotly_scatter(reactVulcanoHighlight())
    # }
    plot <- MetAlyzer:::plotly_scatter(reactLog2FCTbl())
    hide_legend(plot$Plot)
  })
  output$plotScatterLegend <- renderImage({
    req(reactLog2FCTbl()) 
    plot <- MetAlyzer:::plotly_scatter(reactLog2FCTbl())
    legend <- plot$Legend
    # A temp file to save the output.
    # This file will be removed later by renderImage
    outfile <- tempfile(fileext='.svg')
    ggsave(file=outfile, plot=legend)
    
    # Return a list containing the filename
    list(src = normalizePath(outfile),
         contentType = 'image/svg+xml',
         width = 571.675,
         height = 400,
         alt = "",
         deleteFile = TRUE)
  }, deleteFile = TRUE)
  
  # Network plot
  network_data <- reactive({
    # Ensure the base data is available before proceeding
    req(reactLog2FCTbl())

    # Call the expensive function just one time
    MetAlyzer:::plotly_network(
      reactLog2FCTbl(),
      values_col_name = input$networkValueColumn,
      exclude_pathways = input$networkExcludePathways,
      metabolite_node_size = input$networkMetaboliteNodeSize,
      connection_width = input$networkConnectionWidth,
      pathway_text_size = input$networkPathwayTextSize,
      pathway_width = input$networkPathwayWidth,
      plot_height = input$networkPlotHeight*100,
      color_scale = input$networkColorScale
    )
  })
  output$plotNetwork <- plotly::renderPlotly({
    # Access the pre-calculated plot from the reactive expression
    req(network_data())
    network_data()$Plot
  })
  
  ### --- Network Nodes ---
  # Update choices for node stats
  observeEvent(reactLog2FCTbl(), {
    node_table <- network_data()$Table
    node_choices <- unique(node_table$Label_nFeatures)
    if ('T14 -- 14 feature(s)' %in% node_choices) {
      updateSelectInput(session, inputId = "selectedNodesVulcano",
                        choices = node_choices, selected = 'T14 -- 14 feature(s)')
    } else {
      updateSelectInput(session, inputId = "selectedNodesVulcano",
                        choices = node_choices, selected = node_choices[1])
    }
  })
  # Prepare stats table for selected nodes to visualize and download
  selected_nodes_data <- reactive({
    req(network_data(), input$selectedNodesVulcano)
    
    # Access the pre-calculated table from the reactive expression
    Table_nodes <- network_data()$Table
    
    # Filter and process the data based on user selection
    Table_nodes %>%
      dplyr::filter(.data$Label_nFeatures %in% input$selectedNodesVulcano) %>%
      tidyr::separate_rows(.data$Metabolites, .data$log2FC, .data$qval, .data$pval, .data$Class, sep = "\\s*;\\s*") %>%
      dplyr::mutate(
        log2FC = as.numeric(.data$log2FC),
        qval = as.numeric(.data$qval),
        pval = as.numeric(.data$pval),
        Metabolite = Metabolites,
        Class = gsub("\"", "", .data$Class),
      ) %>%
      dplyr::filter(.data$Class != "NA") %>%
      dplyr::mutate(Class = as.factor(Class)) %>%
      dplyr::select(-c(x, y, collapsed_count, Label_nFeatures, Shape))
  })
  # Display stats of selected nodes via volcano plot
  output$plotVolcanoNodes <- plotly::renderPlotly({
    req(selected_nodes_data())
    MetAlyzer:::plotly_vulcano(selected_nodes_data(),
                               x_cutoff = input$plotVolcanoLog2FCCutoff,
                               y_cutoff = as.numeric(input$plotVolcanoPValCutoff))
  })
  # Download stats table of selected nodes
  output$downloadNodesExcel <- downloadHandler(
    filename = function() {
      paste0("Volcano_Data_", input$selectedNodesVulcano, "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      writexl::write_xlsx(selected_nodes_data(), path = file)
    }
  )

  # Revert changed network plot style parameters back to default
  observeEvent(input$defaultNetworkPlotStyles, {
    updateSliderInput(session, 'networkPlotHeight', value = 10)
    updateSliderInput(session, 'networkMetaboliteNodeSize', value = 11)
    updateSliderInput(session, 'networkConnectionWidth', value = 1.25)
    updateSliderInput(session, 'networkPathwayTextSize', value = 20)
    updateSliderInput(session, 'networkPathwayWidth', value = 10)
    updateSelectInput(session, 'networkValueColumn', selected = "log2FC")
    updateSelectInput(session, 'networkExcludePathways', selected = character(0))
    updateSelectInput(session, "networkColorScale", selected = "Viridis")
    updateSelectInput(session, "formatNetwork", selected = "html")
  })
  
  
  # Download data overview
  # Data distribution
  output$downloadDatDist <- downloadHandler(
    filename = function() {
      paste0("data_distribution_", Sys.Date(), '.', input$formatDatDist)
    },
    content = function(file) {
      if (input$formatDatDist == "html") {
        final_plot <- ggplotly(reactOverviewPlots$datDist)
        htmlwidgets::saveWidget(final_plot, file, selfcontained = TRUE)
      } else {
        final_plot <- reactOverviewPlots$datDist
        ggsave(filename = file, plot = final_plot, device = input$formatDatDist,
               dpi = 400, units = "cm", width = 32.0, height = 21.0)
      }
    }
  )
  
  # Missing pattern
  output$downloadDatComplete <- downloadHandler(
    filename = function() {
      paste0("missing_pattern_", Sys.Date(), '.', input$formatDatComplete)
    },
    content = function(file) {
      if (input$formatDatComplete == "html") {
        final_plot <- ggplotly(reactOverviewPlots$datComplete)
        htmlwidgets::saveWidget(final_plot, file, selfcontained = TRUE)
      } else {
        final_plot <- reactOverviewPlots$datComplete
        ggsave(filename = file, plot = final_plot, device = input$formatDatComplete,
               dpi = 400, units = "cm", width = 32.0, height = 21.0)
      }
    }
  )
  
  # Quantification status
  output$downloadQuanStatus <- downloadHandler(
    filename = function() {
      paste0("quant_status_", Sys.Date(), '.', input$formatQuanStatus)
    },
    content = function(file) {
      if (input$formatQuanStatus == "html") {
        final_plot <- ggplotly(reactOverviewPlots$quanStatus)
        htmlwidgets::saveWidget(final_plot, file, selfcontained = TRUE)
      } else {
        final_plot <- reactOverviewPlots$quanStatus
        ggsave(filename = file, plot = final_plot, device = input$formatQuanStatus,
               dpi = 400, units = "cm", width = 32.0, height = 21.0)
      }
    }
  )
  
  # Download log2(FC) visuals
  # Vulcano plot
  output$downloadVulcanoPlot <- downloadHandler(
    filename = function() {
      paste0("vulcano_plot_", Sys.Date(), '.', input$formatVulcano)
    },
    content = function(file) {
      if (input$formatVulcano == "html") {
        if (input$highlightVulcano) {
          req(reactVulcanoHighlight())
          final_plot <- MetAlyzer:::plotly_vulcano(reactVulcanoHighlight(), 
                                                   x_cutoff = input$plotVolcanoLog2FCCutoff,
                                                   y_cutoff = as.numeric(input$plotVolcanoPValCutoff))
        } else {
          final_plot <- MetAlyzer:::plotly_vulcano(reactLog2FCTbl(), 
                                                   x_cutoff = input$plotVolcanoLog2FCCutoff,
                                                   y_cutoff = as.numeric(input$plotVolcanoPValCutoff))
        }
        htmlwidgets::saveWidget(final_plot, file, selfcontained = TRUE)
      } else {
        if (input$highlightVulcano) {
          req(reactVulcanoHighlight())
          final_plot <- MetAlyzer:::plot_vulcano(reactVulcanoHighlight(), 
                                                 x_cutoff = input$plotVolcanoLog2FCCutoff,
                                                 y_cutoff = as.numeric(input$plotVolcanoPValCutoff),
                                                 show_labels_for = input$metabChoicesVulcano)
        } else {
          final_plot <- MetAlyzer:::plot_vulcano(reactLog2FCTbl(), 
                                                 x_cutoff = input$plotVolcanoLog2FCCutoff,
                                                 y_cutoff = as.numeric(input$plotVolcanoPValCutoff))
        }
        ggsave(filename = file, plot = final_plot, device = input$formatVulcano,
               dpi = 400, units = "cm", width = 38.0, height = 21.0)
      }
    }
  )
  
  # Scatter plot
  output$downloadScatterPlot <- downloadHandler(
    filename = function() {
      paste0("scatter_plot_", Sys.Date(), '.', input$formatScatter)
    },
    content = function(file) {
      if (input$formatScatter == "html") {
        htmlwidgets::saveWidget(MetAlyzer:::plotly_scatter(reactLog2FCTbl())$Plot, file, selfcontained = TRUE)
      } else {
        #### Make theme similar to vulcano plot for better visualization
        ggsave(filename = file, plot = MetAlyzer::plot_scatter(reactLog2FCTbl()),
               device = input$formatScatter, dpi = 400, units = "cm", width = 38.0, height = 21.0)
      }
    }
  )
  
  # Network diagram
  output$downloadNetworkPlot <- downloadHandler(
    filename = function() {
      paste0("network_plot_", Sys.Date(), '.', input$formatNetwork)
    },
    content = function(file) {
      if (input$formatNetwork == "html") {
        final_plot <- network_data()$Plot
        htmlwidgets::saveWidget(final_plot, file, selfcontained = TRUE)
      } else {
        final_plot <- MetAlyzer::plot_network(reactLog2FCTbl(),
                                              values_col_name = input$networkValueColumn,
                                              exclude_pathways = input$networkExcludePathways,
                                              color_scale = input$networkColorScale
                                              )$Plot
        ggsave(filename = file, plot = final_plot, device = input$formatNetwork,
               dpi = 400, units = "cm", width = 38.0, height = 21.0)
      }
    }
  )

  ### --- Notifications ---
  current_notification_id <- reactiveVal(NULL)
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

  # 1. Notification for the "Process" button
  observeEvent(input$updateProcessing, {
    if (!is.null(current_notification_id())) {
      removeNotification(current_notification_id())
    }
    notification_ui <- div(
      tags$ul(
        style = "padding-left: 20px; font-size: 0.9em;",
        tags$li(tags$b("Removed Samples: "), paste(input$smpChoicesFiltering, collapse = ", ") %||% "None"),
        tags$li(tags$b("Removed Metabolites: "), paste(input$featChoicesFiltering, collapse = ", ") %||% "None"),
        tags$li(tags$b("% Quantified Cutoff: "), input$featCompleteCutoffFiltering),
        tags$li(tags$b("% Valid Cutoff: "), input$featValidCutoffFiltering),
        tags$li(tags$b("Imputation: "), ifelse(input$imputation, "Half-minimum", "None")),
        tags$li(tags$b("Normalization: "), input$normalization)
      )
    )
    new_id <- showNotification(
      ui = notification_ui,
      duration = 10,
      closeButton = TRUE,
      type = "default"
    )
    current_notification_id(new_id)
  })

  # 2. Notification for the "Revert/Default" button
  observeEvent(input$revertProcessing, {
    if (!is.null(current_notification_id())) {
      removeNotification(current_notification_id())
    }
    
    new_id <- showNotification(
      "Processing undone. Parameters reverted to default.",
      duration = 5,
      type = "warning"
    )
    current_notification_id(new_id)
  })

  # Disclaimer for Uploading Data
  #if (Sys.getenv("SHINY_PORT") != "") {
  
    showModal(modalDialog(
        title = tags$b("Data Privacy Disclaimer", style = "color: #d9534f;"),
        HTML("Please ensure that all uploaded data is <b>fully anonymized</b> and contains no sensitive patient information.
            This application is hosted on a shared server (shinyapps.io), and we cannot guarantee the privacy of identifiable data.
            <br><br>
            By proceeding, you confirm that your are aware of this and proceed at your own risk.
            <br><br>
            You can find a guide to set up the application locally on your machine here:
            <br>
            <a href='https://github.com/Lu-Group-UKHD/MetAlyzer/blob/main/vignettes/shiny_app.Rmd'>Github</a>
            <br>
            or in the tutorial video."),
        # Replace the old footer with this new one
        footer = tags$button(
          "I Understand and Accept",
          type = "button",
          class = "btn btn-primary", # This class makes the button blue
          `data-dismiss` = "modal"  # This ensures the button closes the pop-up
        ),
        easyClose = FALSE
    ))
  
  #}
}

shinyApp(ui = ui, server = server)
