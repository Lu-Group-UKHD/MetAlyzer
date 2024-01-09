library(shiny)
library(MetAlyzer)
library(SummarizedExperiment)
library(tidyverse)

# setwd('/Users/qianwu/Desktop/RShiny_Biocrates_DataAnalysis')
# metabObj <- MetAlyzer_dataset(file_path = './data/extraction_data.xlsx', silent = T)

ui <- fluidPage(
  # App title
  titlePanel('Biocrates Metabolomics Analysis'),
  
  tabsetPanel(
    type = 'tabs',
    tabPanel(
      'Preprocessing',
      sidebarLayout(
        sidebarPanel(
          tags$h4('Data uploading', style = 'color:steelblue;font-weight: bold'),
          fileInput('uploadedFile', NULL, multiple = F, accept = '.xlsx',
                    placeholder = 'No .xlsx file selected'),
          # Show filtering options only after file is uploaded
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           tags$h4('Data filtering', style = 'color:steelblue;font-weight: bold'),
                           selectInput('featChoicesFiltering', 'Select feature(s) to remove:',
                                       choices = character(0), multiple = T),
                           selectInput('smpChoicesFiltering', 'Select sample(s) to remove:',
                                       choices = character(0), multiple = T),
                           fluidRow(
                             column(width = 6, actionButton('updateFiltering', 'Filter', width = '100%')),
                             column(width = 6, actionButton('resetFiltering', 'Reset', width = '100%'))
                           ),
                           tags$br(),
                           tags$h4('Log₂(FC) visualization', style = 'color:steelblue;font-weight: bold'),
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
                             column(width = 5, actionButton('computeLog2FC', 'Compute', width = '100%'))
                           ),
                           tags$br(),
                           checkboxGroupInput('plotLog2FC', 'Visualize through:',
                                              choices = c('Volcano plot', 'Scatter plot', 'Network plot'),
                                              inline = T)
          )
        ),
        mainPanel(
          fluidRow(
            column(width = 7, verbatimTextOutput('textConcSumm')),
            column(width = 5, verbatimTextOutput('textQuanSumm'))
          ),
          tableOutput('outputId'),
          plotlyOutput('plotVolcano'),
          fluidRow(
            column(width = 9, plotlyOutput('plotScatter')),
            column(width = 3, plotOutput('plotScatterLegend'))
          ),
          plotOutput('plotNetwork')
        )
      )
    ),
    tabPanel(
      'Something',
      tags$h4('Could be a section for something else')
    )
  )
)

server <- function(input, output, session) {
  # Create reactive objects for storing up-to-date data
  reactMetabObj <- reactiveValues(metabObj = NULL, ori_metabObj = NULL)
  reactLog2FCTab <- reactiveVal()
  reactLog2FCvulcano <- reactiveVal()
  reactP_data <- reactiveVal()
  reactRects_df <- reactiveVal()
  
  # Initialize MetAlyzer SE object
  observeEvent(input$uploadedFile, {
    validUploadedFile <- try(
      metabObj <- MetAlyzer_dataset(file_path = input$uploadedFile$datapath,
                                    sheet = 1, silent = T),
      silent = T)
    if (!is(validUploadedFile, 'try-error')) {
      reactMetabObj$metabObj <- metabObj
      # Make copy of original data for filtering reset
      reactMetabObj$ori_metabObj <- metabObj
      # Monitor whether uploaded file is valid to show further operations on client side
      output$ifValidUploadedFile <- reactive({
        !is.null(reactMetabObj$metabObj)
      })
      outputOptions(output, 'ifValidUploadedFile', suspendWhenHidden = F)
    } else {
      showModal(modalDialog(
        title = 'Uploaded file reading failed...',
        'Is the uploaded .xlsx file generated from MetIDQ™ software?',
        easyClose = T,
        footer = NULL
      ))
      reactMetabObj$metabObj <- NULL
      reactMetabObj$ori_metabObj <- NULL
    }
  })
  
  
  # Prepare choices for feature filtering
  featChoices <- reactive({
    req(reactMetabObj$metabObj)
    list(Class = as.list(unique(rowData(reactMetabObj$metabObj)$metabolic_classes)),
         Metabolite = as.list(rownames(reactMetabObj$metabObj)))
  })
  # Prepare needed information for sample filtering (i.e., choices, dictionary for
  # choice group searches, and list of choice groups containing duplicated choices)
  # and log2(FC) calculation
  smpChoicePack <- reactive({
    req(reactMetabObj$metabObj)
    smpChoiceList <- as.list(colData(reactMetabObj$metabObj))
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
    # (due to design of MetAlyzer::filterMetaData)
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
    # to design of MetAlyzer::filterMetaData)
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
    if ('Metabolism Indicators' %in% unlist(featChoices()) &
        nrow(reactMetabObj$metabObj) == nrow(reactMetabObj$ori_metabObj)) {
      updateSelectInput(session, 'featChoicesFiltering', choices = featChoices(),
                        selected = 'Metabolism Indicators')
    } else {
      updateSelectInput(session, 'featChoicesFiltering', choices = featChoices())
    }
  })
  # Update choices for sample filtering
  observe({
    req(smpChoicePack()$smpChoiceList)
    # Make choices lists so that sole choice in certain choice group can be shown
    smpChoiceList <- lapply(smpChoicePack()$smpChoiceList, function(choices) {
      as.list(choices)
    })
    updateSelectInput(session, 'smpChoicesFiltering', choices = smpChoiceList)
  })
  # Do feature filtering
  observeEvent(input$updateFiltering, {
    if (!is.null(input$featChoicesFiltering)) {
      #### Probably add slider for 'min_percent_valid' and options for 'per_group'
      #### to filter unmet features out by threshold of valid values
      reactMetabObj$metabObj <- MetAlyzer::filterMetabolites(reactMetabObj$metabObj,
                                                             drop_metabolites = input$featChoicesFiltering)
    } else {
      #### Some features will still be filtered out
      reactMetabObj$metabObj <- MetAlyzer::filterMetabolites(reactMetabObj$metabObj,
                                                             drop_metabolites = NULL)
    }
  })
  # Do sample filtering
  observeEvent(input$updateFiltering, {
    req(input$smpChoicesFiltering)
    for (selectedChoice in input$smpChoicesFiltering) {
      selectedChoiceGp <- smpChoicePack()$smpChoices2Gps[selectedChoice]
      # Prepare up-to-date choices of selected choice group to avoid error that
      # metadata of same sample is selected multiple times at once
      choices <- unique(colData(reactMetabObj$metabObj)[[selectedChoiceGp]])
      # Convert NA into character so that it can be searched
      if (any(is.na(choices))) {
        choices <- choices[-which(is.na(choices))]
        choices <- c(choices, 'NA')
      }
      if (!selectedChoiceGp %in% smpChoicePack()$dupSmpChoiceGps) {
        if (selectedChoice %in% choices) {
          if (selectedChoice == 'NA') {
            reactMetabObj$metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                                                !is.na(.data[[selectedChoiceGp]]))
          } else {
            reactMetabObj$metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
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
            reactMetabObj$metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                                                !is.na(.data[[selectedChoiceGp]]))
          } else {
            reactMetabObj$metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                                                !.data[[selectedChoiceGp]] %in% modSelectedChoice)
          }
        }
      }
    }
  })
  # Reset feature and sample filtering
  observeEvent(input$resetFiltering, {
    reactMetabObj$metabObj <- reactMetabObj$ori_metabObj
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
      selectedChoiceGp <- input$smpChoiceGpsLog2FC
      selectedChoices <- c(input$smpChoicesLog2FC_1, input$smpChoicesLog2FC_2)
      if (!'NA' %in% selectedChoices) {
        metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                              .data[[selectedChoiceGp]] %in% selectedChoices)
      } else {
        metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                              is.na(.data[[selectedChoiceGp]]) |
                                                .data[[selectedChoiceGp]] %in% selectedChoices)
      }
      # Use do.call to prepare arguments to be passed due to design of function:
      # deparse(substitute(categorical))
      log2FCTab <- do.call(MetAlyzer::calculate_log2FC, list(metalyzer_se = metabObj,
                                                             categorical = as.symbol(selectedChoiceGp),
                                                             impute_perc_of_min = 0.2, impute_NA = T))

      
      ### Server Side Data Wrangling
      ## Background: Define colors for significance 
      signif_colors=c("#5F5F5F"=1,
                      "#FEBF6E"=0.1,
                      "#EE5C42"=0.05,
                      "#8B1A1A"=0.01)
      
      ## Background: Load polarity data
      polarity_file <- system.file("extdata", "polarity.csv", package = "MetAlyzer")
      
      polarity_df <- utils::read.csv(polarity_file) %>%
        select(.data$Class,
               .data$Polarity) %>%
        mutate(Class = factor(.data$Class),
               Polarity = factor(.data$Polarity, levels = c('LC', 'FIA'))) %>%
        arrange(.data$Polarity)
      
      ## Background: Set class colors
      class_colors <- metalyzer_colors()
      
      names(class_colors) <- levels(polarity_df$Class)
      
      ## Background: Define LC and FIA classes with color
      lc_polarity_df <- filter(polarity_df,
                               .data$Polarity == 'LC',
                               .data$Class %in% log2FCTab$Class)
      lc_colors <- class_colors[which(names(class_colors) %in% lc_polarity_df$Class)]
      fia_polarity_df <- filter(polarity_df,
                                .data$Polarity == 'FIA',
                                .data$Class %in% log2FCTab$Class)
      fia_colors <- class_colors[which(names(class_colors) %in% fia_polarity_df$Class)]
      
      ## Data: Replace NAs
      log2FCTab$log2FC[is.na(log2FCTab$log2FC)] <- 0
      log2FCTab$qval[is.na(log2FCTab$qval)] <- 1
      
      # Data Vulcano: Create Dataframe for vulcano plot
      log2FCvulcano <- log2FCTab
      log2FCvulcano$Class <- as.character(log2FCvulcano$Class)
      log2FCvulcano$Class[log2FCvulcano$qval > 0.05] <- NA
      log2FCvulcano$Class[abs(log2FCvulcano$log2FC) < log2(1.5)] <- NA

      log2FCvulcano$labels <- as.character(log2FCvulcano$Metabolite)
      log2FCvulcano$labels[which(is.na(log2FCvulcano$Class))] <- ""
      
      ## Data Scatter: Add color to data based on significance
      log2FCTab$signif_color <- sapply(log2FCTab$qval, function(q_val) {
        for (t in signif_colors) {
          if (q_val <= t) {
            color <- names(signif_colors)[which(signif_colors == t)]
          }
        }
        return(color)
      })
      
      ## Data Scatter: Add pseudo x-value to data as a order of metabolites
      ordered_classes <- c(names(lc_colors), names(fia_colors))
      p_data <- lapply(ordered_classes, function(class) {
        log2FCTab %>%
          filter(.data$Class == class) %>%
          bind_rows(data.frame(Class = rep(NA, 5)))
      }) %>%
        bind_rows()
      p_data <- bind_rows(data.frame(Class = rep(NA, 5)), p_data)
      p_data$x <- seq(nrow(p_data))
      p_data <- filter(p_data, !is.na(.data$Class))
      
      ## Data Scatter: Determine labels
      signif_p_data <- filter(p_data, .data$signif_color != names(signif_colors)[1])
      
      labels <- sapply(p_data$Metabolite, function(m) {
        m <- as.character(m)
        label <- if_else(m %in% signif_p_data$Metabolite, m, "")
        return(label)
      })
      
      ## Legend Scatter: Significance color
      signif_colors <- sort(signif_colors, decreasing = TRUE)
      signif_labels <- list()
      for (i in seq_along(signif_colors)) {
        t <- signif_colors[i]
        names(t) <- NULL
        if (i < length(signif_colors)) {
          t2 <- signif_colors[i+1]
          names(t2) <- NULL
          label <- bquote(.(t) ~ "\u2265 q-value >" ~ .(t2))
        } else {
          label <- bquote(.(t) ~ "\u2265 q-value")
        }
        signif_labels[[i]] <- label
      }
      
      breaks <- c(names(lc_colors), names(fia_colors))
      values <- c(lc_colors,fia_colors)
      names(values) <- NULL
      
      ## Background Scatter: Create data for background rects
      rects_df <- p_data %>%
        group_by(.data$Class) %>%
        summarise(Start = min(.data$x)-1,
                  End = max(.data$x)+1,
                  Color = class_colors[unique(.data$Class)],
                  n = n())
      rects_df$Class <- factor(rects_df$Class, levels = breaks)
      rects_df$Technique <- sapply(rects_df$Class, function(c) {
        if (c %in% names(lc_colors)) {
          technique <- 'LC'
        } else if (c %in% names(fia_colors)) {
          technique <- 'FIA'
        } else {
          technique <- NA
        }
        return(technique)
      })
      
      ## Background Scatter: Determine border line between last LC and first FIA class
      lc_fia_border <- p_data %>%
        filter(.data$Class %in% names(lc_colors)) %>%
        select(.data$x) %>%
        max()
      
      
      reactP_data(p_data)
      reactRects_df(rects_df)
      
      reactLog2FCvulcano(log2FCvulcano)
      reactLog2FCTab(log2FCTab)
    } else {
      showModal(modalDialog(
        title = 'Log₂(FC) computation failed...',
        'Please select two sample groups.',
        easyClose = T,
        footer = NULL
      ))
    }
  })
  
  
  # Show statistics of data concentration values and quantification statuses 
  output$textConcSumm <- renderPrint({
    req(reactMetabObj$metabObj)
    MetAlyzer::summarizeConcValues(reactMetabObj$metabObj)
  })
  output$textQuanSumm <- renderPrint({
    req(reactMetabObj$metabObj)
    MetAlyzer::summarizeQuantData(reactMetabObj$metabObj)
  })
  
  output$outputId <- renderTable({
    dim(reactMetabObj$metabObj)
  })
  
  # Visualize log2(FC)
  output$plotVolcano <- renderPlotly({
    req(reactLog2FCTab())
    if ('Volcano plot' %in% input$plotLog2FC) {
      ## Plot: Create ggplot object
      p_fc <- ggplot(reactLog2FCvulcano(),
                     aes(x = .data$log2FC,
                         y = -log10(.data$qval),
                         color = .data$Class,
                         label = labels)) +
        geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="black", linetype="dashed") +
        geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") +
        geom_point(size = 1, aes(text = paste0(Metabolite, "\nClass: ", Class, "\nlog2 Fold Change: ", round(log2FC, digits=5), "\np-value: ", round(pval, digits=5)))) +
        scale_color_manual('Classes',
                           breaks = breaks,
                           values = values,
                           drop = FALSE,
                           guide = guide_legend(override.aes = list(size = 2),
                                                order=2, ncol = 2)) +
        theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
              legend.key = element_rect(fill = 'white')) +
        labs(x = 'log2(FC)', y = "-log10(p)")
      
      ## Interactive: Create interactive plot
      p <- ggplotly(p_fc, tooltip = "text")
      p
    }
  })
  output$plotScatter <- renderPlotly({
    req(reactLog2FCTab())
    if ('Scatter plot' %in% input$plotLog2FC) {
      pre_plot <- ggplot(reactP_data(),
                         aes(x = .data$x,
                             y = .data$log2FC,
                             color = .data$signif_color)) +
        geom_point(size = 0.5)
      
      ylims <- ggplot_build(pre_plot)$layout$panel_params[[1]]$y.range
      
      ## Plot: Create ggplot object
      p_fc <- ggplot(reactP_data(),
                     aes(x = .data$x,
                         y = .data$log2FC,
                         color = .data$signif_color)) +
        geom_rect(data = reactRects_df(),
                  inherit.aes = FALSE,
                  aes(xmin = .data$Start, xmax = .data$End,
                      ymin = ylims[1], ymax = ylims[2],
                      fill = .data$Class,
                      text = paste0(Class, "\n", Technique, "\nNumber of Metabolites: ", n)),
                  show.legend = TRUE,
                  alpha = 0.4) +
        geom_vline(xintercept = 0, linewidth = 0.5, color = 'black') +
        geom_vline(xintercept = lc_fia_border+3, linewidth = 0.5, color = 'black', linetype="dotted") +
        geom_hline(yintercept = 0, linewidth = 0.5, color = 'black') +
        geom_point(size = 0.5, aes(text = paste0(Metabolite, 
                                                 "\nClass: ", Class, 
                                                 "\nlog2 Fold Change: ", round(log2FC, digits=5),  
                                                 "\np-value: ", round(pval, digits=5)))) + 
        scale_color_manual(paste0('Significance\n(linear model fit with FDR correction)'),
                           labels = signif_labels,
                           values = names(signif_colors),
                           guide = guide_legend(order=1)) +
        scale_fill_manual('Classes',
                          breaks = breaks,
                          values = values,
                          drop = FALSE,
                          guide = guide_legend(override.aes = list(alpha = 0.5),
                                               order=2, ncol = 2)) +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold.italic', hjust = 0.5),
              legend.key = element_rect(fill = 'white'),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line('#ECECEC'),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_line('#ECECEC'),
              panel.background = element_blank()) +
        labs(x = 'Metabolites')
      
      ## Interactive: Create interactive plot
      p <- ggplotly(p_fc, tooltip = "text", showlegend = FALSE)
      p <- hide_legend(p)
      print("check")
      p
    }
  })
  output$plotScatterLegend <- renderPlot({
    req(reactLog2FCTab())
    if ('Scatter plot' %in% input$plotLog2FC) {
      pre_plot <- ggplot(reactP_data(),
                         aes(x = .data$x,
                             y = .data$log2FC,
                             color = .data$signif_color)) +
        geom_point(size = 0.5)
      
      ylims <- ggplot_build(pre_plot)$layout$panel_params[[1]]$y.range
      
      ## Plot: Create ggplot object
      p_fc <- ggplot(reactP_data(),
                     aes(x = .data$x,
                         y = .data$log2FC,
                         color = .data$signif_color)) +
        geom_rect(data = reactRects_df(),
                  inherit.aes = FALSE,
                  aes(xmin = .data$Start, xmax = .data$End,
                      ymin = ylims[1], ymax = ylims[2],
                      fill = .data$Class,
                      text = paste0(Class, "\n", Technique, "\nNumber of Metabolites: ", n)),
                  show.legend = TRUE,
                  alpha = 0.4) +
        geom_vline(xintercept = 0, linewidth = 0.5, color = 'black') +
        geom_vline(xintercept = lc_fia_border+3, linewidth = 0.5, color = 'black', linetype="dotted") +
        geom_hline(yintercept = 0, linewidth = 0.5, color = 'black') +
        geom_point(size = 0.5, aes(text = paste0(Metabolite, 
                                                 "\nClass: ", Class, 
                                                 "\nlog2 Fold Change: ", round(log2FC, digits=5),  
                                                 "\np-value: ", round(pval, digits=5)))) + 
        scale_color_manual(paste0('Significance\n(linear model fit with FDR correction)'),
                           labels = signif_labels,
                           values = names(signif_colors),
                           guide = guide_legend(order=1)) +
        scale_fill_manual('Classes',
                          breaks = breaks,
                          values = values,
                          drop = FALSE,
                          guide = guide_legend(override.aes = list(alpha = 0.5),
                                               order=2, ncol = 2)) +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold.italic', hjust = 0.5),
              legend.key = element_rect(fill = 'white'),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line('#ECECEC'),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_line('#ECECEC'),
              panel.background = element_blank()) +
        labs(x = 'Metabolites')
      
      # Grab Legend ggplot
      tmp <- ggplot_gtable(ggplot_build(p_fc))
      leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
      legend <- tmp$grobs[[leg]]
      legend <- grid.arrange(legend, ncol=1)
      legend
    }
  })
  
  output$plotNetwork <- renderPlot({
    req(reactLog2FCTab())
    if ('Network plot' %in% input$plotLog2FC) {
      MetAlyzer::plot_network(reactLog2FCTab())
    }
  })
}

shinyApp(ui = ui, server = server)
