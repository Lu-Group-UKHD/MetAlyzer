library(limma)
library(vsn)
library(MetAlyzer)
library(plotly)
library(gridExtra)
library(viridis)
library(viridisLite)
library(SummarizedExperiment)
library(tidyverse)

#' @title Calculate log2 fold change
#' @description This function calculates log2(FC), p-values, and adjusted p-values
#' of the data using limma. Note that the data has to be filtered and log2 transformed
#' already.
#' 
#' @param metalyzer_se A Metalyzer object
#' @param categorical A character specifying the column containing two groups
calc_log2FC <- function(metalyzer_se, categorical) {
  aggregated_data <- metalyzer_se@metadata$aggregated_data
  # Prepare abundance data
  #### Do not know why Metabolite column is grouped by MetAlyzer
  feat_data <- dplyr::ungroup(aggregated_data) %>%
    dplyr::select(Metabolite, ID, Concentration) %>%
    tidyr::pivot_wider(names_from = 'Metabolite', values_from = 'Concentration')
  # Prepare sample metadata of interest
  smp_metadata <- colData(metalyzer_se) %>%
    tibble::as_tibble(rownames = 'ID')
  # Use original column names whose spaces are not replaced with '.'
  colnames(smp_metadata) <- c('ID', colnames(colData(metalyzer_se)))
  smp_metadata <- dplyr::select(smp_metadata, ID, all_of(categorical))
  # Combine abundance data and sample metadata to ensure matched information
  combined_data <- dplyr::left_join(feat_data, smp_metadata, by = 'ID')
  # Retrieve data matrix and sample metadata from combined data to conduct limma
  data_mat <- combined_data[, seq_len(ncol(feat_data))] %>%
    tibble::column_to_rownames('ID') %>%
    t()
  group_vec <- combined_data[, ncol(feat_data)+1, drop = T]
  # Sanity check if specified categorical can split data into two groups
  if (length(unique(group_vec)) != 2) {
    stop("The specified categorical cannot split data into two groups.")
  }
  
  # Compute log2(FC), p-values, and adjusted p-values using limma
  design <- model.matrix(~ group_vec)
  fit <- limma::lmFit(data_mat, design = design)
  fit <- limma::eBayes(fit)
  log2FCRes <- limma::topTable(fit, coef = 2, number = Inf) %>%
    tibble::rownames_to_column('Metabolite') %>%
    dplyr::select(Metabolite, logFC, P.Value, adj.P.Val) %>%
    dplyr::rename(log2FC = logFC, pval = P.Value, qval = adj.P.Val)
  # Combined all information into a table
  group_info <- combined_data[, c(1, ncol(feat_data)+1)]
  log2FCTab <- dplyr::left_join(aggregated_data, group_info, by = 'ID') %>%
    dplyr::left_join(log2FCRes, by = 'Metabolite') %>%
    dplyr::select(Metabolite, Class, log2FC, pval, qval) %>%
    dplyr::distinct(Metabolite, .keep_all = TRUE)
  metalyzer_se@metadata$log2FC <- log2FCTab
  return(metalyzer_se)
}

#' @title Plotly Log2FC Scatter Plot
#' @description This function returns a list with interactive 
#' scatterplot based on log2 fold change data.
#' 
#' @param Log2FCTab A data frame containing log2 fold change data
plotly_scatter <- function(Log2FCTab) {
  ### Data Wrangling
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
                           .data$Class %in% Log2FCTab$Class)
  lc_colors <- class_colors[which(names(class_colors) %in% lc_polarity_df$Class)]
  fia_polarity_df <- filter(polarity_df,
                            .data$Polarity == 'FIA',
                            .data$Class %in% Log2FCTab$Class)
  fia_colors <- class_colors[which(names(class_colors) %in% fia_polarity_df$Class)]
  
  ## Legend: Manage breaks and values
  breaks <- c(names(lc_colors), names(fia_colors))
  values <- c(lc_colors,fia_colors)
  names(values) <- NULL
  
  ## Data: Replace NAs
  Log2FCTab$log2FC[is.na(Log2FCTab$log2FC)] <- 0
  Log2FCTab$qval[is.na(Log2FCTab$qval)] <- 1
  
  ## Data: Add color to data based on significance
  Log2FCTab$signif_color <- sapply(Log2FCTab$qval, function(q_val) {
    for (t in signif_colors) {
      if (q_val <= t) {
        color <- names(signif_colors)[which(signif_colors == t)]
      }
    }
    return(color)
  })
  
  ## Data: Add pseudo x-value to data as a order of metabolites
  ordered_classes <- c(names(lc_colors), names(fia_colors))
  p_data <- lapply(ordered_classes, function(class) {
    Log2FCTab %>%
      filter(.data$Class == class) %>%
      bind_rows(data.frame(Class = rep(NA, 5)))
  }) %>%
    bind_rows()
  p_data <- bind_rows(data.frame(Class = rep(NA, 5)), p_data)
  p_data$x <- seq(nrow(p_data))
  p_data <- filter(p_data, !is.na(.data$Class))
  
  ## Legend: Significance color
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
  
  ## Background: Create data for background rects
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
  
  ## Background: Determine border line between last LC and first FIA class
  lc_fia_border <- p_data %>%
    filter(.data$Class %in% names(lc_colors)) %>%
    select(.data$x) %>%
    max()
  
  ylims <- c(min(Log2FCTab$log2FC) - 0.75, max(Log2FCTab$log2FC) + 0.75)
  
  ## Plot: Create ggplot object
  p_scatter <- ggplot(p_data,
                      aes(x = .data$x,
                          y = .data$log2FC,
                          color = .data$signif_color)) +
    geom_rect(data = rects_df,
              inherit.aes = FALSE,
              aes(xmin = .data$Start, xmax = .data$End,
                  ymin = ylims[1], ymax = ylims[2],
                  fill = .data$Class,
                  text = paste0(Class, "\nTechnique: ", Technique, "\nNumber of metabolites: ", n)),
              show.legend = TRUE,
              alpha = 0.4) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = 'black') +
    geom_vline(xintercept = lc_fia_border+3, linewidth = 0.5, color = 'black', linetype="dotted") +
    geom_hline(yintercept = 0, linewidth = 0.5, color = 'black') +
    geom_point(size = 1, aes(text = paste0(Metabolite, 
                                           "\nClass: ", Class, 
                                           "\nLog2(FC): ", round(log2FC, digits=2),  
                                           "\nAdj. p-value: ", round(qval, digits=4),
                                           "\nP-value: ", round(pval, digits=4)))) + 
    scale_color_manual(paste0('Significance\n(linear model fit with FDR correction)'),
                       labels = signif_labels,
                       values = names(signif_colors),
                       guide = guide_legend(order=1)) +
    scale_fill_manual('Class',
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
    labs(x = 'Metabolite', y = 'Log2(FC)')
  
  ## Interactive: Create interactive plot
  plotly_plot <- ggplotly(p_scatter, tooltip = "text", showlegend = FALSE)
  
  # TODO: Find god to figure out why this is not working
  ## Highlight Metabolites by changing symbol
  #if ("highlight" %in% colnames(p_data)) {
  #  for (i in 1:nrow(p_data)) {
  #    if (isTRUE(p_data$highlight[i])) {
  #      print(p_data[i,])
  #      x_val <- p_data$x[i]
  #      y_val <- p_data$log2FC[i]
  #      color_hex <- p_data$signif_color[i]
  #            
  #      plotly_plot <- plotly_plot %>% add_trace(
  #        x = x_val,
  #        y = y_val,
  #        type = "scatter",
  #        mode = "markers",
  #        marker = list(size = 5)
  #      )
  #    }
  #  }
  #}
  
  # Grab Legend ggplot
  scatter_legend <- ggpubr::get_legend(p_scatter)
  legend <- grid.arrange(scatter_legend, ncol=1)
  
  return(list(Plot = plotly_plot, Legend = legend))
}

#' @title Plotly Log2FC Vulcano Plot
#' @descritpion This function returns a list with interactive 
#' vulcanoplot based on log2 fold change data.
#' 
#' @param Log2FCTab A data frame containing log2 fold change data
#' @param cutoff_y A numeric value specifying the cutoff for q-value
#' @param cutoff_x A numeric value specifying the cutoff for log2 fold change
plotly_vulcano <- function(Log2FCTab, cutoff_y = 0.05, cutoff_x = 1.5) {
  # Make Colors unique for each class
  polarity_file <- system.file("extdata", "polarity.csv", package = "MetAlyzer")
  polarity_df <- utils::read.csv(polarity_file) %>%
  select(.data$Class,
          .data$Polarity) %>%
  mutate(Class = factor(.data$Class),
          Polarity = factor(.data$Polarity, levels = c('LC', 'FIA'))) %>%
  arrange(.data$Polarity)
  class_colors <- metalyzer_colors()
  names(class_colors) <- levels(polarity_df$Class)
  ## Data: Replace NAs
  Log2FCTab$log2FC[is.na(Log2FCTab$log2FC)] <- 0
  Log2FCTab$qval[is.na(Log2FCTab$qval)] <- 1
  
  if("highlight" %in% colnames(Log2FCTab)) {
    # rename factor levels for improving Legend understanding
    Log2FCTab$highlight <- factor(Log2FCTab$highlight, levels = c(FALSE, TRUE), labels = c("Other metabolites", "Highlighted metabolite(s)"))
  
    ## Plot: Create vulcano ggplot object with highlighted points
    p_fc_vulcano_highlighted <- ggplot(Log2FCTab %>%
                                         arrange(desc(highlight)),
                                       aes(x = .data$log2FC,
                                           y = -log10(.data$qval),
                                           color = .data$highlight)) +
      geom_point(size = 1.5, aes(text = paste0(Metabolite, 
                                               "\nClass: ", Class, 
                                               "\nLog2(FC): ", round(log2FC, digits=2), 
                                               "\nAdj. p-value: ", round(qval, digits=4),
                                               "\nP-value: ", round(pval, digits=4)))) +
      geom_vline(xintercept=c(-cutoff_x, cutoff_x), col="black", linetype="dashed") +
      geom_hline(yintercept=-log10(cutoff_y), col="black", linetype="dashed") +
      scale_color_manual('',
                         breaks = c("Other metabolites", "Highlighted metabolite(s)"),
                         values = c("#d3d3d3","#56070C")) +
      theme_bw() +
      labs(x = 'Log2(FC)', y = "-Log10(p-value)")
    
    ## Interactive: Create interactive plot
    p_vulcano_highlighted <- ggplotly(p_fc_vulcano_highlighted, tooltip = "text")

    return(p_vulcano_highlighted)
  } else {
    # Data Vulcano: Prepare Dataframe for vulcano plot
    Log2FCTab$Class <- as.character(Log2FCTab$Class)
    Log2FCTab$Class[Log2FCTab$qval > cutoff_y] <- "Not Significant"
    Log2FCTab$Class[abs(Log2FCTab$log2FC) < cutoff_x] <- "Not Significant"

    breaks <- unique(Log2FCTab$Class)
    values <- class_colors[names(class_colors) %in% Log2FCTab$Class]

    ## Plot: Create vulcano ggplot object
    p_fc_vulcano <- ggplot(Log2FCTab,
                           aes(x = .data$log2FC,
                               y = -log10(.data$qval),
                               color = .data$Class)) +
      geom_vline(xintercept=c(-cutoff_x, cutoff_x), col="black", linetype="dashed") +
      geom_hline(yintercept=-log10(cutoff_y), col="black", linetype="dashed") +
      geom_point(size = 1.5, aes(text = paste0(Metabolite, 
                                               "\nClass: ", Class, 
                                               "\nLog2(FC): ", round(log2FC, digits=2), 
                                               "\nAdj. p-value: ", round(qval, digits=4),
                                               "\nP-value: ", round(pval, digits=4)))) +
      scale_color_manual('Class',
                         breaks = breaks,
                         values = values,
                         drop = FALSE) +
      theme_bw() +
      labs(x = 'Log2(FC)', y = "-Log10(p-value)")
    
    ## Interactive: Create interactive plot
    p_vulcano <- ggplotly(p_fc_vulcano, tooltip = "text")
    
    return(p_vulcano)
  }
}

#' @title Plotly Log2FC Network Plot
#' @description This function returns a list with interactive 
#' networkplot based on log2 fold change data.
#' 
#' @param Log2FCTab A data frame containing log2 fold change data
#' @param q_value A numeric value specifying the cutoff for q-value
#' @param metabolite_node_size The text size of the metabolite Nodes
#' @param connection_width The line width of connections between metabolites
#' @param pathway_text_size The text size of pathway annotations
#' @param pathway_width The line width of pathway-specific connection coloring
#' @param plot_height The height of the Plot in pixel [px]
plotly_network <- function(Log2FCTab,
                           q_value=0.05,
                           metabolite_node_size=11,
                           connection_width=1.25,
                           pathway_text_size=20,
                           pathway_width=10,
                           plot_height=800) {
  pathway_file <- get_data_file_path("Pathway_120924.xlsx")
  ## Read network nodes, edges and annotations
  pathways <- read_named_region(pathway_file, "Pathways_Header")
  invalid_annotations <- which(
    is.na(pathways$Label) |
      duplicated(pathways$Label) |
      is.na(pathways$x) |
      is.na(pathways$y) |
      is.na(pathways$Color)
  )
  if (length(invalid_annotations) > 0) {
    # print warning and remove
    cat("Warning: Removing", length(invalid_annotations), "invalid pathways.\n")
    pathways <- pathways[-invalid_annotations, ]
  }
  rownames(pathways) <- pathways$Label
  
  nodes <- read_named_region(pathway_file, "Metabolites_Header")
  nodes$Pathway[is.na(nodes$Pathway)] <- ""
  invalid_nodes <- which(
    is.na(nodes$Label) |
      duplicated(nodes$Label) |
      is.na(nodes$x) |
      is.na(nodes$y) |
      !nodes$Pathway %in% c(rownames(pathways), "")
  )
  if (length(invalid_nodes) > 0) {
    # print warning and remove
    cat("Warning: Removing", length(invalid_nodes), "invalid nodes.\n")
    nodes <- nodes[-invalid_nodes, ]
  }
  rownames(nodes) <- nodes$Label
  # Remove #1 at the end
  nodes$Label <- gsub("#[0-9]+", "", nodes$Label)
  
  edges <- read_named_region(pathway_file, "Connections_Header")
  invalid_edges <- which(
    !edges$Node1 %in% rownames(nodes) |
      !edges$Node2 %in% rownames(nodes) |
      edges$Node1 == edges$Node2
  )
  if (length(invalid_edges) > 0) {
    # print warning and remove
    cat("Warning: Removing", length(invalid_edges), "invalid connections.\n")
    edges <- edges[-invalid_edges, ]
  }
  
  edges$x_start <- nodes[edges$Node1, "x"]
  edges$y_start <- nodes[edges$Node1, "y"]
  edges$x_end <- nodes[edges$Node2, "x"]
  edges$y_end <- nodes[edges$Node2, "y"]
  edges$Color <- sapply(rownames(edges), function(rowname) {
    from <- edges[rowname, "Node1"]
    to <- edges[rowname, "Node2"]
    from_pathway <- nodes[from, "Pathway"]
    to_pathway <- nodes[to, "Pathway"]
    color <- NA
    if (from_pathway == to_pathway & from_pathway != "") {
      color <- pathways[from_pathway, "Color"]
    }
    return(color)
  })
  
  ## Add log2FC to nodes_df
  signif_df <- filter(Log2FCTab,
                      !is.na(.data$log2FC),
                      !is.na(.data$qval),
                      .data$qval <= q_value)
  
  nodes$FC_thresh <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
    tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
    if (nrow(tmp_df) > 0) {
      # Alteast 1 significantly changed
      l2fc <- sum(tmp_df$log2FC) / nrow(tmp_df)
    } else if (any(m_vec %in% Log2FCTab$Metabolite)) {
      # Not significantly changed but measured
      l2fc <- 0
    } else {
      # Not measured
      l2fc <- NA
    }
    return(l2fc)
  })
  
  ## Add p-value to nodes_df
  nodes$q_value <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
    tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
    if (nrow(tmp_df) > 0) {
      # Alteast 1 significantly changed
      qval <- sum(tmp_df$qval) / nrow(tmp_df)
    } else if (any(m_vec %in% Log2FCTab$Metabolite)) {
      # Not significantly changed but measured
      qval <- sum(Log2FCTab$qval[which(Log2FCTab$Metabolite %in% m_vec)]) / length(m_vec)
    } else {
      # Not measured
      qval <- NA
    }
    return(qval)
  })
  
  ## Draw network
  # Create a plot of the network using ggplotly
  
  # Preparing Hexcodes for Annotation Colors
  nodes$color <- sapply(nodes$FC_thresh, function(value) {
    if (is.na(value)) {
      return("grey")
    } else {
      # Using the viridis color scale, adjust 'option' based on your preference
      color_scale <- viridis(10, option = "D")
      nodes_range <- na.omit(nodes$FC_thresh)
      
      color_index <- findInterval(value, seq(min(nodes_range), max(nodes_range)+0.1, length.out = length(color_scale) + 1))
      return(color_scale[color_index])
    }
  })
  
  # Prepare the Edges List
  area_shapes <- list()
  for (i in seq_along(edges$Color)) {
    if (is.na(edges$Color[i])) {
      edges$Color[i] <- "white"   # Assign edges without pathway white
    }
  }
  # Create background of edges
  for (i in 1:nrow(edges)) {
    area_shape <- list(
      type = "line",
      x0 = edges$x_start[i],
      y0 = edges$y_start[i],
      x1 = edges$x_end[i],
      y1 = edges$y_end[i],
      line = list(
        color = edges$Color[i],
        width = pathway_width
      )
    )
    area_shapes[[i]] <- area_shape
  }
  # Create the edges
  edge_shapes <- list()
  for (i in 1:nrow(edges)) {
    edge_shape <- list(
      type = "line",
      x0 = edges$x_start[i],
      y0 = edges$y_start[i],
      x1 = edges$x_end[i],
      y1 = edges$y_end[i],
      line = list(
        color = "grey",
        width = connection_width
      )
    )
    edge_shapes[[i]] <- edge_shape
  }
  edges_area_combined <- c(area_shapes, edge_shapes)
  
  # Create the nodes
  network <- plot_ly(nodes,
                     x = nodes$x,
                     y = nodes$y,
                     type = "scatter",
                     mode = "markers",
                     marker = list(color = nodes$FC_thresh, 
                                   colorbar = list(title = ""),
                                   colorscale='Viridis',
                                   showscale = TRUE),
                     height = plot_height)
  
  # Add the edges
  p_network <- layout(
    network,
    title = 'Log2(FC) with FDR Correction',
    shapes = edges_area_combined,
    xaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
    yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
    hovermode = FALSE) 
  
  # Add annotations over the nodes
  for (i in 1:nrow(nodes)) {
    p_network <- p_network %>% add_annotations(
      text = nodes$Label[i],
      x = nodes$x[i],
      y = nodes$y[i],
      arrowhead = 0,
      font = list(size = metabolite_node_size, color = "white"),
      ax = 0,
      ay = 0,
      bgcolor = nodes$color[i],
      opacity = 1,
      hovertext = paste0("log2 Fold Change: ", round(nodes$FC_thresh[i], 5),
                         "\nPathway: ", nodes$Pathway[i],
                         "\nadj. p-value: ", round(nodes$q_value[i], 5))
    )
  }
  for (i in 1:nrow(nodes)) {
    p_network <- p_network %>% add_annotations(
      text = nodes$Label[i],
      x = nodes$x[i],
      y = nodes$y[i],
      arrowhead = 0,
      font = list(size = metabolite_node_size, color = "white"),
      ax = 0,
      ay = 0,
      bgcolor = nodes$color[i],
      opacity = 1,
      hovertext = paste0("log2 Fold Change: ", round(nodes$FC_thresh[i], 5),
                         "\nPathway: ", nodes$Pathway[i],
                         "\nadj. p-value: ", round(nodes$q_value[i], 5))
    )
  }
  
  # Add annotations for the pathways
  for (i in 1:nrow(pathways)) {
    p_network <- p_network %>% add_annotations(
      text = pathways$Pathway[i],
      x = pathways$x[i],
      y = pathways$y[i],
      xref = "x",
      yref = "y",
      showarrow = FALSE,
      font = list(size = pathway_text_size, color = pathways$Color[i])
    )
  }
  
  return(p_network)
}


#' @title Read Named Regions
#' @description This function reads in the named regions of an excel file.
#'
#' @param file_path The file path of the file
#' @param named_region The region name u want to read in
#' @keywords internal
read_named_region <- function(file_path, named_region) {
  full_sheet <- openxlsx::read.xlsx(
    file_path,
    sheet = 1,
    colNames = FALSE,
    skipEmptyRows = FALSE,
    skipEmptyCols = FALSE,
  )
  full_sheet[nrow(full_sheet) + 1, ] <- NA
  header <- colnames(openxlsx::read.xlsx(
    file_path,
    namedRegion = named_region
  ))
  coordinates <- lapply(header, function(col_name) {
    data.frame(which(full_sheet == col_name, arr.ind = TRUE))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(row, col) %>%
    dplyr::group_by(row) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(n == length(header))
  
  header_row <- unique(coordinates$row)
  first_row <- header_row + 1
  cols <- coordinates$col
  df <- full_sheet[
    first_row:nrow(full_sheet),
    cols
  ]
  colnames(df) <- header
  last_row <- min(which(rowSums(is.na(df)) == length(header))) - 1
  df <- df[1:last_row, ]
  
  for (numeric_col in c("x", "y", "Radius")) {
    if (numeric_col %in% header) {
      df[, numeric_col] <- as.numeric(df[, numeric_col])
    }
  }
  for (trim_col in c("Label", "Pathway", "Color", "Node1", "Node2")) {
    if (trim_col %in% header) {
      df[, trim_col] <- stringr::str_trim(df[, trim_col])
    }
  }
  rownames(df) <- NULL
  return(df)
}


#' @title Impute aggregated data in SE MetAlyzer object
#' @description This function imputes zero-valued concentrations (missing values)
#' using half-minimum (HM). If all values are zeros, they are set to NA. The imputed
#' values are stored in column 'Concentration' of aggregated data.
#'
#' @param metalyzer_se A MetAlyzer object
#' @param impute_NA A logical value whether to impute NA values
#' @import dplyr
data_imputation <- function(metalyzer_se) { #impute_NA = F
  aggregated_data <- metalyzer_se@metadata$aggregated_data
  # Prepare HM values for corresponding features
  #### Note that aggregated data is already grouped by Metabolite by MetAlyzer
  HM_feats <- dplyr::filter(aggregated_data, Concentration > 0, !is.na(Concentration)) %>%
    dplyr::summarise(HM = min(Concentration) * 0.5)
  # Do HM imputation, which replaces missing values with half of minimum of observed
  # values in corresponding variables
  aggregated_data <- dplyr::left_join(aggregated_data, HM_feats, by = 'Metabolite') %>%
    #### How we deal with those NA values? Let them stay NA for now
    dplyr::mutate(Concentration = dplyr::case_when(Concentration %in% 0 ~ HM,
                                                   !Concentration %in% 0 ~ Concentration)) %>%
    dplyr::select(-HM) %>%
    dplyr::ungroup()
  metalyzer_se@metadata$aggregated_data <- aggregated_data
  return(metalyzer_se)
}

#' @title Generalized log2 transform data
#' @description This function conducts generalized log2 transformation on data
#'
#' @param x A matrix or a vector containing numerical values
glog2 <- function(x) {
  (asinh(x)-log(2))/log(2)
}

#' @title Normalize aggregated data in SE MetAlyzer object
#' @description This function normalizes concentration values among samples using
#' total ion count (TIC) normalization, variance stabilizing normalization (VSN),
#' or median normalization. The normalized values are stored in column 'Concentration'
#' of aggregated data.
#'
#' @param metalyzer_se A MetAlyzer object
#' @param norm_method A character specifying the normalization method to use, which
#' should be one of 'TIC', 'VSN', or 'median'
#' @import dplyr, limma
data_normalization <- function(metalyzer_se, norm_method) {
  aggregated_data <- metalyzer_se@metadata$aggregated_data
  # Create temporary aggregated data for concatenating needed information later
  #### Do not know why Metabolite column is grouped by MetAlyzer
  tmp_aggregated_data <- dplyr::ungroup(aggregated_data) %>%
    dplyr::select(-Concentration)
  # Prepare data matrix for conducting normalization
  data_mat <- dplyr::ungroup(aggregated_data) %>%
    dplyr::select(Metabolite, ID, Concentration) %>%
    tidyr::pivot_wider(names_from = 'ID', values_from = 'Concentration') %>%
    tibble::column_to_rownames('Metabolite') %>%
    as.matrix()
  if (norm_method %in% 'TIC') {
    # Do TIC normalization
    row_sums <- rowSums(data_mat, na.rm = T)
    median_rowSums <- median(row_sums)
    norm_data <- apply(data_mat, 2, function(smp_conc) {
      smp_conc/row_sums * median_rowSums
    }) %>%
      glog2()
  } else if (norm_method %in% 'VSN') {
    # Do vsn normalization
    fit <- vsnMatrix(data_mat)
    norm_data <- predict(fit, data_mat)
  } else if (norm_method %in% 'median') {
    # Do median normalization that conducts median scaling on log2 transformed data
    # Use generalized log2 transformation to avoid -Inf
    log2_data <- glog2(data_mat)
    median_smpConc <- apply(log2_data, 2, function(smp_conc) {median(smp_conc, na.rm = T)})
    median_allVals <- median(log2_data, na.rm = T)
    norm_data <- sapply(seq_len(ncol(log2_data)), function(i) {
      log2_data[, i] - median_smpConc[i] + median_allVals
    })
    colnames(norm_data) <- colnames(data_mat)
    
    # Substantial negative values exist in log2 transformed data due to values smaller
    # than 1 in original data, which will cause problem in limma::normalizeBetweenArrays
    # that somehow also performs log transformation..
    #' Warning message:
    #' In log(apply(x, 2, median, na.rm = TRUE)) : NaNs produced
    # norm_data <- glog2(data_mat) %>%
    #   limma::normalizeBetweenArrays(method = 'scale')
  }
  # Convert matrix to aggregated long data
  aggregated_data <- tibble::as_tibble(norm_data, rownames = 'Metabolite') %>%
    tidyr::pivot_longer(cols = -'Metabolite',
                        names_to = 'ID',
                        values_to = 'Concentration') %>%
    #### Make character columns class factor as original aggregated data
    dplyr::mutate(Metabolite = factor(Metabolite, levels = levels(tmp_aggregated_data$Metabolite)),
                  ID = factor(ID, levels = levels(tmp_aggregated_data$ID))) %>%
    dplyr::left_join(tmp_aggregated_data, by = c('ID', 'Metabolite')) %>%
    dplyr::select(ID, Metabolite, Class, Concentration, Status)
  metalyzer_se@metadata$aggregated_data <- aggregated_data
  return(metalyzer_se)
}
#' @title Get the File Path of a File Inside the Data Folder
#' @description This function returns the file path for a specified file inside the `data` folder of a Shiny app.
#' It works both locally during development and when the app is deployed on shinyapps.io. The function first tries 
#' to access the file using a relative path. If the file is not found, it falls back to the working directory to ensure 
#' compatibility with the Shiny server environment.
#'
#' @param file_name A string representing the name of the file you want to access (e.g., "example.csv").
#' @return A string representing the full path to the specified file.
#' 
#' @details This function is particularly useful in Shiny apps where the file needs to be accessed from both
#' the local environment and the hosted environment on shinyapps.io. It ensures that the file can be correctly 
#' located regardless of where the app is running.
#' 
#' @examples
#' # Get path to the file "Pathway_120924.xlsx" in the data folder
#' file_path <- get_data_file_path("Pathway_120924.xlsx")
#' 
#' # Read the CSV file
#' data <- read.csv(file_path)
#' 
#' @export
get_data_file_path <- function(file_name) {
  data_path <- file.path("data", file_name)
  
  if (file.exists(data_path)) {
    return(data_path)
  } else {
    return(file.path(getwd(), data_path))
  }
}
