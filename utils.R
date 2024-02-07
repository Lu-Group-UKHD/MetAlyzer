library(ggplot2)
library(plotly)
library(dplyr)
library(MetAlyzer)
library(gridExtra)
library(viridis)
library(ggpubr)
library(viridisLite)

#' @title Calculate log2 fold change
#'
#' @description This function calculates the log2 fold change of two groups from
#' plotting_data.
#' @param metalyzer_se A Metalyzer object
#' @param categorical A column specifying the two groups
#'
#' @return A data frame containing the log2 fold change for each metabolite

calculate_log2FC <- function(metalyzer_se, categorical) {
  
  ## Create a new dataframe to calculate the log2FC
  if (!categorical %in% colnames(metalyzer_se@metadata$aggregated_data)) {
    metalyzer_se <- expand_aggregated_data(metalyzer_se,
                                           meta_data_column = categorical)
  }
  
  # Perform Log2 transformation
  aggregated_data <- metalyzer_se@metadata$aggregated_data
  aggregated_data <- mutate(aggregated_data,
                            log2_Conc = transform(.data$mod_Conc, base::log2), ### Column name might change!!!
                            .after = .data$mod_Conc) ### Column name might change!!!
  metalyzer_se@metadata$aggregated_data <- aggregated_data

  df <- metalyzer_se@metadata$aggregated_data %>%
    ungroup(all_of(categorical)) %>%
    mutate(Value = .data$log2_Conc,
        Group = !!sym(categorical)) %>%
    select(.data$Metabolite,
           .data$Class,
           .data$Group,
           .data$Value)

  ## Check if already factor
  group_vec <- df$Group
  if (!methods::is(group_vec, "factor")) {
    group_vec <- factor(group_vec)
    df$Group <- group_vec
    cat("Warning: No order was given for categorical!\n")
  } else {
    group_vec <- droplevels(group_vec)
  }
  group_levels <- levels(group_vec)
  if (length(group_levels) > 2) {
    cat("Warning: More than two levels were given! Dropping",
        paste(group_levels[3:length(group_levels)], collapse = ", "), "\n")
  }
  cat("Info: Calculating log2 fold change from ", group_levels[1], " to ",
      group_levels[2], " (column: ", categorical, ").\n", sep = "")
  df <- filter(df, .data$Group %in% group_levels[1:2])
  df$Group <- droplevels(df$Group)

  ## Check for further grouping
  grouping_vars <- as.character(groups(df))

  if (!"Metabolite" %in% grouping_vars) {
    grouping_vars[length(grouping_vars)+1] <- "Metabolite"
  }
  
  ## Calculate log2FC and p-val
  cat("Info: Calculating log2 fold change groupwise (",
      paste(grouping_vars, collapse = " * "),
      ") using a linear model...  ", sep = "")
  options(warn = -1)
  change_df <- df %>%
    group_modify(~ apply_linear_model(df = .x)) %>%
    ungroup(.data$Metabolite) %>%
    mutate(qval = qvalue::qvalue(.data$pval, pi0 = 1)$qvalues)
  options(warn = 0)
  cat("finished!\n")
  
  metalyzer_se@metadata$log2FC <- change_df
  return(metalyzer_se)
}

#' @title Calculate log2 fold change
#'
#' @description This function applies a linear model to calculate the log2 fold change and
#' its significance
#' @param df A subset data frame
#' @param ...

apply_linear_model <- function(df, ...) {
  class <- df$Class[1]
  df <- df %>%
    filter(!is.na(.data$Value)) %>%
    droplevels()
  if (length(levels(df$Group)) != 2) {
    l2fc <- NA
    pval <- NA
  } else {
    fit1 <- stats::lm(Value ~ Group, data = df)
    l2fc <- fit1$coefficients[2]
    fit_dim <- dim(summary(fit1)$coefficients)
    pval <- summary(fit1)$coefficients[fit_dim[1],fit_dim[2]]
  }
  output_df <- data.frame(Class = class,
                          log2FC = l2fc,
                          pval = pval,
                          row.names = NULL)
  return(output_df)
}

#' @title Transformation
#'
#' @description This function performs transformation of concentration values.
#'
#' @param vec a vector of concentration values
#' @param func A function for transformation
#'
#' @keywords internal
transform <- function(vec, func) {
  vec[vec > 0 & !is.na(vec)] <- func(vec[vec > 0 & !is.na(vec)])
  return(vec)
}


#' Plotly Log2FC Scatter Plot
#'
#' This function returns a list with interactive 
#' scatterplot based on log2 fold change data. 

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
  
  ## Data: Determine labels
  signif_p_data <- filter(p_data, .data$signif_color != names(signif_colors)[1])
  
  labels <- sapply(p_data$Metabolite, function(m) {
    m <- as.character(m)
    label <- if_else(m %in% signif_p_data$Metabolite, m, "")
    return(label)
  })
  
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
  
  ylims <- c(min(Log2FCTab$log2FC) - 0.75, max(Log2FCTab$log2FC) + 0,75)

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
  plotly_scatter <- ggplotly(p_scatter, tooltip = "text", showlegend = FALSE)
  plotly_scatter <- hide_legend(plotly_scatter)
  
  # Grab Legend ggplot
  scatter_legend <- ggpubr::get_legend(p_scatter)
  legend <- grid.arrange(scatter_legend, ncol=1)
  
  return(list(Plot = plotly_scatter, Legend = legend))
}

#' Plotly Log2FC Vulcano Plot
#'
#' This function returns a list with interactive 
#' vulcanoplot based on log2 fold change data. 

plotly_vulcano <- function(Log2FCTab, cutoff_y = 0.05, cutoff_x = 1.5) {
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

  # Data Vulcano: Create Dataframe for vulcano plot
  log2FCvulcano <- Log2FCTab
  log2FCvulcano$Class <- as.character(log2FCvulcano$Class)
  log2FCvulcano$Class[log2FCvulcano$qval > cutoff_y] <- NA
  log2FCvulcano$Class[abs(log2FCvulcano$log2FC) < log2(cutoff_x)] <- NA
  
  log2FCvulcano$labels <- as.character(log2FCvulcano$Metabolite)
  log2FCvulcano$labels[which(is.na(log2FCvulcano$Class))] <- ""
  
  if("highlight_metabolites" %in% colnames(log2FCvulcano)) {
    ## Plot: Create vulcano ggplot object with highlighted points
    p_fc_vulcano_highlighted <- ggplot(log2FCvulcano %>%
                                      arrange(desc(highlight_metabolites)),
                                    aes(x = .data$log2FC,
                                        y = -log10(.data$qval),
                                        color = .data$highlight_metabolites,
                                        label = labels)) +
      geom_point(size = 1, aes(text = paste0(Metabolite, "\nClass: ", Class, "\nlog2 Fold Change: ", round(log2FC, digits=5), "\np-value: ", round(pval, digits=5)))) +
      geom_vline(xintercept=c(-log2(cutoff_x), log2(cutoff_x)), col="black", linetype="dashed") +
      geom_hline(yintercept=-log10(cutoff_y), col="black", linetype="dashed") +
      scale_color_manual('',
                        breaks = c("Other Metabolites", "Highlighted Metabolite(s)"),
                        values = c("#d3d3d3","#56070C"),
                        drop = FALSE,
                        guide = guide_legend(override.aes = list(size = 2),
                                              order=2, ncol = 2)) +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            legend.key = element_rect(fill = 'white')) +
      labs(x = 'log2(FC)', y = "-log10(p)")
    
    ## Interactive: Create interactive plot
    p_vulcano_highlighted <- ggplotly(p_fc_vulcano_highlighted, tooltip = "text")

    return(p_vulcano_highlighted)
  } else {
    ## Plot: Create vulcano ggplot object
    p_fc_vulcano <- ggplot(log2FCvulcano,
                          aes(x = .data$log2FC,
                              y = -log10(.data$qval),
                              color = .data$Class,
                              label = labels)) +
      geom_vline(xintercept=c(-log2(cutoff_x), log2(cutoff_x)), col="black", linetype="dashed") +
      geom_hline(yintercept=-log10(cutoff_y), col="black", linetype="dashed") +
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
    p_vulcano <- ggplotly(p_fc_vulcano, tooltip = "text")

    return(p_vulcano)
  }
}

#' Plotly Log2FC Network Plot
#'
#' This function returns a list with interactive 
#' networkplot based on log2 fold change data. 

plotly_network <- function(Log2FCTab, q_value=0.05) {
  pathway_file <- MetAlyzer::pathway()
    
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
    if (length(m_vec) > 1) {
      # Nodes with more than 1 metabolite assigned
      tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
      if (nrow(tmp_df) > 0) {
        # At least one of the metabolites is significantly change
        # -> take the mean log2 fold change
        l2fc <- sum(tmp_df$log2FC) / length(m_vec)
      } else {
        if (any(tmp_df$Metabolite %in% levels(Log2FCTab$Metabolite))) {
          # At least one metabolite was measured but none are significantly changed
          l2fc <- 0
        } else {
          # None of the metabolites were measured
          l2fc <- NA
        }
      }
    } else {
      # Nodes with 0 or 1 metabolite assigned
      if (m_vec %in% signif_df$Metabolite) {
        # Metabolite is significantly changed
        l2fc <- signif_df$log2FC[which(signif_df$Metabolite == m_vec)]
      } else if (m_vec %in% levels(Log2FCTab$Metabolite)) {
        # Metabolite was measured but is not significantly changed
        l2fc <- 0
      } else {
        # Metabolite was not measured
        l2fc <- NA
      }
    }
    return(l2fc)
  })

  ## Add p-value to nodes_df
  nodes$p_value <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
    if (length(m_vec) > 1) {
      # Nodes with more than 1 metabolite assigned
      tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
      if (nrow(tmp_df) > 0) {
        # At least one of the metabolites is significantly change
        # -> How to calculate the combined p value???
        pval<- sum(tmp_df$pval) / length(m_vec)
      } else {
        if (any(tmp_df$Metabolite %in% levels(Log2FCTab$Metabolite))) {
          # At least one metabolite was measured but none are significantly changed
          pval <- 0
        } else {
          # None of the metabolites were measured
          pval <- NA
        }
      }
    } else {
      # Nodes with 0 or 1 metabolite assigned
      if (m_vec %in% signif_df$Metabolite) {
        # Metabolite is significantly changed
        pval <- signif_df$pval[which(signif_df$Metabolite == m_vec)]
      } else if (m_vec %in% levels(Log2FCTab$Metabolite)) {
        # Metabolite was measured but is not significantly changed
        pval <- 0
      } else {
        # Metabolite was not measured
        pval <- NA
      }
    }
    return(pval)
  })
  
  ## Draw network
  # Create a plot of the network using ggplotly
  label_size <- 20
  area_size <- 10
  edge_size <- 1.25
  annotation_size <- 11

  # Create the coloured area behind the edges
  # Turn edges with no pathway into white
  area_shapes <- list()
  for (i in seq_along(edges$Color)) {
    if (is.na(edges$Color[i])) {
      edges$Color[i] <- "white"
    }
  }

  for (i in 1:nrow(edges)) {
    area_shape <- list(
      type = "line",
      x0 = edges$x_start[i],
      y0 = edges$y_start[i],
      x1 = edges$x_end[i],
      y1 = edges$y_end[i],
      line = list(
        color = edges$Color[i],
        width = area_size
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
        width = edge_size
      )
    )
    edge_shapes[[i]] <- edge_shape
  }
  edges_area_combined <- c(area_shapes, edge_shapes)

  # Function to convert numerical values to color codes
  value_to_color <- function(value) {
    if (is.na(value)) {
      return("grey")
    } else {
      # Using the viridis color scale, adjust 'option' based on your preference
      color_scale <- viridis(10, option = "D")
      nodes_range <- na.omit(nodes$FC_thresh)

      color_index <- findInterval(value, seq(min(nodes_range), max(nodes_range)+0.1, length.out = length(color_scale) + 1))
      return(color_scale[color_index])
    }
  }

  # Adding a new column 'color' with color codes based on the numerical values
  nodes$color <- sapply(nodes$FC_thresh, value_to_color)

  # Remove axis and add title
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

  # Create the nodes
  network <- plot_ly(nodes,
                    x = nodes$x,
                    y = nodes$y,
                    type = "scatter",
                    mode = "markers",
                    marker = list(color = nodes$FC_thresh, 
                            colorbar = list(title = "log2FC with FDR correction"),
                            colorscale='Viridis',
                            showscale = TRUE))

  p_network <- layout(
    network,
    title = '',
    shapes = edges_area_combined,
    xaxis = axis,
    yaxis = axis,
    hovermode = FALSE) 

  # Add annotations for the nodes
  for (i in 1:nrow(nodes)) {
    p_network <- p_network %>% add_annotations(
      text = nodes$Label[i],
      x = nodes$x[i],
      y = nodes$y[i],
      arrowhead = 0,
      font = list(size = annotation_size, color = "white"),
      ax = 0,
      ay = 0,
      bgcolor = nodes$color[i],
      opacity = 1,
      hovertext = paste0("log2 Fold Change: ", round(nodes$FC_thresh[i], 5),
                        "\nPathway: ", nodes$Pathway[i],
                        "\np value: ", round(nodes$p_value[i], 5))
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
      font = list(size = label_size, color = pathways$Color[i])
    )
  }

  return(p_network)
}

#' @title Read Named Regions
#'
#' @description This function reads in the named regions of an excel file.
#'
#' @param file_path The file path of the file
#' @param named_region The region name u want to read in

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