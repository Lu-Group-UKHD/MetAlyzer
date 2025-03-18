#' Plot Pathway Network
#'
#' @description This function plots the log2 fold change for each metabolite and visualizes it, in a pathway network.
#'
#' @param metalyzer_se A Metalyzer object
#' @param q_value The q-value threshold for significance
#' @param metabolite_text_size The text size of metabolite labels
#' @param connection_width The line width of connections between metabolites
#' @param pathway_text_size The text size of pathway annotations
#' @param pathway_width The line width of pathway-specific connection coloring
#' @param plot_column_name Column name in the Log2FC dataframe to plot; 
#' @param exlude_pathways Pathway names that are exluded from plotting
#' @return ggplot object
#' 
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import SummarizedExperiment
#' @importFrom rlang .data
#' @export
#' 
#' @examples
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_mutation_data_xl())
#' metalyzer_se <- filterMetabolites(
#'   metalyzer_se,
#'   drop_metabolites = "Metabolism Indicators"
#' )
#' metalyzer_se <- renameMetaData(
#'   metalyzer_se,
#'   Mutant_Control = "Sample Description"
#' )
#' 
#' metalyzer_se <- calculate_log2FC(
#'   metalyzer_se,
#'   categorical = "Mutant_Control",
#'   impute_perc_of_min = 0.2,
#'   impute_NA = FALSE
#' )
#' 
#' network <- plot_network(metalyzer_se, q_value = 0.05)

plot_network <- function(
    metalyzer_se,
    q_value=0.05,
    metabolite_text_size=3,
    connection_width=0.75,
    pathway_text_size=6,
    pathway_width=3,
    plot_column_name="log2FC",
    exclude_pathways=NULL 
  ) {
  log2FC_df <- metalyzer_se@metadata$log2FC
  
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
  pathways <- filter(pathways, !Label %in% exclude_pathways)
  excluded_labels <- nodes$Label[nodes$Pathway %in% exclude_pathways]
  edges <- edges %>%
    filter(!(Node1 %in% excluded_labels | Node2 %in% excluded_labels))
  
  nodes <- filter(nodes, !Pathway %in% exclude_pathways)

  ## Add log2FC to nodes_df
  signif_df <- filter(log2FC_df,
                      !is.na(.data$log2FC),
                      !is.na(.data$qval),
                      !is.na(.data$pval),
                      .data$qval <= q_value)

  
  nodes$log2FC <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
    tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
    if (nrow(tmp_df) > 0) {
      # Alteast 1 significantly changed
      l2fc <- sum(tmp_df$log2FC) / nrow(tmp_df)
    } else if (nrow(tmp_df) == 0 && any(m_vec %in% log2FC_df$Metabolite)) {
      # Not significantly changed but measured
      l2fc <- sum(log2FC_df$log2FC[which(log2FC_df$Metabolite %in% m_vec)]) / length(m_vec)
    } else {
      # Not measured
      l2fc <- NA
    }
    return(l2fc)
  })
  
  ## Add p-value to nodes_df
  nodes$qval <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
    tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
    if (nrow(tmp_df) > 0) {
      # Alteast 1 significantly changed
      qval <- sum(tmp_df$qval) / nrow(tmp_df)
    } else if (nrow(tmp_df) == 0 && any(m_vec %in% log2FC_df$Metabolite)) {
      # Not significantly changed but measured
      qval <- sum(log2FC_df$qval[which(log2FC_df$Metabolite %in% m_vec)]) / length(m_vec)
    } else {
      # Not measured
      qval <- NA
    }
    return(qval)
  })

  ## Add p-value to nodes_df
  nodes$pval <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
    tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
    if (nrow(tmp_df) > 0) {
      # Alteast 1 significantly changed
      pval <- sum(tmp_df$pval) / nrow(tmp_df)
    } else if (nrow(tmp_df) == 0 && any(m_vec %in% log2FC_df$Metabolite)) {
      # Not significantly changed but measured
      pval <- sum(log2FC_df$pval[which(log2FC_df$Metabolite %in% m_vec)]) / length(m_vec)
    } else {
      # Not measured
      pval <- NA
    }
    return(pval)
  })

  nodes$add_col <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
    tmp_df <- filter(log2FC_df, .data$Metabolite %in% m_vec)
    if (nrow(tmp_df) > 0) {
      # Alteast 1 value per Node
      value <- sum(tmp_df[plot_column_name]) / nrow(tmp_df)
    } else {
      # Not measured
      value <- NA
    }
    return(value)
  })

  # In case the plotted column is log2FC or qval, as they are calculated differently
  if (plot_column_name == "log2FC") {
    nodes$add_col = nodes$log2FC
  } else if (plot_column_name == "qval") {
    nodes$add_col = nodes$qval
  } else if (plot_column_name == "pval") {
    nodes$add_col = nodes$pval
  }

  ## Draw network

  # Create a plot of the network using ggplot2 and ggrepel
  network <- ggplot()
  for (radius in unique(edges$Radius)) {
    rad_edges <- filter(edges, .data$Radius == radius)
    pathway_edges <- filter(rad_edges, !is.na(.data$Color))

    network <- network +
      # Add the round area behind the edges
      geom_curve(
        data = pathway_edges,
        aes(
          x = .data$x_start,
          y = .data$y_start,
          xend = .data$x_end,
          yend = .data$y_end,
          color = .data$Color
        ),
        linewidth = pathway_width,
        alpha = 0.3,
        curvature = radius,
        show.legend = FALSE
      ) +
      # Or add the edges as curved lines
      geom_curve(
        data = rad_edges,
        aes(
          x = nodes[.data$Node1, "x"],
          y = nodes[.data$Node1, "y"],
          xend = nodes[.data$Node2, "x"],
          yend = nodes[.data$Node2, "y"]
        ),
        color = "grey",
        linewidth = connection_width,
        curvature = radius
      )
  }
  network <- network +
    # Add labels at the position of the nodes
    geom_label(
      data = nodes,
      aes(
        x = .data$x,
        y = .data$y,
        label = .data$Label,
        fill = .data$add_col
      ),
      size = metabolite_text_size,
      color = "white"
    ) +
    scale_fill_viridis_c(option = "D", name = paste0(plot_column_name)) +

    # Add annotations
    geom_text(
      data = pathways,
      aes(
        x = .data$x,
        y = .data$y,
        label = .data$Label,
        color = .data$Color
      ),
      size = pathway_text_size,
      show.legend = FALSE
    ) +
    # # Set the x and y axis limits
    # xlim(0, 10) +
    # ylim(0, 10) +
    theme_void() +
    # Add a title and remove the x and y axis labels
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5))
  network

  #ggsave("network.pdf", network, width = 15, height = 10, bg = "white")
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
