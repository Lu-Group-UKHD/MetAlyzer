#' Plot Pathway Network
#'
#' @description This function plots the log2 fold change for each metabolite and visualizes it, in a pathway network.
#'
#' @param log2fc_df A dataframe with log2FC, qval, additional columns
#' @param q_value The q-value threshold for significance
#' @param values_col_name Column name of a column that holds numeric values, to be plotted \strong{Default = "log2FC"}
#' @param stat_col_name Columnname that holds numeric stat values that are used for significance \strong{Default = "qval"}
#' @param metabolite_col_name Columnname that holds the Metabolites
#' @param metabolite_text_size The text size of metabolite labels
#' @param connection_width The line width of connections between metabolites
#' @param pathway_text_size The text size of pathway annotations
#' @param pathway_width The line width of pathway-specific connection coloring
#' @param exclude_pathways Pathway names that are exluded from plotting
#' @param color_scale A string specifying the color scale to use. Options include `"Viridis"`, `"Plasma"`, `"Magma"`, `"Inferno"`, `"Cividis"`, `"Rocket"`, `"Mako"`, and `"Turbo"`, which use the `viridis` color scales. If `"gradient"` is selected, a custom gradient is applied based on `gradient_colors`.
#' @param gradient_colors A vector of length 2 or 3 specifying the colors for a custom gradient. If two colors are provided (`c(low, high)`), `scale_fill_gradient()` is used. If three colors are provided (`c(low, mid, high)`), `scale_fill_gradient2()` is used. If `NULL` or incorrectly specified, the viridis color scale is applied.
#' @param save_as \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "NULL"}
#' @param folder_name Name of the folder where the plot will be saved. Special characters will be removed automatically. \strong{Default = date}
#' @param folder_path \emph{Optional: } User-defined path where the folder should be created. 
#' If not provided, results will be saved in `MetAlyzer_results` within the working directory. \strong{Default = NULL}
#' @param file_name Name of the output file (without extension). \strong{Default = "network"}
#' @param format File format for saving the plot (e.g., "png", "pdf", "svg"). \strong{Default = "pdf"}
#' @param width Width of the saved plot in specified units. \strong{Default = 29.7}
#' @param height Height of the saved plot in specified units. \strong{Default = 21.0}
#' @param units Units for width and height (e.g., "in", "cm", "mm"). \strong{Default = "cm"}
#' @param overwrite Logical: If `TRUE`, overwrite existing files without asking. If `FALSE`, prompt user before overwriting. \strong{Default = FALSE}
#' @return list with ggplot object and table of node summaries
#' 
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import viridis
#' @import viridisLite
#' @importFrom rlang .data !!
#' @importFrom magrittr %>%
#' @export
#' 
#' @examples
#' log2fc_df <- readRDS(MetAlyzer::toy_diffres())
#' network <- MetAlyzer::plot_network(log2fc_df, q_value = 0.05)
#' network$Plot
#' network$Table

plot_network <- function(log2fc_df,
                         q_value = 0.05,
                         metabolite_col_name = "Metabolite",
                         values_col_name = "log2FC",
                         stat_col_name = "qval",
                         metabolite_text_size = 3,
                         connection_width = 0.75,
                         pathway_text_size = 6,
                         pathway_width = 3,
                         exclude_pathways = NULL,
                         color_scale = "Viridis",
                         gradient_colors = NULL,
                         save_as = NULL,
                         folder_name = format(Sys.Date(), "%Y-%m-%d"),
                         folder_path = NULL,
                         file_name = "network",
                         format = "pdf",
                         width = 29.7,
                         height = 21.0,
                         units = "cm",
                         overwrite = FALSE) {
  ### Checks
  if (!(metabolite_col_name %in% colnames(log2fc_df))) {
    stop(paste0("Column '", metabolite_col_name, "' is missing. Please select one of: ", paste(colnames(log2fc_df), collapse = ", ")))
  }

  if (!(values_col_name %in% colnames(log2fc_df))) {
    stop(paste0("Column '", values_col_name, "' is missing. Please select one of: ", paste(colnames(log2fc_df), collapse = ", ")))
  }

  if (!(stat_col_name %in% colnames(log2fc_df))) {
    stop(paste0("Column '", stat_col_name, "' is missing. Please select one of: ", paste(colnames(log2fc_df), collapse = ", ")))
  }

  if (!is.numeric(log2fc_df[[values_col_name]])) {
    stop(paste0("Column '", values_col_name, "' must contain numerical values."))
  }

  if (!is.numeric(log2fc_df[[stat_col_name]])) {
    stop(paste0("Column '", stat_col_name, "' must contain numerical values."))
  }
  
  network_file <- MetAlyzer::pathway()

  ### Read in Excel file
  pathways <- MetAlyzer:::read_pathways(network_file)
  nodes <- MetAlyzer:::read_nodes(network_file, pathways)
  edges <- MetAlyzer:::read_edges(network_file, nodes, pathways)
  
  pathways <- dplyr::filter(pathways, !Pathway %in% exclude_pathways)
  nodes <- dplyr::filter(nodes, !Pathway %in% exclude_pathways)
  edges <- dplyr::filter(edges, Node1 %in% nodes$Label | Node2 %in% nodes$Label)
  
  nodes_separated <- tidyr::separate_rows(nodes, Metabolites, sep = "\\s*;\\s*")
  
  nodes_joined <- dplyr::left_join(nodes_separated, log2fc_df, by = c("Metabolites" = metabolite_col_name))

  updated_nodes_list <- MetAlyzer:::calculate_node_aggregates_conditional(nodes_sep_df = nodes_joined, nodes_orig_df = nodes, q_value = q_value, stat_col_name = stat_col_name, c(values_col_name))

  # --- Create the dataframe for excel export ---
  nodes_separated_processed <- updated_nodes_list$nodes_separated

  nodes_separated_shortend <- nodes_separated_processed %>%
    dplyr::filter(!is.na(values_col_name)) 
  
  summary_log2fc <- nodes_separated_shortend %>%
    dplyr::group_by(Label) %>%
    dplyr::summarise(
      values_collapsed = paste(.data[[values_col_name]], collapse = "; "),
      stat_collapsed = paste(.data[[stat_col_name]], collapse = "; "),
      .groups = 'drop'
    )

  cols_to_summarise_unique <- setdiff(names(nodes_separated_shortend), c("Label", values_col_name, stat_col_name))

  summary_others <- nodes_separated_shortend %>%
    dplyr::group_by(.data$Label) %>%
    dplyr::summarise(
      collapsed_count = dplyr::n(),
      dplyr::across(
        .cols = all_of(cols_to_summarise_unique), # Use the identified columns
        .fns = ~ paste(., collapse = "; ")
      ),
      .groups = 'drop'
    )

  nodes_collapsed <- left_join(summary_others, summary_log2fc, by = "Label") %>%
    dplyr::rename(values = values_collapsed, stat = stat_collapsed,) %>%
    dplyr::mutate(Pathway = if_else(Pathway == "", NA_character_, Pathway))
    
  # --- The dataframe for plotting ---
  nodes_original_processed <- updated_nodes_list$nodes

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
      data = nodes_original_processed,
      aes(
        x = .data$x,
        y = .data$y,
        label = .data$Label,
        fill = .data[[values_col_name]]
      ),
      size = metabolite_text_size,
      color = "white"
    ) +
    switch(color_scale,
           "gradient" = if (!is.null(gradient_colors) && length(gradient_colors) %in% c(2, 3)) {
             if (length(gradient_colors) == 2) {
               scale_fill_gradient(low = gradient_colors[1], high = gradient_colors[2], name = values_col_name)
             } else {
               scale_fill_gradient2(low = gradient_colors[1], mid = gradient_colors[2], high = gradient_colors[3], 
                                    midpoint = 0, name = values_col_name)
             }
           } else {
             message("gradient_colors is NULL or incorrectly specified. Falling back to viridis scale.")
             scale_fill_viridis(option = "D", name = values_col_name)  # default fallback
           },
           scale_fill_viridis(option = MetAlyzer:::create_viridis_style(color_scale, type = "initial"), name = values_col_name)  # default to viridis
    ) +

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


    MetAlyzer:::save_plot(network,
              folder_name = folder_name,
              folder_path = folder_path,
              file_name = file_name,
              width = width,
              height = height,
              units = units,
              format = save_as,
              overwrite = overwrite)
  
  return(list("Plot" = network, "Table" = nodes_collapsed))

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

# --- Sub-Functions ---

#' Read and Validate Pathway Annotations
#'
#' Reads pathway data from a specified named region in the pathway file,
#' validates entries, removes invalid ones, and sets row names.
#'
#' @param network_file Path to the input file containing pathway data.
#' @param region_name The named region or sheet containing pathway header info.
#' @return A data frame of validated pathway annotations with labels as row names.
read_pathways <- function(network_file, region_name = "Pathways_Header") {
  # Assuming MetAlyzer::pathway() provides the file path
  # and read_named_region is available in the environment
  pathways <- read_named_region(network_file, region_name)
  
  invalid_annotations <- which(
    is.na(pathways$Label) |
    duplicated(pathways$Label) |
    is.na(pathways$x) |
    is.na(pathways$y) |
    is.na(pathways$Color)
  )
  if (length(invalid_annotations) > 0) {
    cat("Warning: Removing", length(invalid_annotations), "invalid pathways.\n")
    pathways <- pathways[-invalid_annotations, ]
  }
  if (nrow(pathways) > 0) {
     rownames(pathways) <- pathways$Label
  } else {
     cat("Warning: No valid pathways found.\n")
     # Return an empty data frame with expected columns if needed downstream
     # return(data.frame(Label=character(), x=numeric(), y=numeric(), Color=character(), stringsAsFactors=FALSE))
  }
  return(pathways)
}

#' Read and Validate Network Nodes (Metabolites)
#'
#' Reads node data from a specified named region, validates entries against
#' pathway information, cleans labels, removes invalid nodes, and sets row names.
#'
#' @param network_file Path to the input file containing node data.
#' @param pathways A data frame of validated pathways (output of read_pathways).
#' @param region_name The named region or sheet containing metabolite header info.
#' @return A data frame of validated nodes with labels as row names.
read_nodes <- function(network_file, pathways, region_name = "Metabolites_Header") {
  nodes <- read_named_region(network_file, region_name)
  nodes$Pathway[is.na(nodes$Pathway)] <- ""

  valid_pathway_names <- character(0) # Initialize empty vector
  if (nrow(pathways) > 0) {
      valid_pathway_names <- c(rownames(pathways), "")
  } else {
      valid_pathway_names <- "" # Only allow unassigned pathway if no pathways exist
  }

  invalid_nodes <- which(
    is.na(nodes$Label) |
    duplicated(nodes$Label) |
    is.na(nodes$x) |
    is.na(nodes$y) |
    !nodes$Pathway %in% valid_pathway_names
  )

  if (length(invalid_nodes) > 0) {
    cat("Warning: Removing", length(invalid_nodes), "invalid nodes.\n")
    nodes <- nodes[-invalid_nodes, ]
  }

  if (nrow(nodes) > 0) {
    # Use the original label for rownames before cleaning
    rownames(nodes) <- nodes$Label
    # Remove #suffix from labels after setting rownames
    nodes$Label <- gsub("#[0-9]+$", "", nodes$Label)
  } else {
     cat("Warning: No valid nodes found.\n")
     # Return an empty data frame with expected columns if needed downstream
     # return(data.frame(Label=character(), x=numeric(), y=numeric(), Pathway=character(), stringsAsFactors=FALSE))
  }
  return(nodes)
}

#' Read and Validate Network Edges (Connections)
#'
#' Reads edge data from a specified named region, validates that connected
#' nodes exist and are not self-loops, and removes invalid edges.
#'
#' @param network_file Path to the input file containing edge data.
#' @param nodes A data frame of validated nodes (output of read_nodes).
#' @param pathways A data frame of validated pathways (output of read_pathways).
#' @param region_name The named region or sheet containing connections header info.
#' @return A data frame of validated edges.
read_edges <- function(network_file, nodes, pathways, region_name = "Connections_Header") {
  edges <- read_named_region(network_file, region_name)
  
  valid_node_names <- character(0) # Initialize empty vector
   if (nrow(nodes) > 0) {
       valid_node_names <- rownames(nodes)
   } else {
       # If there are no valid nodes, all edges are invalid
       cat("Warning: No valid nodes exist, removing all connections.\n")
       return(edges[0, ]) # Return empty dataframe with same columns
   }

  invalid_edges <- which(
    is.na(edges$Node1) | # Check for NA node names first
    is.na(edges$Node2) |
    !edges$Node1 %in% valid_node_names |
    !edges$Node2 %in% valid_node_names |
    edges$Node1 == edges$Node2
  )

  if (length(invalid_edges) > 0) {
    cat("Warning: Removing", length(invalid_edges), "invalid connections.\n")
    edges <- edges[-invalid_edges, ]
  }
   if (nrow(edges) == 0) {
     cat("Info: No valid edges found after validation.\n")
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
  return(edges)
}

#' Creates a viridis color style for Plotly plots.
#'
#' This function can generate either a Plotly-compatible colorscale for a
#' color bar or a vector of hex color codes for manual coloring.
#'
#' @param color_scale The name of the palette (e.g., "Magma", "Viridis").
#' @param type The desired output type: "scale" (for a color bar), "hex" or "inital"
#'   (for a color scalde, a vector of hex codes, or the correct scale for viridis package). Defaults to "scale".
#' @param data The data frame containing the values. Only required if type = "hex".
#' @param values_col_name The name of the column with numeric values. Only
#'   required if type = "hex".
#'
#' @return A data frame if type is "scale", or a character vector if type is "hex".
create_viridis_style <- function(color_scale,
                                 type = "scale",
                                 data = NULL,
                                 values_col_name = NULL) {
  option <- switch(color_scale,
                   "Magma"   = "A",
                   "Inferno" = "B",
                   "Plasma"  = "C",
                   "Viridis" = "D",
                   "Cividis" = "E",
                   "Rocket"  = "F",
                   "Mako"    = "G",
                   "Turbo"   = "H",
                   "D"
  )
  if (type == "scale") {
    n_colors <- 11 # A small number of steps is efficient for a Plotly scale
    palette_colors <- viridis(n_colors, option = option)
    stop_points <- seq(0, 1, length.out = n_colors)
    return(data.frame(stop = stop_points, color = palette_colors))

  } else if (type == "hex") {
    if (is.null(data) || is.null(values_col_name)) {
      stop("For type = 'hex', you must provide 'data' and 'values_col_name'.")
    }

    n_colors <- 256
    palette <- viridis(n_colors, option = option)

    values <- data[[values_col_name]]
    valid_values <- na.omit(as.numeric(values))

    if (length(valid_values) == 0) return(rep("grey", length(values)))

    value_range <- range(valid_values)

    if (diff(value_range) == 0) {
      middle_color <- palette[n_colors / 2]
      return(ifelse(is.na(values), "grey", middle_color))
    }

    breaks <- seq(value_range[1], value_range[2], length.out = n_colors + 1)

    return(sapply(values, function(val) {
      if (is.na(val)) {
        "grey"
      } else {
        color_index <- findInterval(val, breaks, all.inside = TRUE)
        palette[color_index]
      }
    }))
  } else if (type == "initial") {
    return(option)
  } else {
    stop("Invalid 'type' specified. Please choose 'scale' or 'hex'.")
  }
}

#' Calculate Node-Level Aggregate Statistics (Conditional on Significance)
#'
#' Calculates mean log2FC, p-value, and q-value for each node (Label),
#' prioritizing significant metabolites (qval <= q_value). If none are significant,
#' uses all measured metabolites for the node. Adds results to both dataframes.
#'
#' @param nodes_sep_df Dataframe with metabolites separated (e.g., 'nodes_final').
#'   Must contain Label, log2FC, qval.
#' @param nodes_orig_df Original dataframe with potentially semi-colon separated metabolites.
#'   Must contain Label.
#' @param q_value Significance threshold for q-values (e.g., 0.05).
#' @param stat_col_name p value column name
#'
#' @import dplyr
#' @importFrom rlang .data 
#' @importFrom magrittr %>%
#' 
#' @return A list containing two dataframes:
#'   $nodes_separated: Input nodes_sep_df with 2 new columns:
#'     node_values, node_stat
#'   $nodes: Input nodes_orig_df with 2 new columns:
#'     node_values, node_stat
#'
calculate_node_aggregates_conditional <- function(nodes_sep_df,
                                                  nodes_orig_df,
                                                  q_value,
                                                  stat_col_name,
                                                  ...) {
  # --- Capture the dynamic column names ---
  value_col_names <- c(...)

  # --- Input Validation ---
  if (!"Label" %in% names(nodes_sep_df) || !"Label" %in% names(nodes_orig_df)) {
    stop("Both dataframes must contain a 'Label' column for grouping.")
  }
  # Combine all columns that need processing and check if they exist
  all_cols_to_process <- c(value_col_names, stat_col_name)
  if (length(all_cols_to_process) == 0) {
    stop("Please provide at least one column name to be processed.")
  }
  if (!all(all_cols_to_process %in% names(nodes_sep_df))) {
    stop("nodes_sep_df must contain all specified columns: ", paste(all_cols_to_process, collapse=", "))
  }
  if (missing(q_value) || !is.numeric(q_value) || length(q_value) != 1) {
    stop("Please provide a single numeric value for q_value threshold.")
  }

  # --- 1. Calculate Node-Level Aggregates using `across()` ---

  node_summary <- nodes_sep_df %>%
    dplyr::group_by(.data$Label) %>%
    # Use across() to apply the same function to multiple columns
    dplyr::summarise(
      dplyr::across(
        # Target all specified value columns and the stat column
        dplyr::all_of(all_cols_to_process),
        # Apply the calculation to each column (.x)
        ~ MetAlyzer:::calculate_conditional_mean(.x, .data[[stat_col_name]], q_value)
      ),
      .groups = "drop"
    )

  # --- 2. Update nodes_separated Dataframe ---

  # Create a prefixed version of the summary for joining, to avoid name clashes
  node_summary_prefixed <- node_summary %>%
    dplyr::rename_with(~ paste0("node_", .), .cols = -Label)

  nodes_sep_updated <- nodes_sep_df %>%
    dplyr::left_join(node_summary_prefixed, by = "Label")

  # --- 3. Update Original Nodes Dataframe ---
  
  nodes_orig_updated <- nodes_orig_df %>%
      # Remove original columns to prevent name clashes before joining the new aggregated values
      dplyr::select(-any_of(all_cols_to_process)) %>%
      # Join the new summary data; columns are already correctly named
      dplyr::left_join(node_summary, by = "Label")

  # --- 4. Return Updated Dataframes ---
  return(list(nodes_separated = nodes_sep_updated, nodes = nodes_orig_updated))
}
#' Helper function to calculate conditional mean based on significance
#'
#' Calculates mean of a metric vector based on q-values. Prioritizes values
#' where qval <= q_thresh. If none exist, uses all non-NA values.
#'
#' @param metric_vec Numeric vector of the metric to average (e.g., log2FC).
#' @param qval_vec Numeric vector of corresponding q-values.
#' @param q_thresh Significance threshold for q-values.
#' @return The calculated conditional mean, or NA.
calculate_conditional_mean <- function(metric_vec, qval_vec, q_thresh) {

  # Validate inputs are reasonable (basic check)
  if(length(metric_vec) != length(qval_vec)) {
    stop("metric_vec and qval_vec must have the same length.")
  }

  # Identify rows that are BOTH significant AND measured for the specific metric
  is_significant_and_measured <- !is.na(qval_vec) & qval_vec <= q_thresh & !is.na(metric_vec)

  # Identify rows that are measured for the metric (regardless of significance)
  is_measured <- !is.na(metric_vec)

  if (any(is_significant_and_measured)) {
    # Case 1: Significant & measured values exist -> Use mean of these values
    mean_val <- mean(metric_vec[is_significant_and_measured], na.rm = FALSE) # NAs excluded by definition
  } else if (any(is_measured)) {
    # Case 2: No significant values, but measured values exist -> Use mean of all measured values
    mean_val <- mean(metric_vec[is_measured], na.rm = FALSE) # NAs excluded by definition
  } else {
    # Case 3: No measured values exist for this metric in the group
    mean_val <- NA_real_
  }
  return(mean_val)
}
