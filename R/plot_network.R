#' Plot Pathway Network
#'
#' @description This function plots the log2 fold change for each metabolite and visualizes it, in a pathway network.
#'
#' @param log2fc_df A dataframe with log2FC, pval, qval, additional columns
#' @param q_value The q-value threshold for significance
#' @param values_col_name Column name of a column that holds numeric values, to be plotted \strong{Default = "log2FC"}
#' @param stat_col_name Columnname that holds numeric stat values that are used for significance \strong{Default = "qval"}
#' @param metabolite_col_name Columnname that holds the Metabolites
#' @param metabolite_text_size The text size of metabolite labels
#' @param connection_width The line width of connections between metabolites
#' @param pathway_text_size The text size of pathway annotations
#' @param pathway_width The line width of pathway-specific connection coloring
#' @param exlude_pathways Pathway names that are exluded from plotting
#' @param color_scale A string specifying the color scale to use. Options include `"viridis"`, `"plasma"`, `"magma"`, `"inferno"`, `"cividis"`, `"rocket"`, `"mako"`, and `"turbo"`, which use the `viridis` color scales. If `"gradient"` is selected, a custom gradient is applied based on `gradient_colors`.
#' @param gradient_colors A vector of length 2 or 3 specifying the colors for a custom gradient. If two colors are provided (`c(low, high)`), `scale_fill_gradient()` is used. If three colors are provided (`c(low, mid, high)`), `scale_fill_gradient2()` is used. If `NULL` or incorrectly specified, the viridis color scale is applied.
#' @param save_as \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "NULL"}
#' @param folder_name Name of the folder where the plot will be saved. Special characters will be removed automatically. \strong{Default = date}
#' @param folder_path \emph{Optional: } User-defined path where the folder should be created. 
#' If not provided, results will be saved in `MetAlyzer_results` within the working directory. \strong{Default = NULL}
#' @param file_name Name of the output file (without extension). \strong{Default = "network"}
#' @param width Width of the saved plot in specified units. \strong{Default = 29.7}
#' @param height Height of the saved plot in specified units. \strong{Default = 21.0}
#' @param units Units for width and height (e.g., "in", "cm", "mm"). \strong{Default = "cm"}
#' @param overwrite Logical: If `TRUE`, overwrite existing files without asking. If `FALSE`, prompt user before overwriting. \strong{Default = FALSE}
#' @return ggplot object
#' 
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import viridis
#' @import viridisLite
#' @importFrom rlang .data !!
#' @export
#' 
#' @examples
#' log2fc_df <- readRDS(toy_diffres())
#' network <- plot_network(log2fc_df, q_value = 0.05)

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
                         color_scale = "viridis",
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
  
  network_file <- MetAlyzer::pathway()

  ### Read in Excel file
  pathways <- MetAlyzer:::read_pathways(network_file)
  nodes <- MetAlyzer:::read_nodes(network_file, pathways)
  edges <- MetAlyzer:::read_edges(network_file, nodes, pathways)
  
  nodes <- dplyr::filter(nodes, !Pathway %in% exclude_pathways)
  
  nodes_separated <- tidyr::separate_rows(nodes, Metabolites, sep = "\\s*;\\s*")
  
  nodes_joined <- dplyr::left_join(nodes_separated, log2fc_df, by = c("Metabolites" = metabolite_col_name))

  updated_nodes_list <- MetAlyzer:::calculate_node_aggregates_conditional(nodes_joined, nodes, q_value, values_col_name, stat_col_name)

  # --- Create the dataframe for excel export ---
  nodes_separated_processed <- updated_nodes_list$nodes_separated

  nodes_separated_shortend <- nodes_separated_processed %>%
    dplyr::filter(!is.na(values_col_name)) %>%  
    dplyr::select(-c(Class, x, y, Shape))  
  
  summary_log2fc <- nodes_separated_shortend %>%
    dplyr::group_by(Label) %>%
    dplyr::summarise(
      values_collapsed = paste(values_col_name, collapse = "; "),
      stat_collapsed = paste(stat_col_name, collapse = "; "),
      .groups = 'drop'
    )

  cols_to_summarise_unique <- setdiff(names(nodes_separated_shortend), c("Label", values_col_name, stat_col_name))

  summary_others <- nodes_separated_shortend %>%
    dplyr::group_by(Label) %>%
    dplyr::summarise(
      collapsed_count = dplyr::n(),
      dplyr::across(
        .cols = all_of(cols_to_summarise_unique), # Use the identified columns
        .fns = ~ paste(unique(.), collapse = "; ")
      ),
      .groups = 'drop'
    )

  nodes_collapsed <- left_join(summary_others, summary_log2fc, by = "Label") %>%
    rename(values = values_collapsed, stat = stat_collapsed,) %>%
    mutate(Pathway = if_else(Pathway == "", NA_character_, Pathway))
    
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
           scale_fill_viridis(option = get_color_option(color_scale), name = values_col_name)  # default to viridis
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


    save_plot(network,
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

#' Get the color option based on the color scale
#'
#' This function maps a color scale name to a corresponding letter option.
#' The function returns a letter representing the color scale from a predefined mapping.
#' If the color scale provided is not recognized, the function fallsback to the viridis scale.
#'
#' @param color_scale A character string representing the name of the color scale.
#'        It should be one of: "magma", "inferno", "plasma", "viridis", 
#'        "cividis", "rocket", "mako", or "turbo".
#'
#' @return A character string representing the color option corresponding to the 
#'         provided color scale. If the color scale is not recognized, it returns `NA`.
#'
get_color_option <- function(color_scale) {
  color_map <- list(
    magma = "A",
    inferno = "B",
    plasma = "C",
    viridis = "D",
    cividis = "E",
    rocket = "F",
    mako = "G",
    turbo = "H"
  )
  
  return(ifelse(color_scale %in% names(color_map), color_map[[color_scale]], "D"))
}


#' save_plot is a helper function to save plots
#'
#' @description This function saves a given ggplot object to a specified folder and file format.
#' It ensures that the folder structure exists and cleans the folder name to remove special characters.
#'
#' @param plot A ggplot object to be saved.
#' @param folder_name Name of the folder where the plot will be saved. Special characters will be removed automatically. \strong{Default = date}
#' @param folder_path \emph{Optional: } User-defined path where the folder should be created. 
#' If not provided, results will be saved in `MetAlyzer_results` within the working directory. \strong{Default = NULL}
#' @param file_name Name of the output file (without extension). \strong{Default = "network"}
#' @param format File format for saving the plot (e.g., "png", "pdf", "svg"). \strong{Default = "pdf"}
#' @param width Width of the saved plot in specified units. \strong{Default = 29.7}
#' @param height Height of the saved plot in specified units. \strong{Default = 21.0}
#' @param units Units for width and height (e.g., "in", "cm", "mm"). \strong{Default = "cm"}
#' @param overwrite Logical: If `TRUE`, overwrite existing files without asking. If `FALSE`, prompt user before overwriting. \strong{Default = FALSE}
#'
#' @return The function does not return anything but saves the plot to the specified directory.
#'
#' @keywords save, plot, ggplot
#' @import ggplot2

save_plot <- function(plot,
                      folder_name = format(Sys.Date(), "%Y-%m-%d"),
                      folder_path = NULL,
                      file_name = "network",
                      format = "pdf",
                      units = "cm",
                      height = 21.0,
                      width = 29.7,
                      overwrite = FALSE) {
  
  # Don't save plot 
  if(is.null(format)) {
    return(invisible(NULL))
  }
  ##############
  ### CHECKS ###
  ##############
  # Check for invalid folder names (NULL, TRUE, FALSE)
  if (is.null(folder_name) || folder_name == "" || is.logical(folder_name)) {
    message("Invalid folder_name provided. Using today's date as folder name.")
    folder_name <- format(Sys.Date(), "%Y-%m-%d")
  }

  # Set default path if none is provided
  if (is.null(folder_path)) {
    folder_path <- file.path(getwd(), "MetAlyzer_results")
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
  } else if (!dir.exists(folder_path)) {
    message("Provided `folder_path` does not exist. Using default: ", folder_path)
    folder_path <- getwd()
  }

  # Check for invalid file names (NULL, TRUE, FALSE)
  if (is.null(file_name) || file_name == "" || is.logical(file_name)) {
    message("Invalid folder_name provided. Using default folder name: 'network'")
    file_name <- "network"
  }

  # Check for valid format
  valid_formats <- c("pdf", "png", "svg")
  if (!(format %in% valid_formats)) {
    message("Invalid format provided. Please choose from: 'pdf', 'png', 'svg'. Fallback to default 'pdf'.")
    format <- "pdf"
  }

  # Check for valid units
  valid_units <- c("cm", "in", "mm", "px")
  if (!(units %in% valid_units)) {
    message("Invalid units provided. Please choose from: 'cm', 'in', 'mm', 'px'. Fallback to default: 'cm'.")
    units <- "cm"
  }

  # Check for valid height and width (should be numeric and greater than 0)
  if (!is.numeric(height) || height <= 0) {
    message("Invalid height provided. Height must be a positive numeric value. Fallback to default: 21.0.")
    height <- 21.0
  }
  if (!is.numeric(width) || width <= 0) {
    message("Invalid width provided. Width must be a positive numeric value. Fallback to default: 29.7.")
    width <- 29.7
  }

  # Check for valid overwrite value
  if(!is.logical(overwrite)) {
    message("Invalid overwrite value provided. Overwrite must be a boolean value. Fallback to default: FALSE")
    overwrite <- FALSE
  }
  ##############
  # Remove special characters from folder name
  cleaned_folder_name <- gsub("[^a-zA-Z0-9 ]", "", folder_name)
  if (folder_name != cleaned_folder_name) {
    message("Special characters were removed from `folder_name`.")
  }

  # Create subdirectory for results
  results_folder <- file.path(folder_path, cleaned_folder_name)
  if (!dir.exists(results_folder)) {
    dir.create(results_folder)
  }

  # Construct file path
  file_path <- file.path(results_folder, paste0(file_name, ".", format))

  # Check if the file already exists
  if (file.exists(file_path)) {
    if (!overwrite) {
      response <- readline(prompt = paste("File", file_path, "already exists. Overwrite? (y/n): "))
      if (tolower(response) != "y") {
        message("File was not overwritten. Saving process canceled.")
        return(invisible(NULL))  # Exit function without saving
      }
    }
    message("Overwriting existing file: ", file_path)
  }

  # Save the plot
  ggplot2::ggsave(filename = file_path, plot = plot, width = width, height = height, units = units)
  
  message("Plot saved at: ", file_path)
}

#' Calculate Node-Level Aggregate Statistics (Conditional on Significance)
#'
#' Calculates mean log2FC, p-value, and q-value for each node (Label),
#' prioritizing significant metabolites (qval <= q_value). If none are significant,
#' uses all measured metabolites for the node. Adds results to both dataframes.
#'
#' @param nodes_sep_df Dataframe with metabolites separated (e.g., 'nodes_final').
#'   Must contain Label, log2FC, pval, qval.
#' @param nodes_orig_df Original dataframe with potentially semi-colon separated metabolites.
#'   Must contain Label.
#' @param q_value Significance threshold for q-values (e.g., 0.05).
#' @param values_col_name plotted column
#' @param stat_col_name p value column name
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
                                                  values_col_name, 
                                                  stat_col_name) {

  # --- Input Validation ---
  if (!"Label" %in% names(nodes_sep_df) || !"Label" %in% names(nodes_orig_df)) {
    stop("Both dataframes must contain a 'Label' column for grouping.")
  }
  required_cols <- c(values_col_name, stat_col_name)
  if (!all(required_cols %in% names(nodes_sep_df))) {
    stop("nodes_sep_df must contain columns: ", paste(required_cols, collapse=", "))
  }
  if (missing(q_value) || !is.numeric(q_value) || length(q_value) != 1) {
    stop("Please provide a single numeric value for q_value threshold.")
  }

  # --- 1. Calculate Node-Level Aggregates using Conditional Logic ---

  # Group by Label and apply the conditional mean logic for each metric
  node_summary <- nodes_sep_df %>%
    group_by(Label) %>%
    summarise(
      # Use the helper function for conditional mean calculation
      node_values = calculate_conditional_mean(.data[[values_col_name]], .data[[stat_col_name]], q_value),
      # Apply the same logic to qval itself: average significant qvals if present, else average measured qvals.
      node_stat   = calculate_conditional_mean(.data[[stat_col_name]],   .data[[stat_col_name]], q_value),

      .groups = "drop" # Drop grouping after summarise
    )

  # --- 2. Update nodes_separated Dataframe ---

  nodes_sep_updated <- nodes_sep_df %>%
    # Join the calculated node-level aggregates back
    left_join(node_summary, by = "Label")

  # --- 3. Update Original Nodes Dataframe ---

  nodes_orig_updated <- nodes_orig_df %>%
      dplyr::select(-any_of(c("node_values", "node_stat"))) %>%
      dplyr::left_join(node_summary, by = "Label") %>%
      dplyr::rename(!!values_col_name := node_values, !!stat_col_name := node_stat)

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
