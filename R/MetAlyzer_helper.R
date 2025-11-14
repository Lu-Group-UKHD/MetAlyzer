#' @title Save plots
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

#' @title Half-minimum imputation
#' @description Impute NA concentrations using half-minimum (HM) and update the
#' `Concentration` column in the aggregated data in the input `SummarizedExperiment`
#' (SE) object. If all values of a feature are NA, they stay NA.
#'
#' @param metalyzer_se An SE object output from \code{\link[MetAlyzer]{read_webidq()}}.
#' @returns An SE object with the imputed `Concentration` column in the aggregated
#' data accessible via `metalyzer_se@metadata$aggregated_data`.
#' 
#' @importFrom dplyr select filter mutate case_when summarise ungroup
#' 
#' @keywords internal
data_imputation <- function(metalyzer_se) {
  aggregated_data <- metalyzer_se@metadata$aggregated_data
  # Prepare HM values for corresponding features
  #### Note that aggregated data is already grouped by Metabolite by MetAlyzer
  HM_feats <- dplyr::filter(aggregated_data, Concentration >= 0, !is.na(Concentration)) %>%
    dplyr::summarise(HM = min(Concentration) * 0.5)
  # Do HM imputation, which replaces missing values with half of minimum of observed
  # values in corresponding variables
  aggregated_data <- dplyr::left_join(aggregated_data, HM_feats, by = 'Metabolite') %>%
    #### Let NaN stay NaN for now
    dplyr::mutate(Concentration = dplyr::case_when(Concentration %in% NA ~ HM,
                                                   !Concentration %in% NA ~ Concentration)) %>%
    dplyr::select(-HM) %>%
    dplyr::ungroup()
  metalyzer_se@metadata$aggregated_data <- aggregated_data
  return(metalyzer_se)
}

#' @title glog2 transformation
#' @description Perform generalized log2 transformation on data
#'
#' @param x A matrix or vector containing numerical values.
#' @returns A glog2-transformed matrix or vector.
#' 
#' @keywords internal
glog2 <- function(x) {
  (asinh(x)-log(2))/log(2)
}

#' @title Normalization
#' @description Normalize concentration values among samples using glog2 transformation,
#' median normalization, or total ion count (TIC) normalization and update the `Concentration`
#' column in the aggregated data in the input `SummarizedExperiment` (SE) object.
#'
#' @param metalyzer_se An SE object output from \code{\link[MetAlyzer]{read_webidq()}}.
#' @param norm_method A character specifying the normalization method to use, which
#' should be one of 'log2' (default), 'median', or 'TIC'.
#' @returns An SE object with the normalized `Concentration` column in the aggregated
#' data accessible via `metalyzer_se@metadata$aggregated_data`.
#' 
#' @importFrom dplyr select mutate left_join ungroup
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble as_tibble column_to_rownames
#' 
#' @keywords internal
data_normalization <- function(metalyzer_se, norm_method = 'log2') {
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
      MetAlyzer:::glog2()
  } else if (norm_method %in% 'median') {
    # Do median normalization that conducts median scaling on log2 transformed data
    # Use generalized log2 transformation to avoid -Inf
    log2_data <- MetAlyzer:::glog2(data_mat)
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
  } else if (norm_method %in% 'log2') {
    # Do log2 transformation
    norm_data <- MetAlyzer:::glog2(data_mat)
  }
  #### Exclude vsn for now to avoid any issue as it is not being used
  # else if (norm_method %in% 'VSN') {
  #   # Do vsn normalization
  #   fit <- vsn::vsnMatrix(data_mat)
  #   norm_data <- vsn::predict(fit, data_mat)
  # }
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


#' @title Get sample labels from metadata
#' @description Extracts sample labels for plotting.
#' Looks for columns "Sample ID", "Sample Identification", or "Identification"
#' (case-insensitive, ignoring punctuation). Falls back to ID column if not found.
#'
#' @param smpMetadatTbl A data frame with an 'ID' column and optional label column.
#' @return Character vector of sample labels, same length as nrow(smpMetadatTbl).
#'
#' @keywords internal
get_sample_labels <- function(smpMetadatTbl) {
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  norm_cols <- sapply(colnames(smpMetadatTbl), norm, USE.NAMES = FALSE)
  match_idx <- which(grepl("sample.*id|identif", norm_cols))
  
  if (length(match_idx) > 0) {
    labels <- as.character(smpMetadatTbl[[match_idx[1]]])
    labels[is.na(labels) | labels == ""] <- as.character(smpMetadatTbl$ID[is.na(labels) | labels == ""])
    return(labels)
  }
  as.character(smpMetadatTbl$ID)
}

#' @title Differential analysis
#' @description Perform differential analysis and add the results table to the input
#' `SummarizedExperiment` (SE) object. The analysis is conducted using the \pkg{limma}
#' package, yielding log2 fold changes, p-values, and adjusted p-values. Note that
#' the input data must already be log2 transformed, and should be subset if the
#' variable of interest contains more than two groups.
#' 
#' @param metalyzer_se An SE object output from \code{\link[MetAlyzer]{read_webidq()}}.
#' @param group A character specifying the sample metadata column containing two
#' groups that will be compared.
#' @param group_level A length-2 vector of characters specifying the group members
#' in `group`, which decides the direction of comparisons. For example, c('A', 'B')
#' compares Group A to B, and vice versa. Default is NULL (alphabetically).
#' @returns An SE object with an added table of differential analysis results, accessible
#' via `metalyzer_se@metadata$log2FC`.
#' 
#' @keywords internal
calc_log2FC <- function(metalyzer_se, group, group_level = NULL) {
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
  smp_metadata <- dplyr::select(smp_metadata, ID, all_of(group))
  # Combine abundance data and sample metadata to ensure matched information
  combined_data <- dplyr::left_join(feat_data, smp_metadata, by = 'ID')
  # Retrieve data matrix and sample metadata from combined data to conduct limma
  data_mat <- combined_data[, seq_len(ncol(feat_data))] %>%
    tibble::column_to_rownames('ID') %>%
    t()
  group_vec <- combined_data[, ncol(feat_data)+1, drop = T]
  # Sanity check if specified 'group' splits data into two groups
  if (length(unique(group_vec)) != 2) {
    stop("The specified 'group' does not split data into two groups.")
  }
  # Level groups to which samples belong to determine direction of comparison
  if (!is.null(group_level)) {
    group_vec <- relevel(factor(group_vec), ref = group_level[2])
  }
  
  # Compute log2(FC), p-values, and adjusted p-values using limma
  design <- model.matrix(~ group_vec)
  fit <- limma::lmFit(data_mat, design = design)
  fit <- limma::eBayes(fit)
  log2FCRes <- limma::topTable(fit, coef = 2, number = Inf) %>%
    tibble::rownames_to_column('Metabolite') %>%
    dplyr::select(Metabolite, logFC, t, P.Value, adj.P.Val) %>%
    dplyr::rename(log2FC = logFC, tval = t, pval = P.Value, qval = adj.P.Val)
  # Combined all information into a table
  group_info <- combined_data[, c(1, ncol(feat_data)+1)]
  log2FCTab <- dplyr::left_join(aggregated_data, group_info, by = 'ID') %>%
    dplyr::left_join(log2FCRes, by = 'Metabolite') %>%
    dplyr::select(Metabolite, Class, log2FC, tval, pval, qval) %>%
    dplyr::distinct(Metabolite, .keep_all = TRUE)
  metalyzer_se@metadata$log2FC <- log2FCTab
  return(metalyzer_se)
}

#' @title Vulcano plot - ggplotly
#' @description Create an interactive vulcano plot using `ggplotly()`.
#' 
#' @param Log2FCTab A data frame containing the differential analysis results table
#' accessible via `metalyzer_se@metadata$log2FC` where `metalyzer_se` is an SE object
#' output from \code{\link[MetAlyzer]{read_webidq()}} and has gone through `calc_log2FC()`.
#' @param x_cutoff,y_cutoff Numerical values specifying the cutoffs for log2 fold
#' changes and q-values. Default is 1.5 and 0.05, respectively.
#' @returns A `plotly` object.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom plotly ggplotly
#' @importFrom dplyr pull
#' 
#' @keywords internal
plotly_vulcano <- function(Log2FCTab, x_cutoff = 1.5, y_cutoff = 0.05) {
  # Remove metabolites with missing stats
  Log2FCTab <- Log2FCTab[!(is.na(Log2FCTab$log2FC) | is.na(Log2FCTab$qval)),]
  
  if(!"highlight" %in% colnames(Log2FCTab)) {
    # Define metabolic class colors using pre-specified color set
    polarity_file <- system.file("extdata", "polarity.csv", package = "MetAlyzer")
    metab_classes <- utils::read.csv(polarity_file) %>%
      dplyr::pull(Class) %>%
      factor()
    class_colors <- MetAlyzer::metalyzer_colors()
    names(class_colors) <- levels(metab_classes)
    # Add color for insignificant metabolites
    class_colors <- c(class_colors, c(`Not Significant` = 'grey50'))
    
    # Define color class for each metabolite
    Log2FCTab$ClassColor <- as.character(Log2FCTab$Class) #Class was factor, so 'Not Significant' could not be assigned.
    Log2FCTab$ClassColor[Log2FCTab$qval > y_cutoff] <- "Not Significant"
    Log2FCTab$ClassColor[abs(Log2FCTab$log2FC) < x_cutoff] <- "Not Significant"
    # Level ClassColor so that 'Not Significant' is at end of legend
    Log2FCTab$ClassColor <- factor(Log2FCTab$ClassColor, levels = c(levels(Log2FCTab$Class), 'Not Significant'))
    
    # Create static vulcano ggplot object
    vulcano <- ggplot(Log2FCTab, aes(x=log2FC, y=-log10(qval), color=ClassColor)) +
      geom_point(size = 2, aes(text = paste0(Metabolite,
                                             "\nClass: ", Class,
                                             "\nLog2(FC): ", round(log2FC, digits = 2),
                                             "\np-value: ", round(pval, digits = 4),
                                             "\nAdj. p-value: ", round(qval, digits = 4)))) +
      geom_vline(xintercept = c(-x_cutoff, x_cutoff), col = "black", linetype = "dashed") +
      geom_hline(yintercept = -log10(y_cutoff), col = "black", linetype = "dashed") +
      scale_color_manual('Class',
                         values = class_colors[names(class_colors) %in% Log2FCTab$ClassColor]) + #drop = F may cause confusion?
      labs(x = 'Log2(FC)', y = "-Log10(q-value)") +
      theme_bw()
    
    # Create interactive vulcano plotly object
    final_vulcano <- ggplotly(vulcano, tooltip = "text")
  } else {
    # Level highlights to make highlighted metabolites on top of plot and rename
    # factor levels for Legend understanding
    Log2FCTab$highlight <- factor(Log2FCTab$highlight, levels = c(FALSE, TRUE),
                                  labels = c("Other metabolites", "Highlighted metabolite(s)"))
    
    # Create static vulcano ggplot object with highlighted metabolites
    vulcano_highlighted <- ggplot(Log2FCTab, aes(x=log2FC, y=-log10(qval), color=highlight)) +
      geom_point(size = 2, aes(text = paste0(Metabolite,
                                             "\nClass: ", Class,
                                             "\nLog2(FC): ", round(log2FC, digits = 2),
                                             "\np-value: ", round(pval, digits = 4),
                                             "\nAdj. p-value: ", round(qval, digits = 4)))) +
      geom_vline(xintercept = c(-x_cutoff, x_cutoff), col = "black", linetype = "dashed") +
      geom_hline(yintercept = -log10(y_cutoff), col = "black", linetype = "dashed") +
      scale_color_manual('',
                         breaks = c("Other metabolites", "Highlighted metabolite(s)"),
                         values = c("grey50", "#56070C")) +
      labs(x = 'Log2(FC)', y = "-Log10(q-value)") +
      theme_bw()
    
    # Create interactive vulcano plotly object with highlighted metabolites
    final_vulcano <- ggplotly(vulcano_highlighted, tooltip = "text")
  }
  
  return(final_vulcano)
}

#' @title Vulcano plot - ggplot
#' @description Create a static vulcano plot using `ggplot()`.
#'
#' @param Log2FCTab A data frame containing the differential analysis results table
#' accessible via `metalyzer_se@metadata$log2FC` where `metalyzer_se` is an SE object
#' output from \code{\link[MetAlyzer]{read_webidq()}} and has gone through `calc_log2FC()`.
#' @param x_cutoff,y_cutoff Numerical values specifying the cutoffs for log2 fold
#' changes and q-values. Default is 1.5 and 0.05, respectively.
#' @param show_labels_for A vector of characters specifying the names of metabolites
#' or classes to label. Note that the metabolites of a specified class will be labeled.
#' Default is NULL.
#' @returns A `ggplot` object.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr pull
#' 
#' @keywords internal
plot_vulcano <- function(Log2FCTab, x_cutoff = 1.5, y_cutoff = 0.05, show_labels_for = NULL) {
  # Define metabolic class colors using pre-specified color set
  polarity_file <- system.file("extdata", "polarity.csv", package = "MetAlyzer")
  metab_classes <- utils::read.csv(polarity_file) %>%
    dplyr::pull(Class) %>%
    factor()
  class_colors <- MetAlyzer::metalyzer_colors()
  names(class_colors) <- levels(metab_classes)
  # Add color for insignificant metabolites
  class_colors <- c(class_colors, c(`Not Significant` = 'grey50'))
  
  #### Users should be interested in whatever they specify for 'show_labels_for',
  #### even if some of them are not significant
  ## Data: only color classes that are significantly differentially expressed
  # Log2FCTab$Class[Log2FCTab$qval > y_cutoff] <- NA
  # Log2FCTab$Class[abs(Log2FCTab$log2FC) < x_cutoff] <- NA
  
  # Define color class for each metabolite
  Log2FCTab$ClassColor <- as.character(Log2FCTab$Class) #Class was factor, so 'Not Significant' could not be assigned.
  Log2FCTab$ClassColor[Log2FCTab$qval > y_cutoff] <- "Not Significant"
  Log2FCTab$ClassColor[abs(Log2FCTab$log2FC) < x_cutoff] <- "Not Significant"
  # Level ClassColor so that 'Not Significant' is at end of legend
  Log2FCTab$ClassColor <- factor(Log2FCTab$ClassColor, levels = c(levels(Log2FCTab$Class), 'Not Significant'))
  
  # Prepare labels
  Log2FCTab$Label <- NA
  if (!is.null(show_labels_for)) {
    metabLabels <- Log2FCTab$Metabolite %in% show_labels_for
    metabClassLabels <- as.character(Log2FCTab$Class) %in% show_labels_for
    if (sum(metabLabels) > 0 || sum(metabClassLabels) > 0) {
      Log2FCTab$Label[metabLabels] <- Log2FCTab$Metabolite[metabLabels]
      Log2FCTab$Label[metabClassLabels] <- Log2FCTab$Metabolite[metabClassLabels]
    }
    
    # Report missing specified labels
    missLabels <- show_labels_for[!show_labels_for %in% c(Log2FCTab$Metabolite[metabLabels],
                                                          as.character(Log2FCTab$Class)[metabClassLabels])]
    if (length(missLabels) > 0) {
      print(paste("Warning: The following Metabolites/Classes were not found in the data:",
                  paste(missLabels, collapse = ", ")))
    }
  }
  
  # Create static vulcano ggplot object
  # Arrange order of metabolites in row based on significance, so that significant
  # metabolites are shown on top of insignificant ones
  Log2FCTab <- Log2FCTab[c(which(Log2FCTab$ClassColor %in% 'Not Significant'),
                           which(!Log2FCTab$ClassColor %in% 'Not Significant')),]
  vulcano <- ggplot(Log2FCTab, aes(x=log2FC, y=-log10(qval), color=ClassColor, label=Label)) +
    geom_point(size = 4) +
    geom_vline(xintercept = c(-x_cutoff, x_cutoff), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(y_cutoff), col = "black", linetype = "dashed") +
    geom_label_repel(size = 2, color = 'black', box.padding = 0.6, point.padding = 0,
                     min.segment.length = 0, max.overlaps = Inf, force = 10) +
    scale_color_manual('Class',
                       values = class_colors[names(class_colors) %in% Log2FCTab$ClassColor]) + #drop = F may cause confusion?
    labs(x = 'Log2(FC)', y = "-Log10(q-value)") +
    theme_bw()
  
  return(vulcano)
}

#' @title Scatter plot - ggplotly
#' @description Return a list containing an interactive scatter plot made using
#' `ggplotly()` and its static legend using `ggplot()`. The x-axis represents metabolic
#' classes, the y-axis shows log2 fold changes, and each point corresponds to a metabolite.
#' 
#' @param Log2FCTab A data frame containing the differential analysis results table
#' accessible via `metalyzer_se@metadata$log2FC` where `metalyzer_se` is an SE object
#' output from \code{\link[MetAlyzer]{read_webidq()}} and has gone through `calc_log2FC()`.
#' @returns A list containing a `plotly` (plot) and `ggplot` (legend) object.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom plotly ggplotly
#' @importFrom dplyr select mutate filter bind_rows group_by summarise n pull
#' @importFrom cowplot get_legend
#' @importFrom gridExtra grid.arrange
#' 
#' @keywords internal
plotly_scatter <- function(Log2FCTab) {
  # Remove metabolites with missing stats
  Log2FCTab <- Log2FCTab[!(is.na(Log2FCTab$log2FC) | is.na(Log2FCTab$qval)),]
  
  # Load polarity data containing metabolic classes and analytical techniques
  polarity_file <- system.file("extdata", "polarity.csv", package = "MetAlyzer")
  polarity_df <- utils::read.csv(polarity_file) %>%
    dplyr::select(Class, Polarity) %>%
    dplyr::mutate(Class = factor(Class))
  
  # Define metabolic class colors using pre-specified color set
  class_colors <- MetAlyzer::metalyzer_colors()
  names(class_colors) <- levels(polarity_df$Class)
  # Prepare color lists for LC- and FIA-derived classes
  lc_classes <- polarity_df[polarity_df$Polarity %in% 'LC' & polarity_df$Class %in% unique(Log2FCTab$Class), 'Class']
  lc_colors <- class_colors[names(class_colors) %in% as.character(lc_classes)]
  fia_classes <- polarity_df[polarity_df$Polarity %in% 'FIA' & polarity_df$Class %in% unique(Log2FCTab$Class), 'Class']
  fia_colors <- class_colors[names(class_colors) %in% as.character(fia_classes)]
  
  # Assign metabolites with colors indicating significance
  signif_colors = c("#5F5F5F" = 1,
                    "#FEBF6E" = 0.1,
                    "#EE5C42" = 0.05,
                    "#8B1A1A" = 0.01)
  Log2FCTab$SignifColor <- sapply(Log2FCTab$qval, function(q) {
    for (t in signif_colors) {
      if (q <= t) {
        color <- names(signif_colors)[signif_colors %in% t]
      }
    }
    return(color)
  })
  # Prepare labels for significance colors in legend
  signif_labels <- c()
  # signif_colors <- sort(signif_colors, decreasing = TRUE)
  for (i in seq_along(signif_colors)) {
    t <- signif_colors[i]
    names(t) <- NULL
    if (i < length(signif_colors)) {
      t2 <- signif_colors[i+1]
      names(t2) <- NULL
      # label <- bquote(.(t) ~ "\u2265 q >" ~ .(t2))
      label <- paste0(t, ' \u2265 q > ', t2)
    } else {
      # label <- bquote(.(t) ~ "\u2265 q")
      label <- paste0(t, ' \u2265 q')
    }
    signif_labels <- c(signif_labels, label)
  }
  
  # Create pseudo x-axis values for metabolites to determine position of each metabolite
  # and width of each metabolic class rectangle
  ordered_classes <- c(names(lc_colors), names(fia_colors))
  viz_df <- lapply(ordered_classes, function(cla) {
    Log2FCTab %>%
      dplyr::filter(Class == cla) %>%
      dplyr::bind_rows(data.frame(Class = rep(NA, 5))) #gaps between rects
  }) %>%
    dplyr::bind_rows()
  viz_df <- dplyr::bind_rows(data.frame(Class = rep(NA, 5)), viz_df)
  viz_df$x <- seq_len(nrow(viz_df))
  viz_df <- dplyr::filter(viz_df, !is.na(Class))
  
  # Create data for background class rectangles
  rects_df <- dplyr::group_by(viz_df, Class) %>%
    dplyr::summarise(Start = min(x) - 1,
                     End = max(x) + 1,
                     n = dplyr::n())
  # Relevel class based on analytical techniques (LC and FIA)
  rects_df$Class <- factor(rects_df$Class, levels = ordered_classes)
  rects_df$Technique <- sapply(rects_df$Class, function(cla) {
    if (cla %in% names(lc_colors)) {
      'LC'
    } else if (cla %in% names(fia_colors)) {
      'FIA'
    } else {
      NA
    }
  })
  # Determine position of border line between LC and FIA
  lc_fia_border <- dplyr::filter(viz_df, Class %in% names(lc_colors)) %>%
    dplyr::pull(x) %>%
    max()
  lc_fia_border <- lc_fia_border + 3
  # Determine y-axis limits of rectangle
  ylims <- c(min(Log2FCTab$log2FC) - 0.5, max(Log2FCTab$log2FC) + 0.5)
  
  # Create static scatter ggplot object
  scatter <- ggplot(viz_df, aes(x=x, y=log2FC, color=SignifColor)) +
    geom_rect(data = rects_df,
              inherit.aes = FALSE,
              aes(xmin = Start, xmax = End, ymin = ylims[1], ymax = ylims[2], fill = Class,
                  text = paste0(Class, "\nTechnique: ", Technique, "\nn(Metabolite): ", n)),
              show.legend = TRUE,
              alpha = 0.4) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_vline(xintercept = 0, color = 'black') +
    geom_vline(xintercept = lc_fia_border, linewidth = 1, color = 'black', linetype = "dotted") +
    geom_point(size = 2, aes(text = paste0(Metabolite, 
                                           "\nClass: ", Class,
                                           "\nLog2(FC): ", round(log2FC, digits = 2),
                                           "\np-value: ", round(pval, digits = 4),
                                           "\nAdj. p-value: ", round(qval, digits = 4)))) +
    scale_color_manual(paste0('Adjusted p-value (BH)'),
                       labels = signif_labels,
                       values = names(signif_colors),
                       guide = guide_legend(order = 1)) +
    scale_fill_manual('Class',
                      values = c(lc_colors, fia_colors),
                      # drop = FALSE,
                      guide = guide_legend(override.aes = list(alpha = 0.5),
                                           order = 2, ncol = 2)) +
    labs(x = 'Metabolite', y = 'Log2(FC)') +
    theme(axis.title = element_text(size = 11),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 11*0.8),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line('#ECECEC'),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_line('#ECECEC'),
          panel.background = element_blank())
  
  # Create interactive scatter plotly object
  final_scatter <- ggplotly(scatter, tooltip = "text", showlegend = FALSE)
  
  # TODO: Find god to figure out why this is not working
  ## Highlight Metabolites by changing symbol
  #if ("highlight" %in% colnames(viz_df)) {
  #  for (i in 1:nrow(viz_df)) {
  #    if (isTRUE(viz_df$highlight[i])) {
  #      print(viz_df[i,])
  #      x_val <- viz_df$x[i]
  #      y_val <- viz_df$log2FC[i]
  #      color_hex <- viz_df$SignifColor[i]
  #            
  #      final_scatter <- final_scatter %>% add_trace(
  #        x = x_val,
  #        y = y_val,
  #        type = "scatter",
  #        mode = "markers",
  #        marker = list(size = 5)
  #      )
  #    }
  #  }
  #}
  
  # Grab legend from ggplot object
  scatter_legend <- cowplot::get_legend(scatter) %>%
    gridExtra::grid.arrange(ncol = 1)
  
  return(list(Plot = final_scatter, Legend = scatter_legend))
}

#' @title Network diagram - ggplotly
#' @description This function returns a list with interactive
#' networkplot based on log2 fold change data.
#' 
#' @param Log2FCTab A dataframe with log2FC, qval, additional columns
#' @param q_value The q-value threshold for significance
#' @param values_col_name Column name of a column that holds numeric values, to be plotted \strong{Default = "log2FC"}
#' @param stat_col_name Columnname that holds numeric stat values that are used for significance \strong{Default = "qval"}
#' @param metabolite_col_name Columnname that holds the Metabolites
#' @param exclude_pathways Pathway names that are exluded from plotting
#' @param metabolite_node_size The text size of metabolite nodes
#' @param connection_width The line width of connections between metabolites
#' @param pathway_text_size The text size of pathway annotations
#' @param pathway_width The line width of pathway-specific connection coloring
#' @param plot_height The height of the Plot in pixel [px]
#' @param color_scale A string specifying the color scale to use. Options include `"viridis"`, `"plasma"`, `"magma"`, `"inferno"`, `"cividis"`, `"rocket"`, `"mako"`, and `"turbo"`, which use the `viridis` color scales.
#' 
#' @keywords internal
plotly_network <- function(Log2FCTab,
                           q_value=0.05,
                           metabolite_col_name = "Metabolite",
                           values_col_name = "log2FC",
                           stat_col_name = "qval",
                           exclude_pathways = NULL,
                           metabolite_node_size = 11,
                           connection_width = 1.25,
                           pathway_text_size = 20,
                           pathway_width = 10,
                           plot_height = 800,
                           color_scale = "viridis") {
  network_file <- MetAlyzer::pathway()
  
  ### Read in Excel file
  pathways <- MetAlyzer:::read_pathways(network_file)
  nodes <- MetAlyzer:::read_nodes(network_file, pathways)
  edges <- MetAlyzer:::read_edges(network_file, nodes, pathways)
  
  pathways <- dplyr::filter(pathways, !Pathway %in% exclude_pathways)
  nodes <- dplyr::filter(nodes, !Pathway %in% exclude_pathways)
  edges <- dplyr::filter(edges, Node1 %in% nodes$Label & Node2 %in% nodes$Label)
  
  nodes <- dplyr::filter(nodes, !Pathway %in% exclude_pathways)
  
  nodes_separated <- tidyr::separate_rows(nodes, Metabolites, sep = "\\s*;\\s*")
  
  nodes_joined <- dplyr::left_join(nodes_separated, Log2FCTab, by = c("Metabolites" = metabolite_col_name))
  
  updated_nodes_list <- MetAlyzer:::calculate_node_aggregates_conditional(nodes_sep_df = nodes_joined, nodes_orig_df = nodes, q_value = q_value, stat_col_name = stat_col_name, c("log2FC", "pval", "qval", "tval"))
  
  ### --- Create the dataframe for excel export ---
  nodes_separated_processed <- updated_nodes_list$nodes_separated
  
  nodes_separated_shortend <- nodes_separated_processed %>%
    dplyr::filter(!is.na(values_col_name))
  
  summary_all <- nodes_separated_shortend %>%
    dplyr::group_by(Label) %>%
    dplyr::filter(.data$Class != "NA") %>%
    dplyr::summarise(
      collapsed_count = dplyr::n(),
      dplyr::across(
        .cols = all_of(c("Pathway", "x", "y", "Shape", "node_log2FC", "node_pval", "node_qval", "node_tval")),
        .fns = ~ paste(unique(.), collapse = "; ")
      ),
      dplyr::across(
        .cols = !all_of(c("Pathway", "x", "y", "Shape")),
        .fns = ~ paste(., collapse = "; ")
      ),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      Label_nFeatures = paste0(Label, " -- ", collapsed_count, " feature(s)")
    )
  
  # --- The dataframe for plotting ---
  nodes_original_processed <- updated_nodes_list$nodes
  
  ## Draw network
  # Preparing Hexcodes for Annotation Colors
  nodes_original_processed$color <- MetAlyzer:::create_viridis_style(color_scale,
                                                                     type = "hex",
                                                                     data = nodes_original_processed,
                                                                     values_col_name = values_col_name)
  
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
  network <- plotly::plot_ly(nodes_original_processed,
                             x = nodes_original_processed$x,
                             y = nodes_original_processed$y,
                             type = "scatter",
                             mode = "markers",
                             marker = list(
                               color = nodes_original_processed[[values_col_name]], colorscale = MetAlyzer:::create_viridis_style(color_scale, type = "scale"),
                               showscale = TRUE,
                               colorbar = list(
                                 title = values_col_name
                               )
                             ),
                             height = plot_height)
  
  # Add the edges
  p_network <- plotly::layout(
    network,
    # title = 'Log2(FC) with FDR Correction', ####
    shapes = edges_area_combined,
    xaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
    yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
    hovermode = FALSE)
  
  # Add annotations over the nodes
  for (i in 1:nrow(nodes_original_processed)) {
    p_network <- p_network %>% plotly::add_annotations(
      text = nodes_original_processed$Label[i],
      x = nodes_original_processed$x[i],
      y = nodes_original_processed$y[i],
      arrowhead = 0,
      font = list(size = metabolite_node_size, color = "white"),
      ax = 0,
      ay = 0,
      bgcolor = nodes_original_processed$color[i],
      opacity = 1,
      hovertext = paste0("log2 Fold Change: ", round(nodes_original_processed$log2FC[i], 5),
                         "\nPathway: ", nodes_original_processed$Pathway[i],
                         "\nadj. p-value: ", round(nodes_original_processed$qval[i], 5),
                         "\np-value: ", round(nodes_original_processed$pval[i], 5),
                         "\nt-value: ", round(nodes_original_processed$tval[i], 5))
    )
  }
  
  # Add annotations for the pathways
  for (i in 1:nrow(pathways)) {
    p_network <- p_network %>% plotly::add_annotations(
      text = pathways$Pathway[i],
      x = pathways$x[i],
      y = pathways$y[i],
      xref = "x",
      yref = "y",
      showarrow = FALSE,
      font = list(size = pathway_text_size, color = pathways$Color[i])
    )
  }
  return(list("Plot" = p_network, "Table" = summary_all))
}