# === Display MetAlyzer class ===

#' @title Summarize concentration values
#'
#' @description This function prints quantiles and NAs of raw data.
#'
#' @param metalyzer_se SummarizedExperiment
#' @importFrom stats quantile
#' @importFrom SummarizedExperiment rowData assay
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#'
#' MetAlyzer::summarize_conc_values(metalyzer_se)
summarize_conc_values <- function(metalyzer_se) {
  conc_values <- SummarizedExperiment::assay(
    metalyzer_se, "conc_values"
  )
  nas <- sum(is.na(conc_values))
  total <- nrow(conc_values) * ncol(conc_values)
  cat("\nMeasured concentration values:\n")
  cat("------------------------------\n")
  print(stats::quantile(conc_values, na.rm = TRUE))
  cat(paste0("\nNAs: ", nas, " (", round(nas / total * 100, 2), "%)\n"))
  classes <- SummarizedExperiment::rowData(metalyzer_se)$metabolic_classes
  if ("Metabolism Indicators" %in% classes) {
    cat("Note: 'Metabolism Indicators' are frequently NA!")
  }
  cat("\n")
}

#' @title Summarize quantification status
#'
#' @description This function lists the number of each quantification status and
#' its percentage.
#'
#' @param metalyzer_se SummarizedExperiment
#' @importFrom SummarizedExperiment assay 
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#'
#' MetAlyzer::summarize_quant_data(metalyzer_se)
summarize_quant_data <- function(metalyzer_se) {
  # Print number of quantification status
  print_number <- function(quant_status, status, total) {
    number <- sum(quant_status == status, na.rm = TRUE)
    cat(paste0(
      status, ": ",
      number, " (", round(number / total * 100, 2), "%)\n"
    ))
  }
  quant_status <- SummarizedExperiment::assay(
    metalyzer_se, "quant_status"
  )

  status_vec <- names(metalyzer_se@metadata$status_list)
  every_measured_status <- status_vec[
    which(status_vec %in% unlist(quant_status))
  ]
  nas <- sum(is.na(quant_status))
  total <- nrow(quant_status) * ncol(quant_status)
  cat("\nMeasured quantification status:\n")
  cat("-------------------------------\n")
  for (status in every_measured_status) {
    print_number(quant_status, status, total)
  }
  cat(paste0("NAs: ", nas, " (", round(nas / total * 100, 2), "%)\n"))
  cat("\n")
}


# === Handle Meta Data ===

#' @title Filter meta data
#'
#' @description This function updates the "Filter" column in meta_data to
#' filter out samples.
#'
#' @param metalyzer_se SummarizedExperiment
#' @param ... Use ´col_name´ and condition to filter selected variables.
#' @param inplace If FALSE, return a copy. Otherwise, do operation inplace and
#' return None.
#' @return An updated SummarizedExperiment
#' @importFrom dplyr mutate select filter
#' @importFrom rlang .data enquos exprs
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment colData
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#'
#' metalyzer_se <- MetAlyzer::filter_meta_data(metalyzer_se, `Sample Description` %in% 1:6)
filter_meta_data <- function(metalyzer_se, ..., inplace = FALSE) {
  # Get the parent environment
  env <- parent.frame()
  # Get the name of the metalyzer_se object
  var_name <- deparse(substitute(metalyzer_se))

  # Get meta_data and aggregated_data
  meta_data <- as.data.frame(SummarizedExperiment::colData(metalyzer_se))
  colnames(meta_data) <- gsub("\\.", " ", colnames(meta_data))

  aggregated_data <- metalyzer_se@metadata$aggregated_data

  # Filter meta_data
  conditions <- rlang::enquos(...)
  true_samples <- meta_data %>%
    dplyr::mutate(Index = rownames(meta_data)) %>%
    dplyr::filter(!!!rlang::exprs(!!!conditions)) %>%
    dplyr::select(.data$Index) %>%
    unlist()

  # Filter aggregated_data
  aggregated_data <- aggregated_data %>%
    dplyr::filter(.data$ID %in% true_samples) %>%
    droplevels()

  # Update SummarizedExperiment
  subset_condition <- rownames(meta_data) %in% true_samples
  metalyzer_se <- metalyzer_se[, subset_condition]
  metalyzer_se@metadata$aggregated_data <- aggregated_data

  # Print how many samples were removed
  diff <- nrow(meta_data) - length(true_samples)
  if (diff == 1) {
    cat("Info: 1 sample was removed!\n")
  } else if (diff > 1) {
    cat("Info:", paste(diff, "samples were removed!\n"))
  } else {
    cat("Info: No samples were removed!\n")
  }

  if (inplace) {
    # Assign the modified metalyzer_se object back to the parent environment
    assign(var_name, metalyzer_se, envir = env)
  } else {
    return(metalyzer_se)
  }
}

#' @title Update meta data
#'
#' @description This function adds another column to filtered meta_data.
#'
#' @param metalyzer_se SummarizedExperiment
#' @param ... Use ´new_col_name = new_column´ to rename selected variables
#' @param inplace If FALSE, return a copy. Otherwise, do operation inplace
#' and return None.
#' @return An updated SummarizedExperiment
#' @importFrom SummarizedExperiment colData 
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#'
#' metalyzer_se <- MetAlyzer::update_meta_data(
#'   metalyzer_se,
#'   Date = Sys.Date(), Analyzed = TRUE
#' )

update_meta_data <- function(metalyzer_se, ..., inplace = FALSE) {
  # Get the parent environment
  env <- parent.frame()
  # Get the name of the metalyzer_se object
  var_name <- deparse(substitute(metalyzer_se))

  meta_data <- SummarizedExperiment::colData(metalyzer_se)
  new_cols <- list(...)

  for (col_name in names(new_cols)) {
    new_col <- new_cols[[col_name]]
    if (inherits(new_col, "factor")) {
      levels <- levels(new_col)
    } else {
      levels <- unique(new_col)
    }
    meta_data[, col_name] <- factor(new_col, levels = levels)
  }
  SummarizedExperiment::colData(metalyzer_se) <- meta_data

  if (inplace) {
    # Assign the modified metalyzer_se object back to the parent environment
    assign(var_name, metalyzer_se, envir = env)
  } else {
    return(metalyzer_se)
  }
}

#' @title Rename meta data
#'
#' @description This function renames a column of meta_data.
#'
#' @param metalyzer_se SummarizedExperiment
#' @param ... Use new_name = old_name to rename selected variables
#' @param inplace If FALSE, return a copy. Otherwise, do operation inplace
#' and return None.
#' @return An updated SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr rename
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#'
#' metalyzer_se <- MetAlyzer::rename_meta_data(
#'   metalyzer_se,
#'   Method = `Sample Description`
#' )
rename_meta_data <- function(metalyzer_se, ..., inplace = FALSE) {
  # Get the parent environment
  env <- parent.frame()
  # Get the name of the metalyzer_se object
  var_name <- deparse(substitute(metalyzer_se))

  meta_data <- SummarizedExperiment::colData(metalyzer_se)
  meta_data <- as.data.frame(meta_data)
  colnames(meta_data) <- gsub("\\.", " ", colnames(meta_data))
  meta_data <- dplyr::rename(meta_data, ...)
  meta_data <- S4Vectors::DataFrame(meta_data)
  colnames(meta_data) <- gsub("\\.", " ", colnames(meta_data))

  SummarizedExperiment::colData(metalyzer_se) <- meta_data

  if (inplace) {
    # Assign the modified metalyzer_se object back to the parent environment
    assign(var_name, metalyzer_se, envir = env)
  } else {
    return(metalyzer_se)
  }
}

# === Manage Metabolites ===

#' @title Filter metabolites
#'
#' @description This function filters out certain classes or metabolites
#' of the metabolites vector. If aggregated_data is not empty,
#' metabolites and class will also be filtered here.
#'
#' @param metalyzer_se SummarizedExperiment
#' @param drop_metabolites A character vector defining metabolite classes
#' or individual metabolites to be removed
#' @param drop_NA_concentration A boolean whether to drop metabolites which have
#' any NAs in their concentration value
#' @param drop_quant_status A character, vector of characters or list of
#' characters specifying which quantification status to remove. Metabolites with
#' at least one quantification status of this vector will be removed.
#' @param min_percent_valid A numeric lower threshold between 0 and 1 (t less than or equal to x) to
#' remove invalid metabolites that do not meet a given percentage of valid
#' measurements per group (default per Metabolite).
#' @param valid_status A character vector that defines which quantification
#' status is considered valid.
#' @param per_group A character vector of column names from meta_data that will
#' be used to split each metabolite into groups. The threshold
#' `min_percent_valid` will be applied for each group. The selected columns from
#' meta_data will be added to aggregated_data.
#' @param inplace If FALSE, return a copy. Otherwise, do operation inplace
#' and return None.
#' @return An updated SummarizedExperiment
#' @importFrom dplyr mutate group_by_at arrange_at select group_by filter 
#' @importFrom data.table :=
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment rowData assay colData
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#'
#' drop_metabolites <- c("C0", "C2", "C3", "Metabolism Indicators",
#'   inplace = TRUE
#' )
#' metalyzer_se <- MetAlyzer::filter_metabolites(metalyzer_se, drop_metabolites)
filter_metabolites <- function(metalyzer_se,
                              drop_metabolites = c("Metabolism Indicators"),
                              drop_NA_concentration = FALSE,
                              drop_quant_status = NULL,
                              min_percent_valid = NULL,
                              valid_status = c("Valid", "LOQ"),
                              per_group = NULL,
                              inplace = FALSE) {
  # Get the parent environment
  env <- parent.frame()
  # Get the name of the metalyzer_se object
  var_name <- deparse(substitute(metalyzer_se))

  # Get metabolites and aggregated_data
  metabolites_df <- SummarizedExperiment::rowData(metalyzer_se)
  metabolites <- as.vector(rownames(metabolites_df))
  classes <- as.vector(metabolites_df$metabolic_classes)
  aggregated_data <- metalyzer_se@metadata$aggregated_data

  # Vector of metabolites that will be removed
  rm_metabolites <- c()

  # Get all metabolites that given by their name or their class
  if (!is.null(drop_metabolites)) {
    drop_metabolites <- unlist(as.vector(drop_metabolites))
    rm_metabolites <- c(
      rm_metabolites,
      drop_metabolites[
        which(drop_metabolites %in% metabolites)
      ]
    )
    rm_classes <- drop_metabolites[
      which(drop_metabolites %in% classes)
    ]
    rm_metabolites <- c(
      rm_metabolites,
      metabolites[
        which(classes %in% rm_classes)
      ]
    )
  }

  # Get all metabolites that have the quantification status at least once
  if (!is.null(drop_quant_status)) {
    quant_status <- SummarizedExperiment::assay(
      metalyzer_se, "quant_status"
    )

    drop_quant_status <- unlist(as.vector(drop_quant_status))
    metabolites_per_status <- lapply(drop_quant_status, function(status) {
      if (status == "NA") {
        n <- rowSums(is.na(quant_status))
      } else {
        n <- rowSums(quant_status == status)
      }
      status_metabolites <- rownames(quant_status)[n > 0]
      return(status_metabolites)
    })
    rm_metabolites <- c(
      rm_metabolites,
      as.vector(unlist(metabolites_per_status))
    )
  }

  # Get all metabolites whose percentage of valid quantification status does
  # not meet a given lower threshold value
  if (!is.null(min_percent_valid)) {
    stopifnot(0 <= min_percent_valid && min_percent_valid <= 1)
    grouping_vars <- c("Metabolite")

    # If additional groups of meta_data are given add them to aggregated data
    if (!is.null(per_group)) {
      meta_data <- SummarizedExperiment::colData(metalyzer_se)
      per_group <- unlist(as.vector(per_group))
      grouping_vars <- c(grouping_vars, per_group)
      for (group in rev(per_group)) {
        mapping_vec <- unlist(meta_data[group])
        names(mapping_vec) <- rownames(meta_data[group])
        aggregated_data <- dplyr::mutate(
          aggregated_data,
          !!group := factor(
            sapply(.data$ID, function(id) {
              mapping_vec[id]
            }),
            levels = unique(mapping_vec)
          ),
          .after = .data$ID
        )
      }
      cat(
        "Info: A group counts as valid, if at least ",
        round(min_percent_valid * 100, 2),
        "% of samples are considered as valid.\n",
        "Group-wise calculation: (",
        paste(grouping_vars, collapse = " * "), ")\n",
        "A metabolite counts as invalid, if all groups are invalid.\n",
        sep = ""
      )
    } else {
      cat(
        "Info: A metabolite counts as valid, if at least ",
        round(min_percent_valid * 100, 2),
        "% of samples are considered as valid.\n",
        sep = ""
      )
    }

    aggregated_data <- aggregated_data %>%
      dplyr::group_by_at(grouping_vars) %>%
      dplyr::arrange_at(c(grouping_vars, "ID")) %>%
      dplyr::mutate(
        valid_status = .data$Status %in% valid_status,
        Valid_Group = min_percent_valid <= sum(valid_status) / n(),
        .after = .data$Status
      ) %>%
      dplyr::select(-valid_status)

    invalid_metabolites <- aggregated_data %>%
      dplyr::group_by(.data$Metabolite) %>%
      dplyr::filter(sum(.data$Valid_Group) == 0) %>%
      dplyr::select(.data$Metabolite) %>%
      unlist()
    rm_metabolites <- c(
      rm_metabolites,
      as.vector(invalid_metabolites)
    )
  }

  # Get all metabolites with at least one concentration being NA
  if (drop_NA_concentration) {
    conc_values <- SummarizedExperiment::assay(
      metalyzer_se, "conc_values"
    )
    n_nas <- rowSums(is.na(conc_values))
    na_metabolites <- rownames(conc_values)[n_nas > 0]

    rm_metabolites <- c(
      rm_metabolites,
      as.vector(na_metabolites)
    )
  }

  # Remove metabolites
  rm_metabolites <- unique(rm_metabolites)
  if (length(rm_metabolites) > 0) {
    aggregated_data <- dplyr::filter(
      aggregated_data,
      !(.data$Metabolite %in% rm_metabolites)
    ) %>%
      droplevels()

    # Update SummarizedExperiment
    subset_condition <- !metabolites %in% rm_metabolites
    metalyzer_se <- metalyzer_se[subset_condition, ]
    metalyzer_se@metadata$aggregated_data <- aggregated_data
  }

  diff <- length(rm_metabolites)
  if (diff == 1) {
    cat("Info: 1 metabolite was removed!\n")
  } else if (diff > 1) {
    cat("Info:", paste(diff, "metabolites were removed!\n"))
  } else {
    cat("Info: No metabolites were removed!\n")
  }

  if (inplace) {
    # Assign the modified metalyzer_se object back to the parent environment
    assign(var_name, metalyzer_se, envir = env)
  } else {
    return(metalyzer_se)
  }
}

# === Handle Aggregated Data ===
#' @title Get Aggregated Data
#'
#' @description This function returns the tibble "aggregated_data".
#'
#' @param metalyzer_se SummarizedExperiment
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#'
#' MetAlyzer::aggregated_data(metalyzer_se)
aggregated_data <- function(metalyzer_se) {
  return(metalyzer_se@metadata$aggregated_data)
}

# === Export data ===

#' @title Export filtered raw data as csv
#'
#' @description This function exports the filtered raw data in the CSV format.
#'
#' @param metalyzer_se SummarizedExperiment
#' @param ... Additional columns from meta_data
#' @param file_path file path
#' @importFrom dplyr bind_cols select
#' @importFrom SummarizedExperiment colData assay 
#' @importFrom utils write.csv
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#' 
#' output_file <- file.path(tempdir(), "metabolomics_data.csv")
#' MetAlyzer::export_conc_values(metalyzer_se,
#'                               `Sample Description`,
#'                               file_path = output_file
#'                               )
#' unlink(output_file)
export_conc_values <- function(metalyzer_se,
                             ...,
                             file_path = "metabolomics_data.csv") {
  meta_data <- as.data.frame(SummarizedExperiment::colData(metalyzer_se))
  colnames(meta_data) <- gsub("\\.", " ", colnames(meta_data))
  conc_values <- SummarizedExperiment::assay(
    metalyzer_se, "conc_values"
  )
  df <- dplyr::bind_cols(
    dplyr::select(meta_data, ...),
    t(conc_values)
  )
  cat("Info: Number of samples:", nrow(meta_data), "\n")
  cat("Info: Number of Metabolites:", nrow(conc_values), "\n")
  utils::write.csv(
    x = df,
    file = file_path,
    row.names = FALSE
  )
}

# === Handle log2FC Data ===
#' @title Get log2FC Data
#'
#' @description This function returns the tibble "log2FC".
#'
#' @param metalyzer_se SummarizedExperiment
#' @export
#' 
#' @examples
#' metalyzer_se <- MetAlyzer::read_metidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#' metalyzer_se@metadata$log2FC <- readRDS(MetAlyzer::toy_diffres())
#' MetAlyzer::log2FC(metalyzer_se)
log2FC <- function(metalyzer_se) {
  return(metalyzer_se@metadata$log2FC)
}

# === Helper functions ===

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
#' @param metalyzer_se An SE object output from \code{\link[MetAlyzer]{read_metidq()}}.
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
#' median normalization, total ion count (TIC) normalization, or variance stabilizing
#' normalization (VSN) and update the `Concentration` column in the aggregated data
#' in the input `SummarizedExperiment` (SE) object.
#'
#' @param metalyzer_se An SE object output from \code{\link[MetAlyzer]{read_metidq()}}.
#' @param norm_method A character specifying the normalization method to use, which
#' should be one of 'log2' (default), 'median', 'TIC', or 'VSN'.
#' @returns An SE object with the normalized `Concentration` column in the aggregated
#' data accessible via `metalyzer_se@metadata$aggregated_data`.
#' 
#' @importFrom dplyr select mutate left_join ungroup
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble as_tibble column_to_rownames
#' @import vsn
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
  } else if (norm_method %in% 'VSN') {
    # Do vsn normalization
    fit <- vsnMatrix(data_mat)
    norm_data <- predict(fit, data_mat)
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

#' @title Differential analysis
#' @description Perform differential analysis and add the results table to the input
#' `SummarizedExperiment` (SE) object. The analysis is conducted using the \pkg{limma}
#' package, yielding log2 fold changes, p-values, and adjusted p-values. Note that
#' the input data must already be log2 transformed, and should be subset if the
#' variable of interest contains more than two groups.
#' 
#' @param metalyzer_se An SE object output from \code{\link[MetAlyzer]{read_metidq()}}.
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