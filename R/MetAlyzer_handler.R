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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
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
#' metalyzer_se <- MetAlyzer::read_metidq(file_path = MetAlyzer::load_demodata_biocrates())
#' metalyzer_se@metadata$log2FC <- readRDS(MetAlyzer::toy_diffres())
#' MetAlyzer::log2FC(metalyzer_se)
log2FC <- function(metalyzer_se) {
  return(metalyzer_se@metadata$log2FC)
}