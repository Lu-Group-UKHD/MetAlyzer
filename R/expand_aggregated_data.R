#' @title Expand aggregated_data
#' @description Add sample metadata into the aggregated_data in the input `SummarizedExperiment`
#' (SE) object.
#'
#' @param metalyzer_se An SE object output from \code{\link[MetAlyzer]{read_webidq()}}.
#' @param metadata_col A vector of character(s) specifying the column in the sample
#' metadata to add, accessible via `SummarizedExperiment::colData(metalyzer_se)`.
#' @return An SE object with the added sample metadata in the aggregated data accessible
#' via `metalyzer_se@metadata$aggregated_data`.
#' 
#' @importFrom dplyr select left_join relocate all_of
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @import SummarizedExperiment
#' 
#' @keywords internal
expand_aggregated_data <- function(metalyzer_se, metadata_col) {
  # Extract aggregated data and sample metadata
  aggregated_data <- metalyzer_se@metadata$aggregated_data
  metadat <- tibble::as_tibble(SummarizedExperiment::colData(metalyzer_se), rownames = 'ID')
  # Use original column names whose spaces are not replaced with '.'
  colnames(metadat) <- c('ID', colnames(SummarizedExperiment::colData(metalyzer_se)))
  
  # Sanity check
  if (!is.character(metadata_col)) {
    stop("Argument for 'metadata_col' has to be of class character.")
  }
  if (any(!metadata_col %in% colnames(metadat))) {
    tmp_col <- metadata_col[!metadata_col %in% colnames(metadat)]
    stop("Column(s) ", paste0("'", tmp_col, "'", collapse = ', '), " do not exist in sample metadata.")
  }
  if (any(metadata_col %in% colnames(aggregated_data))) {
    tmp_col <- metadata_col[metadata_col %in% colnames(aggregated_data)]
    stop("Column(s) ", paste0("'", tmp_col, "'", collapse = ', '), " already exist in aggregated_data.")
  }
  
  # Add column of interest into aggregated data
  added_col <- dplyr::select(metadat, ID, dplyr::all_of(metadata_col))
  aggregated_data <- dplyr::left_join(aggregated_data, added_col, by = 'ID') %>%
    dplyr::relocate(.data[[metadata_col]], .after = ID)
  
  # Update aggregated data
  metalyzer_se@metadata$aggregated_data <- aggregated_data
  
  return(metalyzer_se)
}
