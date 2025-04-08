#' @title Get example extraction data
#'
#' @description This function returns the extraction_data_MxP_Quant_500.xlsx file path.
#'
#' @return extraction_data_MxP_Quant_500.xlsx file path
#' @export
#'
#' @examples
#' fpath <- example_extraction_data()
example_extraction_data <- function() {
  system.file("extdata", "extraction_data_MxP_Quant_500.xlsx", package = "MetAlyzer")
}


#' @title Get example mutation data
#'
#' @description This function returns the mutation_data_MxP_Quant_500_XL.xlsx file path.
#'
#' @return mutation_data_MxP_Quant_500_XL.xlsx file path
#' @export
#'
#' @examples
#' fpath <- example_mutation_data_xl()
example_mutation_data_xl <- function() {
  system.file("extdata", "mutation_data_MxP_Quant_500_XL.xlsx", package = "MetAlyzer")
}

#' @title Get example log2fc data
#'
#' @description This function returns the log2fc dataframe of the mutation_data_MxP_Quant_500_XL.xlsx file, created with the MetaVizPro package.
#'
#' @return example_log2fc_data.rds file path
#' @export
#'
#' @examples
#' fpath <- example_log2fc_data()
example_log2fc_data <- function() {
  system.file("extdata", "example_log2fc_data.rds", package = "MetAlyzer")
}


#' @title Get example meta data
#'
#' @description This function returns the data frame loaded from example_meta_data.RDS.
#'
#' @return data frame loaded from example_meta_data.RDS
#' @export
#'
#' @examples
#' fpath <- example_meta_data()
example_meta_data <- function() {
  readRDS(system.file("extdata", "example_meta_data.RDS", package = "MetAlyzer"))
}

#' @title Get MetAlyzer colors
#'
#' @description This function returns the vector loaded from metalyzer_colors.RDS.
#'
#' @return data frame loaded from metalyzer_colors.RDS
#' @export
#'
#' @examples
#' fpath <- metalyzer_colors()
metalyzer_colors <- function() {
  readRDS(system.file("extdata", "metalyzer_colors.RDS", package = "MetAlyzer"))
}

#' @title Get polarity file path
#'
#' @description This function returns the polarity.csv file path.
#'
#' @return polarity.csv file path
#' @export
#'
#' @examples
#' fpath <- polarity()
polarity <- function() {
  system.file("extdata", "polarity.csv", package = "MetAlyzer")
}


#' @title Get pathway file path
#'
#' @description This function returns the latest pathway.xlsx file path.
#'
#' @return pathway.xlsx file path
#' @export
#'
#' @examples
#' fpath <- pathway()
pathway <- function() {
  system.file("extdata", "Pathway_120325.xlsx", package = "MetAlyzer")
}
