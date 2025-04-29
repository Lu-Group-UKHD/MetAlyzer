#' @title Get example extraction data
#'
#' @description This function returns the extraction_data_MxP_Quant_500.xlsx file path.
#'
#' @return extraction_data_MxP_Quant_500.xlsx file path
#' @export
#'
#' @examples
#' fpath <- load_rawdata_extraction()
load_rawdata_extraction <- function() {
  system.file("extdata", "extraction_data_MxP_Quant_500.xlsx", package = "MetAlyzer")
}


#' @title Get demodata from biocrates
#'
#' @description This function returns the Metalyzer_demo dataset_biocrates MxP Quant 500 XL_2025-04.xlsx file path.
#'
#' @return Metalyzer_demo dataset_biocrates MxP Quant 500 XL_2025-04 file path
#' @export
#'
#' @examples
#' fpath <- load_demodata_biocrates()
load_demodata_biocrates <- function() {
  system.file("extdata", "Metalyzer_demo dataset_biocrates MxP Quant 500 XL_2025-04.xlsx", package = "MetAlyzer")
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
#' @description This function returns the log2fc dataframe of the Metalyzer_demo dataset_biocrates MxP Quant 500 XL_2025-04 file, created with the MetaVizPro package.
#'
#' @return toy_diffres.rds file path
#' @export
#'
#' @examples
#' fpath <- toy_diffres()
toy_diffres <- function() {
  system.file("extdata", "toy_diffres_biocrates.rds", package = "MetAlyzer")
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
