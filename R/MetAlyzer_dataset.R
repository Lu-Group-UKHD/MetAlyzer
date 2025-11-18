#' @title Open file and read data
#'
#' @description This function creates a SummarizedExperiment (SE) from the given
#' 'webidq' output Excel sheet: metabolites (rowData), meta data (colData),
#' concentration data (assay), quantification status(assay)
#' The column "Sample Type" and the row "Class" are used as anchor cells in the
#' Excel sheet and are therefore a requirement.
#'
#' @param file_path A character specifying the file path to the Excel file.
#' @param sheet A numeric index specifying which sheet of the Excel file to use.
#' @param status_list A list of HEX color codes for each quantification status.
#' @param silent If TRUE, mute any print command.
#' @return A Summarized Experiment object
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' Path <- MetAlyzer::load_demodata_biocrates()
#' metalyzer_se <- MetAlyzer::read_webidq(file_path = Path)
read_webidq <- function(
    file_path,
    sheet = 1,
    status_list = list(
      "Valid" = c("#B9DE83", "#00CD66"),
      "LOQ" = c("#B2D1DC", "#7FB2C5", "#87CEEB"), # <LLOQ or > ULOQ
      "LOD" = c("#A28BA3", "#6A5ACD", "#BBA7B9"), # < LOD
      "ISTD Out of Range" = c("#FFF099", "#FFFF33"),
      "Invalid" = "#FFFFCC",
      "Incomplete" = c("#CBD2D7", "#FFCCCC")
    ),
    silent = FALSE) {
  # Input checks
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("`file_path` must be a single character string.", call. = FALSE)
  }
  if (!file.exists(file_path)) {
    stop("`file_path` does not exist: ", file_path, call. = FALSE)
  }
  if (!((is.numeric(sheet) && length(sheet) == 1 && sheet > 0 && floor(sheet) == sheet) ||
        (is.character(sheet) && length(sheet) == 1))) {
    stop("`sheet` must be a single positive integer or a single sheet name (character string).", call. = FALSE)
  }
  if (!is.list(status_list)) {
    stop("`status_list` must be a list.", call. = FALSE)
  }
  if (length(status_list) > 0 && !all(sapply(status_list, function(x) is.character(x) && !any(is.na(x))))) {
    stop("All elements in `status_list` must be character vectors of hex color codes.", call. = FALSE)
  }
  if (!is.logical(silent) || length(silent) != 1) {
    stop("`silent` must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  # Rest of the function code would go here...
  if (!silent) {
    message("Input checks passed. Proceeding with reading webidq file.")
  }
  # Print MetAlyzer logo
  if (silent == FALSE) {
    metalyzer_ascii_logo()
  }

  # Open webidq Excel sheet
  starter_list <- list(
    "file_path" = as.character(file_path),
    "sheet" = as.numeric(sheet)
  )
  full_sheet <- open_file(starter_list)
  data_ranges <- get_data_range(full_sheet)

  # Extract metabolites, meta data and concentration values
  metabolites <- slice_metabolites(full_sheet, data_ranges)
  meta_data <- slice_meta_data(full_sheet, data_ranges)
  conc_values <- slice_conc_values(full_sheet, data_ranges, metabolites)

  # Read quantification status
  quant_status <- read_quant_status(
    starter_list = starter_list,
    sheet_dim = c(nrow(full_sheet), ncol(full_sheet)),
    data_ranges = data_ranges,
    metabolites = metabolites,
    status_list = status_list,
    silent = silent
  )

  # Aggregate data and add it to the metadata of SE object
  aggregated_data <- aggregate_data(
    metabolites = metabolites,
    meta_data = meta_data,
    conc_values = conc_values,
    quant_status = quant_status,
    status_vec = names(status_list)
  )

  # Fill Summarized Experiment
  # rowData: features
  # colData: samples
  # assays: conc_values, quant_status
  # metadata: file_path, sheet, aggregated_data

  rowData <- data.frame(
    "metabolic_classes" = names(metabolites),
    row.names = metabolites
  )
  colData <- meta_data
  assays <- list(
    "conc_values" = t(conc_values),
    "quant_status" = t(quant_status)
  )
  metadata <- list(
    "file_path" = file_path,
    "sheet_index" = sheet,
    "status_list" = status_list,
    "aggregated_data" = aggregated_data
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = colData,
    rowData = rowData,
    metadata = metadata,
  )

  # Print summary of conc_values and quant_status
  if (silent == FALSE) {
    summarize_conc_values(se)
    summarize_quant_data(se)
  }

  return(se)
}

#' @title MetAlyzer logo
#'
#' @description This function prints "MetAlyzer" as an ASCII logo
#'
#' @keywords internal
metalyzer_ascii_logo <- function() {
  line1 <- "\n"
  line2a <- " _____ ______   _______  _________  ________  ___           ___   "
  line2b <- " ___ ________  _______   ________"
  line3a <- "|\\   _ \\  _   \\|\\  ___ \\|\\___   ___\\\\   __  \\|\\  \\     "
  line3b <- "    |\\  \\  /  /|\\_____  \\|\\  ___ \\ |\\   __  \\"
  line4a <- "\\ \\  \\\\\\__\\ \\  \\ \\   __/\\|___ \\  \\_\\ \\  \\|\\  \\ \\"
  line4b <- "  \\        \\ \\  \\/  / /\\|___/  /\\ \\   __/|\\ \\  \\|\\  \\"
  line5a <- " \\ \\  \\\\|__| \\  \\ \\  \\_|/__  \\ \\  \\ \\ \\   __  \\ \\  "
  line5b <- "\\        \\ \\    / /     /  / /\\ \\  \\_|/_\\ \\   _  _\\"
  line6a <- "  \\ \\  \\    \\ \\  \\ \\  \\_|\\ \\  \\ \\  \\ \\ \\  \\ \\  \\"
  line6b <- " \\  \\____    \\/   / /     /  /_/__\\ \\  \\_|\\ \\ \\  \\\\  \\"
  line6c <- "| "
  line7a <- "   \\ \\__\\    \\ \\__\\ \\_______\\  \\ \\__\\ \\ \\__\\ \\__\\ "
  line7b <- "\\_______\\__/   / /     |\\________\\ \\_______\\ \\__\\\\ _\\ "
  line8a <- "    \\|__|     \\|__|\\|_______|   \\|__|  \\|__|\\|__|\\|_______|"
  line8b <- "\\____/ /       \\|_______|\\|_______|\\|__|\\|__|"
  line9 <- "                                                          \\|____|/"
  line10 <- "\n"

  logo_string <- c(
    line1,
    paste(line2a, line2b, sep = ""),
    paste(line3a, line3b, sep = ""),
    paste(line4a, line4b, sep = ""),
    paste(line5a, line5b, sep = ""),
    paste(line6a, line6b, line6c, sep = ""),
    paste(line7a, line7b, sep = ""),
    paste(line8a, line8b, sep = ""),
    line9,
    line10
  )
  cat(logo_string, sep = "\n")
}


#' @title Open Excel file
#'
#' @description This function opens the given webidq output Excel file and reads the full
#' given sheet.
#'
#' @param starter_list contains the file path and the sheet index
#' @importFrom openxlsx read.xlsx
#'
#' @keywords internal
open_file <- function(starter_list) {
  full_sheet <- openxlsx::read.xlsx(starter_list$file_path,
    sheet = starter_list$sheet,
    colNames = FALSE,
    skipEmptyRows = FALSE,
    skipEmptyCols = FALSE
  )
  full_sheet <- as.matrix(full_sheet)
  full_sheet[full_sheet == "NA"] <- NA # all entries are strings
  return(full_sheet)
}


#' @title Get data range
#'
#' @description This function extracts rows and column indices to slice the .full_sheet into
#' metabolites, conc_values and meta_data.
#'
#' @param full_sheet full_sheet
#'
#' @keywords internal
get_data_range <- function(full_sheet) {
  if (!"Class" %in% full_sheet) {
    stop('The cell "Class" is missing. Metabolites could not be detected...')
  }
  # row of header "Class"
  row_class <- which(full_sheet == "Class") %% nrow(full_sheet)
  # column of header "Class"
  col_class <- ceiling(which(full_sheet == "Class") / nrow(full_sheet))

  if (!any(grepl("^\\s*Sample\\s*Type\\s*$", full_sheet, ignore.case = TRUE))) {
    stop('The cell "Sample Type" is missing. Data could not be detected...')
  }

  # row of "Sample Type"
  row_sample_type <- arrayInd(
    which(grepl("^\\s*Sample\\s*Type\\s*$", full_sheet, ignore.case = TRUE)),
    dim(full_sheet)
  )[, 1]
  # column of "Sample Type"
  col_sample_type <- arrayInd(
    which(grepl("^\\s*Sample\\s*Type\\s*$", full_sheet, ignore.case = TRUE)),
    dim(full_sheet)
  )[, 2]

  rows_data <- which(full_sheet[, col_sample_type] == "Sample") # rows
  names(rows_data) <- NULL
  cols_data <- (col_class + 1):ncol(full_sheet) # columns

  data_ranges <- list(
    "class_row" = row_class,
    "class_col" = col_class,
    "sample_type_row" = row_sample_type,
    "sample_type_col" = col_sample_type,
    "data_rows" = rows_data,
    "data_cols" = cols_data
  )
  return(data_ranges)
}


#' @title Slice metabolites
#'
#' @description This function extracts metabolites with their corresponding
#' metabolite class from .full_sheet into metabolites. It robustly finds the
#' metabolite row by searching upwards from the "Class" row and skipping
#' any known intermediate rows (e.g., "Synonym").
#'
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#'
#' @keywords internal
slice_metabolites <- function(full_sheet, data_ranges) {
  # Add any other potential junk row headers to this vector.
  junk_patterns <- c("^\\s*Synonym(s)?\\s*$") # Matches "Synonym" or "Synonyms"
  current_row_idx <- data_ranges$class_row - 1

  while (current_row_idx > 0) {
    cell_to_check <- full_sheet[current_row_idx, data_ranges$class_col]
    if (is.na(cell_to_check)) {
      cell_to_check <- ""
    }
    is_junk <- any(sapply(junk_patterns, grepl, x = cell_to_check, ignore.case = TRUE))
    if (is_junk) {
      current_row_idx <- current_row_idx - 1
    } else {
      break
    }
  }
  metabolite_row_idx <- current_row_idx
  if (metabolite_row_idx <= 0) {
    stop(paste(
      "Could not find metabolite names row.",
      "Searched upwards from 'Class' row (", data_ranges$class_row, ")",
      "and found no valid row."
    ), call. = FALSE)
  }
  metabolites <- full_sheet[
    metabolite_row_idx,
    data_ranges$data_cols
  ]
  classes <- full_sheet[
    data_ranges$class_row,
    data_ranges$data_cols
  ]
  if (any(is.na(metabolites))) {
    warning(paste0(
      "NA values found in the identified metabolite names row (Row ",
      metabolite_row_idx,
      "). This may indicate an incorrect row was parsed."
    ), call. = FALSE)
  }

  names(metabolites) <- classes
  return(metabolites)
}


#' @title Slice concentration values
#'
#' @description This function slices measurements from .full_sheet into
#' conc_values.
#'
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#' @param metabolites metabolites
#'
#' @keywords internal
slice_conc_values <- function(full_sheet, data_ranges, metabolites) {
  conc_values <- full_sheet[
    data_ranges$data_rows,
    data_ranges$data_cols
  ]
  conc_values <- as.data.frame(conc_values)
  colnames(conc_values) <- metabolites
  conc_values[] <- lapply(conc_values[], as.numeric)
  return(conc_values)
}


#' @title Slice meta data
#'
#' @description This function slices meta data from .full_sheet into meta_data.
#'
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#'
#' @keywords internal
slice_meta_data <- function(full_sheet, data_ranges) {
  meta_header <- paste(full_sheet[
    data_ranges[["sample_type_row"]],
    1:data_ranges[["class_col"]]
  ])
  meta_data <- as.data.frame(full_sheet[
    data_ranges$data_rows,
    1:data_ranges[["class_col"]]
  ])
  colnames(meta_data) <- gsub("[^a-zA-Z ]", "", meta_header)
  meta_data <- meta_data[
    colSums(is.na(meta_data)) != nrow(meta_data)
  ] # remove full NA columns
  nas <- grep("NA", colnames(meta_data))
  colnames(meta_data)[nas] <- paste0(
    "X", seq_along(nas)
  ) # replace NAs in colnames with X and consecutive numbers
  return(meta_data)
}


#' @title Unify hex codes
#'
#' @description This function gets hex codes of different length and returns
#' them in a unified format.
#'
#' @param hex A 3, 4, 6 or 8 digit hex code
#' @importFrom stringr str_extract_all
#'
#' @keywords internal
unify_hex <- function(hex) {
  upper_hex <- toupper(hex)
  hex_digits <- stringr::str_extract_all(upper_hex, "[0-9A-F]")[[1]]
  n_digits <- length(hex_digits)
  if (n_digits == 6) {
    hex_code <- paste0("#", paste(hex_digits, collapse = ""))
  } else if (n_digits == 8) {
    hex_code <- paste0("#", paste(hex_digits[3:8], collapse = ""))
  } else if (n_digits == 3) {
    hex_code <- paste0(
      "#", hex_digits[1], hex_digits[1],
      hex_digits[2], hex_digits[2],
      hex_digits[3], hex_digits[3]
    )
  } else if (n_digits == 4) {
    hex_code <- paste0(
      "#", hex_digits[2], hex_digits[2],
      hex_digits[3], hex_digits[3],
      hex_digits[4], hex_digits[4]
    )
  } else {
    hex_code <- ""
  }
  return(hex_code)
}


#' @title Read quantification status
#'
#' @description This function gets the background color of each cell in
#' full_sheet and assigns it to the corresponding quantification status.
#'
#' @param starter_list starter_list
#' @param sheet_dim sheet_dim
#' @param data_ranges data_ranges
#' @param metabolites metabolites
#' @param status_list status_list
#' @param silent silent
#' @importFrom openxlsx loadWorkbook getSheetNames
#'
#' @keywords internal
read_quant_status <- function(
    starter_list,
    sheet_dim,
    data_ranges,
    metabolites,
    status_list,
    silent) {
  wb <- openxlsx::loadWorkbook(file = starter_list$file_path)
  sheet_name <- openxlsx::getSheetNames(
    file = starter_list$file_path
  )[starter_list$sheet]
  styles <- wb$styleObjects
  mat_bg <- matrix(NA, nrow = sheet_dim[1], ncol = sheet_dim[2])
  for (x in styles) { # iterate through different cell styles in the sheet
    if (x$sheet == sheet_name) {
      rgb <- x$style$fill$fillFg
      rgb <- ifelse(length(rgb) == 0, "", rgb)
      hex_code <- unify_hex(rgb)
      matching_status <- names(status_list)[
        sapply(status_list, function(colors) {
          hex_code %in% colors
        })
      ]
      if (length(matching_status) > 0) {
        if (hex_code != rgb && silent == FALSE) {
          cat(
            "Info: Reading color code \"", rgb, "\" as \"", hex_code, "\"",
            "\n",
            sep = ""
          )
        }
        status <- matching_status[1] # take the first matching status
        rows <- x$rows # row indices
        cols <- x$cols # column indices
        for (i in seq_along(rows)) {
          mat_bg[rows[i], cols[i]] <- status
        }
      }
    }
  }
  quant_status <- as.data.frame(mat_bg)[
    data_ranges$data_rows,
    data_ranges$data_cols
  ]
  colnames(quant_status) <- metabolites
  return(quant_status)
}

#' @title Aggregate data
#'
#' @description This function reshapes conc_values, quant_status,
#' metatabolites and sample IDs and combines them into a tibble data frame
#' for filtering with dplyr and plotting with 'ggplot2'. "aggregated_data"
#' is grouped by metabolites.
#'
#' @param metabolites metabolites MetAlyzer object
#' @param meta_data Meta_data of the MetAlyzer object
#' @param conc_values conc_values of a MetAlyzer object
#' @param quant_status quant_status of a MetAlyzer object
#' @param status_vec A vector of quantification status
#' @importFrom dplyr arrange mutate group_by
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @keywords internal
aggregate_data <- function(
    metabolites,
    meta_data,
    conc_values,
    quant_status,
    status_vec) {
  classes <- names(metabolites)

  comb_data <- dplyr::mutate(
    conc_values,
    ID = rownames(meta_data),
    .before = 1
  )
  gathered_data <- tidyr::gather(comb_data,
    key = "Metabolite",
    value = "Concentration",
    -"ID"
  )
  gathered_status <- tidyr::gather(quant_status,
    key = "Metabolite",
    value = "Status"
  )

  aggregated_data <- gathered_data %>%
    dplyr::group_by(.data$Metabolite) %>%
    dplyr::mutate(
      Class = sapply(.data$Metabolite, function(x) {
        classes[metabolites == x]
      }),
      .after = .data$Metabolite
    )

  aggregated_data$ID <- factor(aggregated_data$ID,
    levels = rownames(meta_data)
  )
  aggregated_data$Metabolite <- factor(aggregated_data$Metabolite,
    levels = unique(metabolites)
  )
  aggregated_data$Class <- factor(aggregated_data$Class,
    levels = unique(classes)
  )
  aggregated_data$Status <- factor(gathered_status$Status,
    levels = status_vec
  )
  aggregated_data <- dplyr::arrange(aggregated_data, .data$Metabolite)

  return(droplevels(aggregated_data))
}

#' @title Open file and read data
#' @param ... Declare this function as out of date
#'
#' @description This function was deprecated in version v2.0.0
#'
MetAlyzer_dataset <- function(...) {
  cat("This function was deprecated in v2.0.0, please use read_webidq()\n")
}
