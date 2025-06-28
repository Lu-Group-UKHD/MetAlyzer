#' Scatter Plot Visualization
#'
#' This method creates a scatter plot of the log2 fold change for each metabolite.
#'
#' @param log2fc_df DF with metabolites as row names and columns including log2FC, Class, qval columns.
#' @param show_labels_for Vector with Strings of Metabolite names or classes.
#' @param values_col_name Column name of a column that holds numeric values, to be plotted \strong{Default = "log2FC"}
#' @param stat_col_name Columnname that holds numeric stat values that are used for significance \strong{Default = "qval"}
#' @param show_p_value Boolean Value, to color p-values according to their significance level and add a Legend \strong{Default = TRUE}
#' @param signif_colors Vector assigning significance values different colors
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
#'
#' @return ggplot object
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import SummarizedExperiment
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' log2fc_df <- readRDS(MetAlyzer::toy_diffres())
#' scatter <- MetAlyzer::plot_scatter(log2fc_df)
plot_scatter <- function(log2fc_df,
                         show_labels_for = NULL,
                         values_col_name = "log2FC",
                         stat_col_name = "qval",
                         show_p_value = TRUE,
                         signif_colors = c("#5F5F5F" = 1,
                                           "#FEBF6E" = 0.1,
                                           "#EE5C42" = 0.05,
                                           "#8B1A1A" = 0.01),
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
  
  ## Background: Load polarity data
  polarity_file <- MetAlyzer:::polarity()

  polarity_df <- utils::read.csv(polarity_file) %>%
    dplyr::select(.data$Class,
           .data$Polarity) %>%
    dplyr::mutate(Class = factor(.data$Class),
           Polarity = factor(.data$Polarity, levels = c('LC', 'FIA'))) %>%
    dplyr::arrange(.data$Polarity)

  ## Background: Set class colors

  class_colors <- MetAlyzer:::metalyzer_colors()

  names(class_colors) <- levels(polarity_df$Class)

  ## Background: Define LC and FIA classes with color
  lc_polarity_df <- filter(polarity_df,
                           .data$Polarity == 'LC',
                           .data$Class %in% log2fc_df$Class)
  lc_colors <- class_colors[which(names(class_colors) %in% lc_polarity_df$Class)]
  fia_polarity_df <- filter(polarity_df,
                            .data$Polarity == 'FIA',
                            .data$Class %in% log2fc_df$Class)
  fia_colors <- class_colors[which(names(class_colors) %in% fia_polarity_df$Class)]

  ## Data: Replace NAs
  log2fc_df[[values_col_name]][is.na(log2fc_df[[values_col_name]])] <- 0
  log2fc_df[[stat_col_name]][is.na(log2fc_df[[stat_col_name]])] <- 1

  ## Data: Add color to data based on significance
  if (isTRUE(show_p_value)) {
    log2fc_df$signif_color <- sapply(log2fc_df[[stat_col_name]], function(q_val) {
      for (t in signif_colors) {
          if (q_val <= t) {
            color <- names(signif_colors)[which(signif_colors == t)]
          }
        }
        return(color)
    })
  } else {
    signif_color <- c("black")
    log2fc_df$signif_color  <- "black"
  }

  ## Data: Add pseudo x-value to data as a order of metabolites
  ordered_classes <- c(names(lc_colors), names(fia_colors))
  p_data <- lapply(ordered_classes, function(class) {
    log2fc_df %>%
      filter(.data$Class == class) %>%
      bind_rows(data.frame(Class = rep(NA, 5)))
  }) %>%
    bind_rows()
  p_data <- bind_rows(data.frame(Class = rep(NA, 5)), p_data)
  p_data$x <- seq(nrow(p_data))
  p_data <- filter(p_data, !is.na(.data$Class))

  ## Data: Determine labels
  if (!is.null(show_labels_for)) {
    found_metabolites <- show_labels_for[show_labels_for %in% p_data$Metabolite]
    found_classes <- show_labels_for[show_labels_for %in% p_data$Class]
    
    not_found <- setdiff(show_labels_for, c(found_metabolites, found_classes))
    
    if (length(not_found) > 0) {
      print(paste("Warning: The following values were not found in the dataset:", paste(not_found, collapse = ", ")))
    }
    
    # Only assign labels to the ones in show_labels_for
    p_data$labels[p_data$Metabolite %in% show_labels_for] <- p_data$Metabolite
    p_data$labels[p_data$Class %in% show_labels_for] <- p_data$Metabolite
  }

  labels <- sapply(p_data$Metabolite, function(m) {
    m <- as.character(m)
    label <- ifelse(m %in% p_data$labels, m, "")
    return(label)
  })

  ## Legend: Significance color
  if(isTRUE(show_p_value)) {
    signif_colors <- sort(signif_colors, decreasing = TRUE)
    signif_labels <- list()
    for (i in seq_along(signif_colors)) {
      t <- signif_colors[i]
      names(t) <- NULL
      if (i < length(signif_colors)) {
        t2 <- signif_colors[i+1]
        names(t2) <- NULL
        label <- bquote(.(t) ~ "\u2265 q-value >" ~ .(t2))
      } else {
        label <- bquote(.(t) ~ "\u2265 q-value")
      }
      signif_labels[[i]] <- label
    }
  } else {
    signif_labels <- list()
  }
  

  ## Legend: Manage breaks and values for background rects
  len_diff <- length(lc_colors) - length(fia_colors)
  if (len_diff != 0) {
    blank_names <- sapply(1:abs(len_diff), function(i) {
      paste(rep(' ', i), collapse = '')
    })
    extension <- rep("white", abs(len_diff))
    names(extension) <- blank_names
    if (len_diff > 0) {
      # more classes from lc than fia
      # -> extend fia colors
      fia_colors <- c(fia_colors, extension)
    } else if (len_diff < 0) {
      # more classes from fia than lc
      # -> extend lc colors
      lc_colors <- c(lc_colors, extension)
    }
  }
  breaks <- c('LC:', names(lc_colors), 'FIA:', names(fia_colors))
  values <- c('white', lc_colors, 'white', fia_colors)
  names(values) <- NULL

  ## Background: Create data for background rects
  rects_df <- p_data %>%
    group_by(.data$Class) %>%
    summarise(Start = min(.data$x)-1,
              End = max(.data$x)+1,
              Color = class_colors[unique(.data$Class)])
  rects_df$Class <- factor(rects_df$Class, levels = breaks)

  ## Background: Determine border line between last LC and first FIA class
  lc_fia_border <- p_data %>%
    filter(.data$Class %in% names(lc_colors)) %>%
    select(.data$x) %>%
    max()

  # Create y-axis limits for the rectangles
  ylims <- c(min(log2fc_df[[values_col_name]]) - 0.75, max(log2fc_df[[values_col_name]]) + 0,75)

  ## Plot graph
  scatter <- ggplot(p_data, aes(x = .data$x,
                            y = .data[[values_col_name]],
                            color = .data$signif_color,
                            label = labels)) +
    geom_rect(data = rects_df,
              inherit.aes = FALSE,
              aes(xmin = .data$Start, xmax = .data$End,
                  ymin = ylims[1], ymax = ylims[2],
                  fill = .data$Class),
              show.legend = TRUE,
              alpha = 0.4) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = 'black') +
    geom_vline(xintercept = lc_fia_border + 3, linewidth = 0.5, color = 'black', 
              linetype = "dotted") +
    geom_hline(yintercept = 0, linewidth = 0.5, color = 'black') +
    geom_point(size = 0.5) +
    scale_color_manual(paste0('Significance\n(linear model fit with FDR correction)'),
                labels = signif_labels,
                values = names(signif_colors),
                guide = guide_legend(order = 1)) +
    scale_fill_manual('Classes',
                      breaks = breaks,
                      values = values,
                      drop = FALSE,
                      guide = guide_legend(override.aes = list(alpha = 0.5),
                                          order = 2, ncol = 2)) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          legend.key = element_rect(fill = 'white'),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line('#ECECEC'),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_line('#ECECEC'),
          panel.background = element_blank()) +
    labs(x = 'Metabolites') +
    geom_label_repel(size = 2, color = 'black',
                    box.padding = 0.6,
                    point.padding = 0,
                    min.segment.length = 0,
                    max.overlaps = Inf,
                    force = 10)

  save_plot(scatter,
          folder_name = folder_name,
          folder_path = folder_path,
          file_name = file_name,
          width = width,
          height = height,
          units = units,
          format = save_as,
          overwrite = overwrite)

  return(scatter)
}