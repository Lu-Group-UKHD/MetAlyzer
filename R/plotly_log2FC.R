#' @title Plotly Log2FC Scatter Plot
#'
#' @description This function returns a list with an interactive 
#' scatterplot based on log2 fold change data and a comprehensive Legend.
#' 
#' @param metalyzer_se A Metalyzer object
#' @param signif_colors signif_colors
#' @param class_colors A csv file containing class colors hexcodes
#' 
#' @return plotly object
#'
#' @import dplyr
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom plotly ggplotly
#' @export
#'
#' @examples
#' 
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_mutation_data_xl())
#' metalyzer_se <- filterMetabolites(
#'   metalyzer_se,
#'   drop_metabolites = "Metabolism Indicators"
#' )
#' metalyzer_se <- renameMetaData(
#'   metalyzer_se,
#'   Mutant_Control = "Sample Description"
#' )
#' metalyzer_se <- calculate_log2FC(
#'   metalyzer_se,
#'   categorical = "Mutant_Control",
#'   impute_perc_of_min = 0.2,
#'   impute_NA = TRUE
#' )
#' 
#' p_scatter <- plotly_scatter(metalyzer_se)
plotly_scatter <- function(metalyzer_se, 
    signif_colors=c("#5F5F5F"=1,
                    "#FEBF6E"=0.1,
                    "#EE5C42"=0.05,
                    "#8B1A1A"=0.01),
    class_colors = metalyzer_colors()) {
  
    ### Data Wrangling
    log2FC_df <- metalyzer_se@metadata$log2FC
    ## Background: Load polarity data
    polarity_file <- system.file("extdata", "polarity.csv", package = "MetAlyzer")
    
    polarity_df <- utils::read.csv(polarity_file) %>%
        select(.data$Class,
            .data$Polarity) %>%
        mutate(Class = factor(.data$Class),
            Polarity = factor(.data$Polarity, levels = c('LC', 'FIA'))) %>%
        arrange(.data$Polarity)
    
    ## Background: Set class colors
    class_colors <- metalyzer_colors()
    
    names(class_colors) <- levels(polarity_df$Class)
    
    ## Background: Define LC and FIA classes with color
    lc_polarity_df <- filter(polarity_df,
                            .data$Polarity == 'LC',
                            .data$Class %in% log2FC_df$Class)
    lc_colors <- class_colors[which(names(class_colors) %in% lc_polarity_df$Class)]
    fia_polarity_df <- filter(polarity_df,
                                .data$Polarity == 'FIA',
                                .data$Class %in% log2FC_df$Class)
    fia_colors <- class_colors[which(names(class_colors) %in% fia_polarity_df$Class)]
    
    ## Legend: Manage breaks and values
    breaks <- c(names(lc_colors), names(fia_colors))
    values <- c(lc_colors,fia_colors)
    names(values) <- NULL
    
    ## Data: Replace NAs
    log2FC_df$log2FC[is.na(log2FC_df$log2FC)] <- 0
    log2FC_df$qval[is.na(log2FC_df$qval)] <- 1
    
    ## Data: Add color to data based on significance
    log2FC_df$signif_color <- sapply(log2FC_df$qval, function(q_val) {
        for (t in signif_colors) {
        if (q_val <= t) {
            color <- names(signif_colors)[which(signif_colors == t)]
        }
        }
        return(color)
    })

    ## Data: Add pseudo x-value to data as a order of metabolites
    ordered_classes <- c(names(lc_colors), names(fia_colors))
    p_data <- lapply(ordered_classes, function(class) {
        log2FC_df %>%
        filter(.data$Class == class) %>%
        bind_rows(data.frame(Class = rep(NA, 5)))
    }) %>%
        bind_rows()
    p_data <- bind_rows(data.frame(Class = rep(NA, 5)), p_data)
    p_data$x <- seq(nrow(p_data))
    p_data <- filter(p_data, !is.na(.data$Class))
    
    ## Legend: Significance color
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
    
    ## Background: Create data for background rects
    rects_df <- p_data %>%
        group_by(.data$Class) %>%
        summarise(Start = min(.data$x)-1,
                End = max(.data$x)+1,
                Color = class_colors[unique(.data$Class)],
                n = n())
    rects_df$Class <- factor(rects_df$Class, levels = breaks)
    rects_df$Technique <- sapply(rects_df$Class, function(c) {
        if (c %in% names(lc_colors)) {
        technique <- 'LC'
        } else if (c %in% names(fia_colors)) {
        technique <- 'FIA'
        } else {
        technique <- NA
        }
        return(technique)
    })
    
    ## Background: Determine border line between last LC and first FIA class
    lc_fia_border <- p_data %>%
        filter(.data$Class %in% names(lc_colors)) %>%
        select(.data$x) %>%
        max()

    ylims <- c(min(log2FC_df$log2FC) - 0.75, max(log2FC_df$log2FC) + 0,75)
    
    ## Plot: Create ggplot object
    p_scatter <- ggplot(p_data,
                        aes(x = .data$x,
                            y = .data$log2FC,
                            color = .data$signif_color)) +
        geom_rect(data = rects_df,
                inherit.aes = FALSE,
                aes(xmin = .data$Start, xmax = .data$End,
                    ymin = ylims[1], ymax = ylims[2],
                    fill = .data$Class,
                    text = paste0(.data$Class, "\n", .data$Technique, "\nNumber of Metabolites: ", n)),
                show.legend = TRUE,
                alpha = 0.4) +
        geom_vline(xintercept = 0, linewidth = 0.5, color = 'black') +
        geom_vline(xintercept = lc_fia_border+3, linewidth = 0.5, color = 'black', linetype="dotted") +
        geom_hline(yintercept = 0, linewidth = 0.5, color = 'black') +
        geom_point(size = 0.5, aes(text = paste0(.data$Metabolite, 
                                                "\nClass: ", .data$Class,
                                                "\nlog2 Fold Change: ", round(.data$log2FC, digits=2),
                                                "\nadj. p-value: ", round(.data$qval, digits=4),
                                                "\np-value: ", round(pval, digits=4)))) +
        scale_color_manual(paste0('Significance\n(linear model fit with FDR correction)'),
                        labels = signif_labels,
                        values = names(signif_colors),
                        guide = guide_legend(order=1)) +
        scale_fill_manual('Classes',
                        breaks = breaks,
                        values = values,
                        drop = FALSE,
                        guide = guide_legend(override.aes = list(alpha = 0.5),
                                            order=2, ncol = 2)) +
        theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            legend.key = element_rect(fill = 'white'),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line('#ECECEC'),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_line('#ECECEC'),
            panel.background = element_blank()) +
        labs(x = 'Metabolites')

    ## Interactive: Create interactive plot
    plotly_plot <- ggplotly(p_scatter, tooltip = "text")

    return(plotly_plot)
}
#' @title Plotly Log2FC Vulcano Plot

#'
#' @description This function returns a list with interactive 
#' vulcanoplot based on log2 fold change data.
#' 
#' @param metalyzer_se A Metalyzer object
#' @param cutoff_y A numeric value specifying the cutoff for q-value
#' @param cutoff_x A numeric value specifying the cutoff for log2 fold change
#' @param class_colors A csv file containing class colors hexcodes
#' 
#' @return plotly object
#'
#' @import dplyr
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom rlang .data
#' @importFrom plotly ggplotly
#' @export
#'
#' @examples
#' 
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_mutation_data_xl())
#' metalyzer_se <- filterMetabolites(
#'   metalyzer_se,
#'   drop_metabolites = "Metabolism Indicators"
#' )
#' metalyzer_se <- renameMetaData(
#'   metalyzer_se,
#'   Mutant_Control = "Sample Description"
#' )
#' metalyzer_se <- calculate_log2FC(
#'   metalyzer_se,
#'   categorical = "Mutant_Control",
#'   impute_perc_of_min = 0.2,
#'   impute_NA = TRUE
#' )
#' 
#' p_vulcano <- plotly_vulcano(metalyzer_se, 
#'                        cutoff_y = 0.05,
#'                        cutoff_x = 1.5)
#' 
plotly_vulcano <- function(metalyzer_se,
    cutoff_y = 0.05,
    cutoff_x = 1.5,
    class_colors = metalyzer_colors()) {

    log2FC_df <- metalyzer_se@metadata$log2FC                    
    # Make Colors unique for each class
    polarity_file <- system.file("extdata", "polarity.csv", package = "MetAlyzer")
    polarity_df <- utils::read.csv(polarity_file) %>%
    select(.data$Class,
            .data$Polarity) %>%
    mutate(Class = factor(.data$Class),
            Polarity = factor(.data$Polarity, levels = c('LC', 'FIA'))) %>%
    arrange(.data$Polarity)
    
    names(class_colors) <- levels(polarity_df$Class)
    ## Data: Replace NAs
    log2FC_df$log2FC[is.na(log2FC_df$log2FC)] <- 0
    log2FC_df$qval[is.na(log2FC_df$qval)] <- 1
    
    # Data Vulcano: Prepare Dataframe for vulcano plot
    log2FC_df$Class <- log2FC_df$Class
    log2FC_df$ClassColor <- log2FC_df$Class
    log2FC_df$ClassColor[log2FC_df$qval > cutoff_y] <- "Not Significant"
    log2FC_df$ClassColor[abs(log2FC_df$log2FC) < cutoff_x] <- "Not Significant"
    
    breaks <- unique(log2FC_df$Class)
    values <- class_colors[names(class_colors) %in% log2FC_df$Class]
    
    ## Plot: Create vulcano ggplot object
    p_fc_vulcano <- ggplot(log2FC_df,
                            aes(x = .data$log2FC,
                                y = -log10(.data$qval),
                                color = .data$ClassColor)) +
        geom_vline(xintercept=c(-log2(cutoff_x), log2(cutoff_x)), col="black", linetype="dashed") +
        geom_hline(yintercept=-log10(cutoff_y), col="black", linetype="dashed") +
        geom_point(size = 1, aes(text = paste0(.data$Metabolite, 
                                            "\nClass: ", .data$Class,
                                            "\nlog2 Fold Change: ", round(.data$log2FC, digits=2),
                                            "\nadj. p-value: ", round(.data$qval, digits=4),
                                            "\np-value: ", round(pval, digits=4)))) +
        scale_color_manual('Classes',
                            breaks = breaks,
                            values = values,
                            drop = FALSE) +
        theme_bw() +
        labs(x = 'log2(FC)', y = "-log10(p)")
    
    ## Interactive: Create interactive plot
    p_vulcano <- ggplotly(p_fc_vulcano, tooltip = "text")
    
    return(p_vulcano)
}