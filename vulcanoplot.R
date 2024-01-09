
metalyzer_se <- MetAlyzer_dataset(file_path = example_extraction_data())

metalyzer_se <- renameMetaData(metalyzer_se, Method = 'Sample Description')

log2FC_df <- calculate_log2FC(metalyzer_se, Method, impute_perc_of_min = 0.2, impute_NA = TRUE)

## Background: Define colors for significance 
signif_colors=c("#5F5F5F"=1,
                "#FEBF6E"=0.1,
                "#EE5C42"=0.05,
                "#8B1A1A"=0.01)

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

## Data: Replace NAs
log2FC_df$log2FC[is.na(log2FC_df$log2FC)] <- 0
log2FC_df$qval[is.na(log2FC_df$qval)] <- 1

## Data: only color classes that are significantly differentially expressed
log2FC_df$Class <- as.character(log2FC_df$Class)
log2FC_df$Class[log2FC_df$qval > 0.05] <- NA
log2FC_df$Class[abs(log2FC_df$log2FC) < log2(1.5)] <- NA

## Data: Determine labels
log2FC_df$labels <- as.character(log2FC_df$Metabolite)
log2FC_df$labels[which(is.na(log2FC_df$Class))] <- ""

## Update lc_colors and fia_colors
lc_colors <- lc_colors[which(names(lc_colors) %in% log2FC_df$Class)]
fia_colors <- fia_colors[which(names(fia_colors) %in% log2FC_df$Class)]

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

## Plot: Create ggplot object
p_fc <- ggplot(log2FC_df,
               aes(x = .data$log2FC,
                   y = -log10(.data$qval),
                   color = .data$Class,
                   label = labels)) +
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") +
  geom_point(size = 1, aes(text = paste0(Metabolite, "\nClass: ", Class, "\nlog2 Fold Change: ", round(log2FC, digits=5), "\np-value: ", round(pval, digits=5)))) +
  scale_color_manual('Classes',
                     breaks = breaks,
                     values = values,
                     drop = FALSE,
                     guide = guide_legend(override.aes = list(size = 2),
                                          order=2, ncol = 2)) +
  theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
        legend.key = element_rect(fill = 'white')) +
  labs(x = 'log2(FC)', y = "-log10(p)")

## Interactive: Create interactive plot
p <- ggplotly(p_fc, tooltip = "text")

#htmlwidgets::saveWidget(p, file = 'example_plot.html', selfcontained = FALSE, libdir= "lib/")
