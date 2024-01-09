Log2FCTab <- calculate_log2FC(metalyzer_se, Method, impute_perc_of_min = 0.2, impute_NA = TRUE)

### Server Side Data Wrangling
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
                         .data$Class %in% Log2FCTab$Class)
lc_colors <- class_colors[which(names(class_colors) %in% lc_polarity_df$Class)]
fia_polarity_df <- filter(polarity_df,
                          .data$Polarity == 'FIA',
                          .data$Class %in% Log2FCTab$Class)
fia_colors <- class_colors[which(names(class_colors) %in% fia_polarity_df$Class)]

## Data: Replace NAs
Log2FCTab$log2FC[is.na(Log2FCTab$log2FC)] <- 0
Log2FCTab$qval[is.na(Log2FCTab$qval)] <- 1


## Data: Add color to data based on significance
Log2FCTab$signif_color <- sapply(Log2FCTab$qval, function(q_val) {
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
  Log2FCTab %>%
    filter(.data$Class == class) %>%
    bind_rows(data.frame(Class = rep(NA, 5)))
}) %>%
  bind_rows()
p_data <- bind_rows(data.frame(Class = rep(NA, 5)), p_data)
p_data$x <- seq(nrow(p_data))
p_data <- filter(p_data, !is.na(.data$Class))

## Data: Determine labels
signif_p_data <- filter(p_data, .data$signif_color != names(signif_colors)[1])

labels <- sapply(p_data$Metabolite, function(m) {
  m <- as.character(m)
  label <- if_else(m %in% signif_p_data$Metabolite, m, "")
  return(label)
})

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

breaks <- c(names(lc_colors), names(fia_colors))
values <- c(lc_colors,fia_colors)
names(values) <- NULL

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

## Plot: Pre-plot to determine ylims
pre_plot <- ggplot(p_data,
                   aes(x = .data$x,
                       y = .data$log2FC,
                       color = .data$signif_color)) +
  geom_point(size = 0.5)

ylims <- ggplot_build(pre_plot)$layout$panel_params[[1]]$y.range

## Plot: Create ggplot object
p_fc <- ggplot(p_data,
               aes(x = .data$x,
                   y = .data$log2FC,
                   color = .data$signif_color)) +
  geom_rect(data = rects_df,
            inherit.aes = FALSE,
            aes(xmin = .data$Start, xmax = .data$End,
                ymin = ylims[1], ymax = ylims[2],
                fill = .data$Class,
                text = paste0(Class, "\n", Technique, "\nNumber of Metabolites: ", n)),
            show.legend = TRUE,
            alpha = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = 'black') +
  geom_vline(xintercept = lc_fia_border+3, linewidth = 0.5, color = 'black', linetype="dotted") +
  geom_hline(yintercept = 0, linewidth = 0.5, color = 'black') +
  geom_point(size = 0.5, aes(text = paste0(Metabolite, 
                                           "\nClass: ", Class, 
                                           "\nlog2 Fold Change: ", round(log2FC, digits=5),  
                                           "\np-value: ", round(pval, digits=5)))) + 
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
p <- ggplotly(p_fc, tooltip = "text", showlegend = FALSE)
p <- hide_legend(p)
p

# Grab Legend ggplot
tmp <- ggplot_gtable(ggplot_build(p_fc))
leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
legend <- tmp$grobs[[leg]]
legend <- grid.arrange(legend, ncol=1)
