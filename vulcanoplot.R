
## Data: only color classes that are significantly differentially expressed
log2FC_vulcano$Class <- as.character(log2FC_vulcano$Class)
log2FC_vulcano$Class[log2FC_vulcano$qval > 0.05] <- NA
log2FC_vulcano$Class[abs(log2FC_vulcano$log2FC) < log2(1.5)] <- NA

## Data: Determine labels
log2FC_vulcano$labels <- as.character(log2FC_vulcano$Metabolite)
log2FC_vulcano$labels[which(is.na(log2FC_vulcano$Class))] <- ""

## Plot: Create ggplot object
p_fc <- ggplot(log2FC_vulcano,
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
p

#htmlwidgets::saveWidget(p, file = 'example_plot.html', selfcontained = FALSE, libdir= "lib/")
