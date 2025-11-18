# Scatter Plot Visualization

This method creates a scatter plot of the log2 fold change for each
metabolite.

## Usage

``` r
plot_scatter(
  log2fc_df,
  show_labels_for = NULL,
  values_col_name = "log2FC",
  stat_col_name = "qval",
  show_p_value = TRUE,
  signif_colors = c(`#5F5F5F` = 1, `#FEBF6E` = 0.1, `#EE5C42` = 0.05, `#8B1A1A` = 0.01),
  save_as = NULL,
  folder_name = format(Sys.Date(), "%Y-%m-%d"),
  folder_path = NULL,
  file_name = "network",
  format = "pdf",
  width = 29.7,
  height = 21,
  units = "cm",
  overwrite = FALSE
)
```

## Arguments

- log2fc_df:

  DF with metabolites as row names and columns including log2FC, Class,
  qval columns.

- show_labels_for:

  Vector with Strings of Metabolite names or classes.

- values_col_name:

  Column name of a column that holds numeric values, to be plotted
  **Default = "log2FC"**

- stat_col_name:

  Columnname that holds numeric stat values that are used for
  significance **Default = "qval"**

- show_p_value:

  Boolean Value, to color p-values according to their significance level
  and add a Legend **Default = TRUE**

- signif_colors:

  Vector assigning significance values different colors

- save_as:

  *Optional:* Select the file type of output plots. Options are svg,
  pdf, png or NULL. **Default = "NULL"**

- folder_name:

  Name of the folder where the plot will be saved. Special characters
  will be removed automatically. **Default = date**

- folder_path:

  *Optional:* User-defined path where the folder should be created. If
  not provided, results will be saved in \`MetAlyzer_results\` within
  the working directory. **Default = NULL**

- file_name:

  Name of the output file (without extension). **Default = "network"**

- format:

  File format for saving the plot (e.g., "png", "pdf", "svg"). **Default
  = "pdf"**

- width:

  Width of the saved plot in specified units. **Default = 29.7**

- height:

  Height of the saved plot in specified units. **Default = 21.0**

- units:

  Units for width and height (e.g., "in", "cm", "mm"). **Default =
  "cm"**

- overwrite:

  Logical: If \`TRUE\`, overwrite existing files without asking. If
  \`FALSE\`, prompt user before overwriting. **Default = FALSE**

## Value

ggplot object

## Examples

``` r
log2fc_df <- readRDS(MetAlyzer::toy_diffres())
scatter <- MetAlyzer::plot_scatter(log2fc_df)
```
