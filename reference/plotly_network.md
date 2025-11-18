# Network diagram - ggplotly

This function returns a list with interactive networkplot based on log2
fold change data.

## Usage

``` r
plotly_network(
  Log2FCTab,
  q_value = 0.05,
  metabolite_col_name = "Metabolite",
  values_col_name = "log2FC",
  stat_col_name = "qval",
  exclude_pathways = NULL,
  metabolite_node_size = 11,
  connection_width = 1.25,
  pathway_text_size = 20,
  pathway_width = 10,
  plot_height = 800,
  color_scale = "viridis"
)
```

## Arguments

- Log2FCTab:

  A dataframe with log2FC, qval, additional columns

- q_value:

  The q-value threshold for significance

- metabolite_col_name:

  Columnname that holds the Metabolites

- values_col_name:

  Column name of a column that holds numeric values, to be plotted
  **Default = "log2FC"**

- stat_col_name:

  Columnname that holds numeric stat values that are used for
  significance **Default = "qval"**

- exclude_pathways:

  Pathway names that are exluded from plotting

- metabolite_node_size:

  The text size of metabolite nodes

- connection_width:

  The line width of connections between metabolites

- pathway_text_size:

  The text size of pathway annotations

- pathway_width:

  The line width of pathway-specific connection coloring

- plot_height:

  The height of the Plot in pixel \[px\]

- color_scale:

  A string specifying the color scale to use. Options include
  \`"viridis"\`, \`"plasma"\`, \`"magma"\`, \`"inferno"\`,
  \`"cividis"\`, \`"rocket"\`, \`"mako"\`, and \`"turbo"\`, which use
  the \`viridis\` color scales.
