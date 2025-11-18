# Calculate Node-Level Aggregate Statistics (Conditional on Significance)

Calculates mean log2FC, p-value, and q-value for each node (Label),
prioritizing significant metabolites (qval \<= q_value). If none are
significant, uses all measured metabolites for the node. Adds results to
both dataframes.

## Usage

``` r
calculate_node_aggregates_conditional(
  nodes_sep_df,
  nodes_orig_df,
  q_value,
  stat_col_name,
  ...
)
```

## Arguments

- nodes_sep_df:

  Dataframe with metabolites separated (e.g., 'nodes_final'). Must
  contain Label, log2FC, qval.

- nodes_orig_df:

  Original dataframe with potentially semi-colon separated metabolites.
  Must contain Label.

- q_value:

  Significance threshold for q-values (e.g., 0.05).

- stat_col_name:

  p value column name

- values_col_name:

  plotted column

## Value

A list containing two dataframes: \$nodes_separated: Input nodes_sep_df
with 2 new columns: node_values, node_stat \$nodes: Input nodes_orig_df
with 2 new columns: node_values, node_stat
