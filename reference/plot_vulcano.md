# Vulcano plot - ggplot

Create a static vulcano plot using \`ggplot()\`.

## Usage

``` r
plot_vulcano(
  Log2FCTab,
  x_cutoff = 1.5,
  y_cutoff = 0.05,
  show_labels_for = NULL
)
```

## Arguments

- Log2FCTab:

  A data frame containing the differential analysis results table
  accessible via \`metalyzer_se@metadata\$log2FC\` where
  \`metalyzer_se\` is an SE object output from
  [`read_webidq()`](read_webidq.md) and has gone through
  \`calc_log2FC()\`.

- x_cutoff, y_cutoff:

  Numerical values specifying the cutoffs for log2 fold changes and
  q-values. Default is 1.5 and 0.05, respectively.

- show_labels_for:

  A vector of characters specifying the names of metabolites or classes
  to label. Note that the metabolites of a specified class will be
  labeled. Default is NULL.

## Value

A \`ggplot\` object.
