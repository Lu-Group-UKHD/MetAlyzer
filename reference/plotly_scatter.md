# Scatter plot - ggplotly

Return a list containing an interactive scatter plot made using
\`ggplotly()\` and its static legend using \`ggplot()\`. The x-axis
represents metabolic classes, the y-axis shows log2 fold changes, and
each point corresponds to a metabolite.

## Usage

``` r
plotly_scatter(Log2FCTab)
```

## Arguments

- Log2FCTab:

  A data frame containing the differential analysis results table
  accessible via \`metalyzer_se@metadata\$log2FC\` where
  \`metalyzer_se\` is an SE object output from
  [`read_webidq`](read_webidq.md) and has gone through
  \`calc_log2FC()\`.

## Value

A list containing a \`plotly\` (plot) and \`ggplot\` (legend) object.
