# Differential analysis

Perform differential analysis and add the results table to the input
\`SummarizedExperiment\` (SE) object. The analysis is conducted using
the limma package, yielding log2 fold changes, p-values, and adjusted
p-values. Note that the input data must already be log2 transformed, and
should be subset if the variable of interest contains more than two
groups.

## Usage

``` r
calc_log2FC(metalyzer_se, group, group_level = NULL)
```

## Arguments

- metalyzer_se:

  An SE object output from [`read_webidq()`](read_webidq.md).

- group:

  A character specifying the sample metadata column containing two
  groups that will be compared.

- group_level:

  A length-2 vector of characters specifying the group members in
  \`group\`, which decides the direction of comparisons. For example,
  c('A', 'B') compares Group A to B, and vice versa. Default is NULL
  (alphabetically).

## Value

An SE object with an added table of differential analysis results,
accessible via \`metalyzer_se@metadata\$log2FC\`.
