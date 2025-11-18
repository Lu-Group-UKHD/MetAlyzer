# Expand aggregated_data

Add sample metadata into the aggregated_data in the input
\`SummarizedExperiment\` (SE) object.

## Usage

``` r
expand_aggregated_data(metalyzer_se, metadata_col)
```

## Arguments

- metalyzer_se:

  An SE object output from [`read_webidq()`](read_webidq.md).

- metadata_col:

  A vector of character(s) specifying the column in the sample metadata
  to add, accessible via
  \`SummarizedExperiment::colData(metalyzer_se)\`.

## Value

An SE object with the added sample metadata in the aggregated data
accessible via \`metalyzer_se@metadata\$aggregated_data\`.
