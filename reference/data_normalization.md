# Normalization

Normalize concentration values among samples using glog2 transformation,
median normalization, or total ion count (TIC) normalization and update
the \`Concentration\` column in the aggregated data in the input
\`SummarizedExperiment\` (SE) object.

## Usage

``` r
data_normalization(metalyzer_se, norm_method = "log2")
```

## Arguments

- metalyzer_se:

  An SE object output from [`read_webidq`](read_webidq.md).

- norm_method:

  A character specifying the normalization method to use, which should
  be one of 'log2' (default), 'median', or 'TIC'.

## Value

An SE object with the normalized \`Concentration\` column in the
aggregated data accessible via
\`metalyzer_se@metadata\$aggregated_data\`.
