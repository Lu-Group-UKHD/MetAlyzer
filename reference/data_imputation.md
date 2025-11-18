# Half-minimum imputation

Impute NA concentrations using half-minimum (HM) and update the
\`Concentration\` column in the aggregated data in the input
\`SummarizedExperiment\` (SE) object. If all values of a feature are NA,
they stay NA.

## Usage

``` r
data_imputation(metalyzer_se)
```

## Arguments

- metalyzer_se:

  An SE object output from [`read_webidq()`](read_webidq.md).

## Value

An SE object with the imputed \`Concentration\` column in the aggregated
data accessible via \`metalyzer_se@metadata\$aggregated_data\`.
