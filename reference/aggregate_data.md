# Aggregate data

This function reshapes conc_values, quant_status, metatabolites and
sample IDs and combines them into a tibble data frame for filtering with
dplyr and plotting with 'ggplot2'. "aggregated_data" is grouped by
metabolites.

## Usage

``` r
aggregate_data(metabolites, meta_data, conc_values, quant_status, status_vec)
```

## Arguments

- metabolites:

  metabolites MetAlyzer object

- meta_data:

  Meta_data of the MetAlyzer object

- conc_values:

  conc_values of a MetAlyzer object

- quant_status:

  quant_status of a MetAlyzer object

- status_vec:

  A vector of quantification status
