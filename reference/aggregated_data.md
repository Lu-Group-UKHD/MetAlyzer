# Get Aggregated Data

This function returns the tibble "aggregated_data".

## Usage

``` r
aggregated_data(metalyzer_se)
```

## Arguments

- metalyzer_se:

  SummarizedExperiment

## Examples

``` r
metalyzer_se <- MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#> Error in MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates()): unused argument (conc_file_path = MetAlyzer::load_demodata_biocrates())

MetAlyzer::aggregated_data(metalyzer_se)
#> Error: object 'metalyzer_se' not found
```
