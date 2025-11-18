# Get log2FC Data

This function returns the tibble "log2FC".

## Usage

``` r
log2FC(metalyzer_se)
```

## Arguments

- metalyzer_se:

  SummarizedExperiment

## Examples

``` r
metalyzer_se <- MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#> Error in MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates()): unused argument (conc_file_path = MetAlyzer::load_demodata_biocrates())
metalyzer_se@metadata$log2FC <- readRDS(MetAlyzer::toy_diffres())
#> Error: object 'metalyzer_se' not found
MetAlyzer::log2FC(metalyzer_se)
#> Error: object 'metalyzer_se' not found
```
