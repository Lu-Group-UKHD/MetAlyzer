# Summarize quantification status

This function lists the number of each quantification status and its
percentage.

## Usage

``` r
summarize_quant_data(metalyzer_se)
```

## Arguments

- metalyzer_se:

  SummarizedExperiment

## Examples

``` r
metalyzer_se <- MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#> Error in MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates()): unused argument (conc_file_path = MetAlyzer::load_demodata_biocrates())

MetAlyzer::summarize_quant_data(metalyzer_se)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'assay': object 'metalyzer_se' not found
```
