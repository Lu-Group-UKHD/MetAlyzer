# Filter meta data

This function updates the "Filter" column in meta_data to filter out
samples.

## Usage

``` r
filter_meta_data(metalyzer_se, ..., inplace = FALSE)
```

## Arguments

- metalyzer_se:

  SummarizedExperiment

- ...:

  Use ´col_name´ and condition to filter selected variables.

- inplace:

  If FALSE, return a copy. Otherwise, do operation inplace and return
  None.

## Value

An updated SummarizedExperiment

## Examples

``` r
metalyzer_se <- MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#> Error in MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates()): unused argument (conc_file_path = MetAlyzer::load_demodata_biocrates())

metalyzer_se <- MetAlyzer::filter_meta_data(metalyzer_se, `Sample Description` %in% 1:6)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error in evaluating the argument 'x' in selecting a method for function 'colData': object 'metalyzer_se' not found
```
