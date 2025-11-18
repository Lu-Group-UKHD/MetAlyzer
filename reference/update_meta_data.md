# Update meta data

This function adds another column to filtered meta_data.

## Usage

``` r
update_meta_data(metalyzer_se, ..., inplace = FALSE)
```

## Arguments

- metalyzer_se:

  SummarizedExperiment

- ...:

  Use ´new_col_name = new_column´ to rename selected variables

- inplace:

  If FALSE, return a copy. Otherwise, do operation inplace and return
  None.

## Value

An updated SummarizedExperiment

## Examples

``` r
metalyzer_se <- MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#> Error in MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates()): unused argument (conc_file_path = MetAlyzer::load_demodata_biocrates())

metalyzer_se <- MetAlyzer::update_meta_data(
  metalyzer_se,
  Date = Sys.Date(), Analyzed = TRUE
)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colData': object 'metalyzer_se' not found
```
