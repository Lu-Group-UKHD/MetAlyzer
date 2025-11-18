# Export filtered raw data as csv

This function exports the filtered raw data in the CSV format.

## Usage

``` r
export_conc_values(metalyzer_se, ..., file_path = "metabolomics_data.csv")
```

## Arguments

- metalyzer_se:

  SummarizedExperiment

- ...:

  Additional columns from meta_data

- file_path:

  file path

## Examples

``` r
metalyzer_se <- MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#> Error in MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates()): unused argument (conc_file_path = MetAlyzer::load_demodata_biocrates())

output_file <- file.path(tempdir(), "metabolomics_data.csv")
MetAlyzer::export_conc_values(metalyzer_se,
                              `Sample Description`,
                              file_path = output_file
                              )
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error in evaluating the argument 'x' in selecting a method for function 'colData': object 'metalyzer_se' not found
unlink(output_file)
```
