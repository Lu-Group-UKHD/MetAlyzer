# Filter metabolites

This function filters out certain classes or metabolites of the
metabolites vector. If aggregated_data is not empty, metabolites and
class will also be filtered here.

## Usage

``` r
filter_metabolites(
  metalyzer_se,
  drop_metabolites = c("Metabolism Indicators"),
  drop_NA_concentration = FALSE,
  drop_quant_status = NULL,
  min_percent_valid = NULL,
  valid_status = c("Valid", "LOQ"),
  per_group = NULL,
  inplace = FALSE
)
```

## Arguments

- metalyzer_se:

  SummarizedExperiment

- drop_metabolites:

  A character vector defining metabolite classes or individual
  metabolites to be removed

- drop_NA_concentration:

  A boolean whether to drop metabolites which have any NAs in their
  concentration value

- drop_quant_status:

  A character, vector of characters or list of characters specifying
  which quantification status to remove. Metabolites with at least one
  quantification status of this vector will be removed.

- min_percent_valid:

  A numeric lower threshold between 0 and 1 (t less than or equal to x)
  to remove invalid metabolites that do not meet a given percentage of
  valid measurements per group (default per Metabolite).

- valid_status:

  A character vector that defines which quantification status is
  considered valid.

- per_group:

  A character vector of column names from meta_data that will be used to
  split each metabolite into groups. The threshold \`min_percent_valid\`
  will be applied for each group. The selected columns from meta_data
  will be added to aggregated_data.

- inplace:

  If FALSE, return a copy. Otherwise, do operation inplace and return
  None.

## Value

An updated SummarizedExperiment

## Examples

``` r
metalyzer_se <- MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates())
#> Error in MetAlyzer::read_webidq(conc_file_path = MetAlyzer::load_demodata_biocrates()): unused argument (conc_file_path = MetAlyzer::load_demodata_biocrates())

drop_metabolites <- c("C0", "C2", "C3", "Metabolism Indicators",
  inplace = TRUE
)
metalyzer_se <- MetAlyzer::filter_metabolites(metalyzer_se, drop_metabolites)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rowData': object 'metalyzer_se' not found
```
