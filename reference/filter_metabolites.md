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
metalyzer_se <- MetAlyzer::read_webidq(file_path = MetAlyzer::load_demodata_biocrates())
#> Input checks passed. Proceeding with reading webidq file.
#> 
#> 
#>  _____ ______   _______  _________  ________  ___           ___    ___ ________  _______   ________
#> |\   _ \  _   \|\  ___ \|\___   ___\\   __  \|\  \         |\  \  /  /|\_____  \|\  ___ \ |\   __  \
#> \ \  \\\__\ \  \ \   __/\|___ \  \_\ \  \|\  \ \  \        \ \  \/  / /\|___/  /\ \   __/|\ \  \|\  \
#>  \ \  \\|__| \  \ \  \_|/__  \ \  \ \ \   __  \ \  \        \ \    / /     /  / /\ \  \_|/_\ \   _  _\
#>   \ \  \    \ \  \ \  \_|\ \  \ \  \ \ \  \ \  \ \  \____    \/   / /     /  /_/__\ \  \_|\ \ \  \\  \| 
#>    \ \__\    \ \__\ \_______\  \ \__\ \ \__\ \__\ \_______\__/   / /     |\________\ \_______\ \__\\ _\ 
#>     \|__|     \|__|\|_______|   \|__|  \|__|\|__|\|_______|\____/ /       \|_______|\|_______|\|__|\|__|
#>                                                           \|____|/
#> 
#> 
#> Info: Reading color code "FFCBD2D7" as "#CBD2D7"
#> Info: Reading color code "FFB9DE83" as "#B9DE83"
#> Info: Reading color code "FFB9DE83" as "#B9DE83"
#> Info: Reading color code "FFA28BA3" as "#A28BA3"
#> Info: Reading color code "FFA28BA3" as "#A28BA3"
#> Info: Reading color code "FFB2D1DC" as "#B2D1DC"
#> Info: Reading color code "FF7FB2C5" as "#7FB2C5"
#> Info: Reading color code "FFB2D1DC" as "#B2D1DC"
#> Info: Reading color code "FF7FB2C5" as "#7FB2C5"
#> 
#> Measured concentration values:
#> ------------------------------
#>           0%          25%          50%          75%         100% 
#>     0.000000     0.286299     1.289381     6.308854 12522.000000 
#> 
#> NAs: 762 (3.74%)
#> 
#> 
#> Measured quantification status:
#> -------------------------------
#> Valid: 15419 (75.66%)
#> LOQ: 983 (4.82%)
#> LOD: 3978 (19.52%)
#> NAs: 0 (0%)
#> 

drop_metabolites <- c("C0", "C2", "C3", "Metabolism Indicators",
  inplace = TRUE
)
metalyzer_se <- MetAlyzer::filter_metabolites(metalyzer_se, drop_metabolites)
#> Info: 3 metabolites were removed!
```
