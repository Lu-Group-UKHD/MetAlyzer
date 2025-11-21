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

metalyzer_se <- MetAlyzer::update_meta_data(
  metalyzer_se,
  Date = Sys.Date(), Analyzed = TRUE
)
```
