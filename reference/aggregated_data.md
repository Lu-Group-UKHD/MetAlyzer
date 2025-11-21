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

MetAlyzer::aggregated_data(metalyzer_se)
#> # A tibble: 20,380 × 5
#> # Groups:   Metabolite [1,019]
#>    ID    Metabolite Class          Concentration Status
#>    <fct> <fct>      <fct>                  <dbl> <fct> 
#>  1 14    C0         Acylcarnitines          16.4 Valid 
#>  2 15    C0         Acylcarnitines          15   Valid 
#>  3 16    C0         Acylcarnitines          14.9 Valid 
#>  4 17    C0         Acylcarnitines          16.4 Valid 
#>  5 18    C0         Acylcarnitines          16.2 Valid 
#>  6 19    C0         Acylcarnitines          15.9 Valid 
#>  7 20    C0         Acylcarnitines          15.3 Valid 
#>  8 21    C0         Acylcarnitines          16.2 Valid 
#>  9 22    C0         Acylcarnitines          15.2 Valid 
#> 10 23    C0         Acylcarnitines          16.3 Valid 
#> # ℹ 20,370 more rows
```
