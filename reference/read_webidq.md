# Open file and read data

This function creates a SummarizedExperiment (SE) from the given
'webidq' output Excel sheet: metabolites (rowData), meta data (colData),
concentration data (assay), quantification status(assay) The column
"Sample Type" and the row "Class" are used as anchor cells in the Excel
sheet and are therefore a requirement.

## Usage

``` r
read_webidq(
  file_path,
  sheet = 1,
  status_list = list(Valid = c("#B9DE83", "#00CD66"), LOQ = c("#B2D1DC", "#7FB2C5",
    "#87CEEB"), LOD = c("#A28BA3", "#6A5ACD", "#BBA7B9"), `ISTD Out of Range` =
    c("#FFF099", "#FFFF33"), Invalid = "#FFFFCC", Incomplete = c("#CBD2D7", "#FFCCCC")),
  silent = FALSE
)
```

## Arguments

- file_path:

  A character specifying the file path to the Excel file.

- sheet:

  A numeric index specifying which sheet of the Excel file to use.

- status_list:

  A list of HEX color codes for each quantification status.

- silent:

  If TRUE, mute any print command.

## Value

A Summarized Experiment object

## Examples

``` r
Path <- MetAlyzer::load_demodata_biocrates()
metalyzer_se <- MetAlyzer::read_webidq(file_path = Path)
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
```
