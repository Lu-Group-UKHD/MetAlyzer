# Slice metabolites

This function extracts metabolites with their corresponding metabolite
class from .full_sheet into metabolites. It robustly finds the
metabolite row by searching upwards from the "Class" row and skipping
any known intermediate rows (e.g., "Synonym").

## Usage

``` r
slice_metabolites(full_sheet, data_ranges)
```

## Arguments

- full_sheet:

  full_sheet

- data_ranges:

  data_ranges
