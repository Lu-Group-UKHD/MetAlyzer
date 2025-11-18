# Get sample labels from metadata

Extracts sample labels for plotting. Looks for columns "Sample ID",
"Sample Identification", or "Identification" (case-insensitive, ignoring
punctuation). Falls back to ID column if not found.

## Usage

``` r
get_sample_labels(smpMetadatTbl)
```

## Arguments

- smpMetadatTbl:

  A data frame with an 'ID' column and optional label column.

## Value

Character vector of sample labels, same length as nrow(smpMetadatTbl).
