# Read and Validate Pathway Annotations

Reads pathway data from a specified named region in the pathway file,
validates entries, removes invalid ones, and sets row names.

## Usage

``` r
read_pathways(network_file, region_name = "Pathways_Header")
```

## Arguments

- network_file:

  Path to the input file containing pathway data.

- region_name:

  The named region or sheet containing pathway header info.

## Value

A data frame of validated pathway annotations with labels as row names.
