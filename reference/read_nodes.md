# Read and Validate Network Nodes (Metabolites)

Reads node data from a specified named region, validates entries against
pathway information, cleans labels, removes invalid nodes, and sets row
names.

## Usage

``` r
read_nodes(network_file, pathways, region_name = "Metabolites_Header")
```

## Arguments

- network_file:

  Path to the input file containing node data.

- pathways:

  A data frame of validated pathways (output of read_pathways).

- region_name:

  The named region or sheet containing metabolite header info.

## Value

A data frame of validated nodes with labels as row names.
