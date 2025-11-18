# Read and Validate Network Edges (Connections)

Reads edge data from a specified named region, validates that connected
nodes exist and are not self-loops, and removes invalid edges.

## Usage

``` r
read_edges(network_file, nodes, pathways, region_name = "Connections_Header")
```

## Arguments

- network_file:

  Path to the input file containing edge data.

- nodes:

  A data frame of validated nodes (output of read_nodes).

- pathways:

  A data frame of validated pathways (output of read_pathways).

- region_name:

  The named region or sheet containing connections header info.

## Value

A data frame of validated edges.
