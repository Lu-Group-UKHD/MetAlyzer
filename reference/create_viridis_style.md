# Creates a viridis color style for Plotly plots.

This function can generate either a Plotly-compatible colorscale for a
color bar or a vector of hex color codes for manual coloring.

## Usage

``` r
create_viridis_style(
  color_scale,
  type = "scale",
  data = NULL,
  values_col_name = NULL
)
```

## Arguments

- color_scale:

  The name of the palette (e.g., "Magma", "Viridis").

- type:

  The desired output type: "scale" (for a color bar), "hex" or "inital"
  (for a color scalde, a vector of hex codes, or the correct scale for
  viridis package). Defaults to "scale".

- data:

  The data frame containing the values. Only required if type = "hex".

- values_col_name:

  The name of the column with numeric values. Only required if type =
  "hex".

## Value

A data frame if type is "scale", or a character vector if type is "hex".
