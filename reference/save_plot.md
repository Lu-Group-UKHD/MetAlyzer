# Save plots

This function saves a given ggplot object to a specified folder and file
format. It ensures that the folder structure exists and cleans the
folder name to remove special characters.

## Usage

``` r
save_plot(
  plot,
  folder_name = format(Sys.Date(), "%Y-%m-%d"),
  folder_path = NULL,
  file_name = "network",
  format = "pdf",
  units = "cm",
  height = 21,
  width = 29.7,
  overwrite = FALSE
)
```

## Arguments

- plot:

  A ggplot object to be saved.

- folder_name:

  Name of the folder where the plot will be saved. Special characters
  will be removed automatically. **Default = date**

- folder_path:

  *Optional:* User-defined path where the folder should be created. If
  not provided, results will be saved in \`MetAlyzer_results\` within
  the working directory. **Default = NULL**

- file_name:

  Name of the output file (without extension). **Default = "network"**

- format:

  File format for saving the plot (e.g., "png", "pdf", "svg"). **Default
  = "pdf"**

- units:

  Units for width and height (e.g., "in", "cm", "mm"). **Default =
  "cm"**

- height:

  Height of the saved plot in specified units. **Default = 21.0**

- width:

  Width of the saved plot in specified units. **Default = 29.7**

- overwrite:

  Logical: If \`TRUE\`, overwrite existing files without asking. If
  \`FALSE\`, prompt user before overwriting. **Default = FALSE**

## Value

The function does not return anything but saves the plot to the
specified directory.
