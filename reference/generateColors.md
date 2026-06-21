# Generate Paired Colors for Cell Types

This function assigns paired colors (light and dark) to a list of cell
type names. The colors are selected from various color palettes in the
\`pals\` package.

## Usage

``` r
generateColors(cell_type_names, paired = FALSE)
```

## Arguments

- cell_type_names:

  A character vector of cell type names that need to be assigned colors.

- paired:

  If TRUE, the colored returned should be paired. Default is FALSE.

## Value

A named character vector where the names are the original cell type
names, and the values are the assigned colors.

## Details

The function uses color palettes from the \`pals\` package to generate
colors or pairs of colors (light and dark) for each cell type name
provided. It cycles through different color families (blues, greens,
reds, oranges, purples, purd and greys) to create the colors

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
