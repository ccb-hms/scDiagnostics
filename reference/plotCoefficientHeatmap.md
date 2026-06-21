# Generate Regression Coefficients Heatmap

Creates a heatmap visualization of regression coefficients from
principal component regression analysis, organized by coefficient type
(cell type, batch, interaction) and annotated with significance
indicators.

## Usage

``` r
plotCoefficientHeatmap(x, alpha = 0.05, coefficients_include = NULL, ...)
```

## Arguments

- x:

  An object of class `regressPCObject` containing regression results
  from `regressPC` function.

- alpha:

  Numeric value specifying the significance threshold for p-value
  adjustment. Default is 0.05.

- coefficients_include:

  Character vector specifying which coefficient types to include.
  Options are `c("cell_type", "batch", "interaction")`. Default is
  `NULL`, which includes all available coefficient types.

- ...:

  Additional arguments passed to the plotting function (currently
  unused).

## Value

A `ggplot2` object representing the regression coefficients heatmap with
faceted organization and significance annotations.

## Details

This function generates a comprehensive heatmap showing regression
coefficients for each principal component and model term. The
visualization includes:

- Color-coded coefficient values (blue = negative, red = positive)

- Significance indicators (asterisks) for adjusted p-values below
  threshold

- Faceted organization by coefficient category (Cell Type, Batch,
  Interaction)

- Clean term labels with proper formatting and reference category
  information

The function handles different model types automatically:

- Simple cell type models: `PC ~ cell_type`

- Batch interaction models: `PC ~ cell_type * batch`

- Dataset interaction models: `PC ~ cell_type * dataset`

Term labels are cleaned and formatted for better readability, with
batch/dataset terms converted to consistent "Query Batch" terminology.
The plot includes comprehensive subtitle information showing model
specification, significance threshold, and reference categories.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
