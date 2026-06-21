# Generate Variance Contribution Bar Plot with Component Breakdown

Creates a bar plot visualization of variance contributions for each
principal component, showing how much total dataset variance is
explained by the regression model. When available, displays stacked bars
showing individual component contributions (cell type, batch/dataset,
and interaction effects).

## Usage

``` r
plotVarianceContribution(x, ...)
```

## Arguments

- x:

  An object of class `regressPCObject` containing regression results
  from `regressPC` function.

- ...:

  Additional arguments passed to the plotting function (currently
  unused).

## Value

A `ggplot2` object representing the variance contribution bar plot with
optional component breakdown.

## Details

This function visualizes the variance contribution of each principal
component, calculated as the product of PC variance and R-squared
values. The variance contribution represents the percentage of total
dataset variance explained by the regression model for each PC.

When component decomposition is available, the function creates stacked
bars with:

- Cell type main effect contribution (blue)

- Batch or dataset main effect contribution (orange)

- Interaction effect contribution (green)

The plot subtitle includes the total variance explained across all
principal components. PC labels show the individual variance percentage
for each component. The function automatically selects appropriate PCA
variance values based on analysis type (query-only vs. query+reference).

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
