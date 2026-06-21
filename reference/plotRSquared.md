# Generate R-squared Bar Plot with Component Breakdown

Creates a bar plot visualization of R-squared values for each principal
component, with optional stacked bars showing the contribution of
individual model components (cell type, batch/dataset, and interaction
effects) when component decomposition is available.

## Usage

``` r
plotRSquared(x, ...)
```

## Arguments

- x:

  An object of class `regressPCObject` containing regression results
  from `regressPC` function.

- ...:

  Additional arguments passed to the plotting function (currently
  unused).

## Value

A `ggplot2` object representing the R-squared bar plot with optional
component breakdown.

## Details

This function generates either a simple bar plot or a stacked bar plot
depending on the availability of component-wise R-squared decomposition
in the input object. When component breakdown is available, the bars are
stacked to show:

- Cell type main effect (blue)

- Batch or dataset main effect (orange)

- Interaction effect (green)

Principal component labels include the percentage of total variance
explained by each PC. The total R-squared value is displayed above each
bar. For query-only analyses, the function uses query PCA variance; for
query+reference analyses, it uses reference PCA variance.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
