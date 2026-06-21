# Extract Gene Order from ComplexHeatmap Object

Extracts the hierarchically clustered gene order from a ComplexHeatmap
Heatmap object.

## Usage

``` r
extractGeneOrder(heatmap_object)
```

## Arguments

- heatmap_object:

  A ComplexHeatmap Heatmap object (undrawn or drawn).

## Value

A character vector of gene names in the clustered order (top to bottom
as they appear in heatmap).

## Details

This function initializes a ComplexHeatmap object (without displaying
it) to perform hierarchical clustering, then extracts the resulting gene
order. This allows matching gene ordering between heatmap and other plot
types. The function uses a null graphics device to initialize the
heatmap clustering without actually drawing the plot.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
