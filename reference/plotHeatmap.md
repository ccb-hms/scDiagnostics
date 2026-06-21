# Plot Heatmaps for Top Loading Gene Shifts (Simplified Single Heatmap)

This internal helper function creates a single, hierarchically clustered
heatmap displaying expression data for top loading genes from principal
component analysis. The function handles gene selection, data
preprocessing, and visualization formatting. Optionally includes anomaly
status annotations with proper cell ordering.

## Usage

``` r
plotHeatmap(
  x,
  cell_type,
  available_pcs,
  plot_by,
  n_genes,
  significance_threshold,
  show_anomalies,
  pseudo_bulk = FALSE,
  cluster_cols = FALSE
)
```

## Arguments

- x:

  An object of class `calculateGeneShiftsObject` containing expression
  data and analysis results.

- cell_type:

  A character string specifying the cell type to visualize.

- available_pcs:

  A character vector of principal components to include in analysis.

- plot_by:

  A character string indicating gene selection criterion ("top_loading"
  or "p_adjusted").

- n_genes:

  An integer specifying the number of top genes to display per PC. Can
  be NULL.

- significance_threshold:

  A numeric value between 0 and 1 for significance filtering. Can be
  NULL.

- show_anomalies:

  Logical indicating whether to show anomaly annotations.

- pseudo_bulk:

  Logical indicating whether to create pseudo-bulk profiles instead of
  showing individual cells. When TRUE, expression values are averaged
  within groups (dataset and optionally anomaly status). Not compatible
  with boxplot visualization. Default is FALSE.

- cluster_cols:

  Logical indicating whether to cluster columns in the heatmap when
  \`pseudo_bulk = TRUE\`. When TRUE, columns (pseudo-bulk profiles) will
  be hierarchically clustered. When FALSE, columns maintain their
  original ordering (Query groups followed by Reference groups). Only
  applicable when \`pseudo_bulk = TRUE\` and \`plot_type = "heatmap"\`.
  Default is FALSE.

## Value

A ComplexHeatmap object ready for plotting, or NULL if no genes meet
selection criteria.

## Details

This function generates a ComplexHeatmap visualization showing scaled
gene expression data across cells. Genes are selected based on either
their loading values or statistical significance. The function
automatically handles data filtering, scaling, and visual formatting
including color schemes and annotations.

Cell ordering (left to right): 1. Query Anomalous cells (leftmost) 2.
Query Normal cells 3. Reference Anomalous cells 4. Reference Normal
cells (rightmost)

This ensures clear dataset separation while grouping anomalous cells
together within each dataset for easy visual identification.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
