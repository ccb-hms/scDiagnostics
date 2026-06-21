# Plot Barplots for Top Loading Gene Fold Changes (Pseudo-Bulk)

This internal helper function creates a barplot visualization showing
pseudo-bulk fold changes between different cell groups for top loading
genes.

## Usage

``` r
plotBarplot(
  x,
  cell_type,
  available_pcs,
  plot_by,
  n_genes,
  significance_threshold,
  show_anomalies,
  show_all_query = TRUE
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

  An integer specifying the number of top genes to display per PC.

- significance_threshold:

  A numeric value between 0 and 1 for significance annotation.

- show_anomalies:

  Logical indicating whether to show anomaly-related bars.

- show_all_query:

  Logical indicating whether to show the yellow bar for all query vs
  reference comparison. Default is TRUE. When FALSE, only green and red
  bars are shown.

## Value

A ggplot2 object ready for display, or NULL if no genes meet selection
criteria.

## Details

This function generates a ggplot2 barplot where genes are arranged
vertically in hierarchically clustered order (identical to heatmap gene
ordering), and horizontal bars show log2 fold changes for different
pseudo-bulk comparisons vs reference. The function creates an internal
heatmap object with the same parameters to extract the exact gene
clustering order, ensuring consistency between plot types.

Bar colors:

- Green: Query non-anomaly (pseudo-bulk) vs Reference (pseudo-bulk)

- Yellow: All Query (pseudo-bulk) vs Reference (pseudo-bulk) (shown when
  show_all_query = TRUE and anomaly data available)

- Red: Query anomaly (pseudo-bulk) vs Reference (pseudo-bulk)

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
