# Plot Boxplots for Top Loading Gene Shifts (Two-Panel Summary using ggplot2)

This internal helper function creates a comprehensive two-panel summary
plot displaying gene expression distributions and principal component
loadings. The visualization uses ggplot2 faceting to create side-by-side
panels. Optionally includes anomaly status information using visual
cues.

## Usage

``` r
plotBoxplot(
  x,
  cell_type,
  available_pcs,
  plot_by,
  n_genes,
  significance_threshold,
  show_anomalies
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

  Logical indicating whether to show anomaly annotations using line type
  differences in boxplot borders.

## Value

A ggplot2 object ready for display, or NULL if no genes meet selection
criteria. The returned plot contains:

- Left panel: Expression boxplots with dataset-specific fill colors

- Right panel: PC loading scatter points with PC-specific colors/shapes

- Gene labels on y-axis with significance indicators (\*)

- P-value annotations on secondary y-axis

- Legend showing dataset, PC information, and anomaly status (if
  applicable)

## Details

This function generates a dual-panel ggplot2 visualization where the
left panel shows horizontal boxplots of gene expression distributions
comparing Reference and Query datasets, while the right panel displays
PC loading values as points with adjusted p-values. Gene selection is
based on the union of top genes across specified principal components.

When anomaly information is available and requested, anomalous cells are
distinguished using dashed boxplot borders (normal cells have solid
borders). This approach avoids color conflicts between dataset
identification (fill colors) and PC identification (point
colors/shapes).

Visual encoding:

- Dataset: Fill colors (Reference = blue, Query = red)

- Anomaly status: Line types (Normal = solid, Anomaly = dashed borders)

- PC identity: Point colors and shapes in loading panel

## See also

[`plotHeatmap`](https://ccb-hms.github.io/scDiagnostics/reference/plotHeatmap.md),
[`plot.calculateGeneShiftsObject`](https://ccb-hms.github.io/scDiagnostics/reference/calculateGeneShifts.md)

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
