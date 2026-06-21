# Calculate Top Loading Gene Expression Shifts

This function identifies genes with the highest loadings for specified
principal components and performs statistical tests to detect
distributional differences between query and reference data. It also
calculates the proportion of variance explained by each principal
component within specific cell types. Optionally, it can detect
anomalous cells using isolation forests.

This function creates visualizations showing expression distributions
for top loading genes that exhibit distributional differences between
query and reference datasets. Can display results as elegant complex
heatmaps, information-rich summary boxplots, or pseudo-bulk fold change
barplots. Optionally displays anomaly status when available.

## Usage

``` r
calculateGeneShifts(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  cell_types = NULL,
  pc_subset = 1:5,
  n_top_loadings = 50,
  genes_to_analyze = NULL,
  p_value_threshold = 0.05,
  adjust_method = "fdr",
  assay_name = "logcounts",
  detect_anomalies = FALSE,
  anomaly_comparison = FALSE,
  threshold_method = c("MAD", "absolute"),
  mad_multiplier = 2,
  anomaly_threshold = 0.5,
  n_tree = 500,
  max_cells_query = 5000,
  max_cells_ref = 5000
)

# S3 method for class 'calculateGeneShiftsObject'
plot(
  x,
  cell_type,
  pc_subset = 1:3,
  plot_type = c("heatmap", "barplot", "boxplot"),
  plot_by = c("p_adjusted", "top_loading"),
  n_genes = 10,
  significance_threshold = 0.05,
  show_anomalies = FALSE,
  pseudo_bulk = FALSE,
  cluster_cols = FALSE,
  draw_plot = TRUE,
  show_all_query = TRUE,
  max_cells_ref = NULL,
  max_cells_query = NULL,
  ...
)
```

## Arguments

- query_data:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing numeric expression matrix for the query cells.

- reference_data:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing numeric expression matrix for the reference cells.

- query_cell_type_col:

  The column name in the `colData` of `query_data` that identifies the
  cell types.

- ref_cell_type_col:

  The column name in the `colData` of `reference_data` that identifies
  the cell types.

- cell_types:

  A character vector specifying the cell types to analyze. If NULL, all
  common cell types are used.

- pc_subset:

  A numeric vector specifying which principal components to plot.
  Default is 1:3.

- n_top_loadings:

  Number of top loading genes to analyze per PC. Default is 50.

- genes_to_analyze:

  A character vector specifying genes to analyze. If NULL (default),
  genes are selected based on top loadings from specified principal
  components (see `n_top_loadings`). Default is NULL.

- p_value_threshold:

  P-value threshold for statistical significance. Default is 0.05.

- adjust_method:

  Method for multiple testing correction. Default is "fdr".

- assay_name:

  Name of the assay on which to perform computations. Default is
  "logcounts".

- detect_anomalies:

  Logical indicating whether to perform anomaly detection using
  isolation forests. Default is FALSE.

- anomaly_comparison:

  Logical indicating whether to perform statistical comparisons between
  non-anomalous reference cells and anomalous query cells instead of
  all-vs-all comparisons. When TRUE, only non-anomalous reference cells
  are compared against only anomalous query cells for each cell type.
  Requires detect_anomalies = TRUE. Default is FALSE.

- threshold_method:

  A character string specifying the method to determine anomaly cutoffs.
  Options are `"MAD"` (Median Absolute Deviation) or `"absolute"`.
  Default is `"MAD"`.

- mad_multiplier:

  A numeric value specifying the number of MADs above the reference
  median to use as the cutoff when `threshold_method = "MAD"`. Default
  is 2.

- anomaly_threshold:

  A numeric value specifying the absolute threshold for identifying
  anomalies when `threshold_method = "absolute"`. Default is 0.5.

- n_tree:

  An integer specifying the number of trees for the isolation forest
  when `detect_anomalies` is TRUE. Default is 500.

- max_cells_query:

  Maximum number of query cells to include in the plot. If NULL, all
  available query cells are plotted. Default is NULL.

- max_cells_ref:

  Maximum number of reference cells to include in the plot. If NULL, all
  available reference cells are plotted. Default is NULL.

- x:

  An object of class `calculateGeneShiftsObject`.

- cell_type:

  A character string specifying the cell type to plot (must be exactly
  one).

- plot_type:

  A character string specifying visualization type. Either "heatmap",
  "barplot", or "boxplot". Default is "heatmap".

- plot_by:

  A character string specifying gene selection method when \`n_genes\`
  is not NULL. Either "top_loading" or "p_adjusted". Default is
  "p_adjusted".

- n_genes:

  Number of top genes to show per PC. Can be NULL if
  \`significance_threshold\` is set. Default is 10.

- significance_threshold:

  If not NULL, a numeric value between 0 and 1. Used for gene selection
  or annotation. Default is 0.05.

- show_anomalies:

  Logical indicating whether to display anomaly status annotations.
  Default is FALSE. Requires anomaly results to be present in the
  object.

- pseudo_bulk:

  Logical indicating whether to create pseudo-bulk profiles instead of
  showing individual cells. When TRUE, expression values are averaged
  within groups (dataset and optionally anomaly status). Not compatible
  with boxplot visualization. Required for barplot visualization.
  Default is FALSE.

- cluster_cols:

  Logical indicating whether to cluster columns in the heatmap when
  \`pseudo_bulk = TRUE\`. When TRUE, columns (pseudo-bulk profiles) will
  be hierarchically clustered. When FALSE, columns maintain their
  original ordering (Query groups followed by Reference groups). Only
  applicable when \`pseudo_bulk = TRUE\` and \`plot_type = "heatmap"\`.
  Default is FALSE.

- draw_plot:

  Logical indicating whether to draw the plot immediately (TRUE) or
  return the undrawn plot object (FALSE). For heatmaps, FALSE returns a
  ComplexHeatmap object that can be further customized before drawing.
  Default is TRUE.

- show_all_query:

  Logical indicating whether to show the yellow bar for all query vs
  reference comparison. Default is TRUE. When FALSE, only green and red
  bars are shown.

- ...:

  Additional arguments passed to
  [`draw`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html) or
  not used for other plot types.

## Value

A list containing:

- PC results: Named elements for each PC (e.g., "PC1", "PC2") containing
  data frames with gene-level analysis results.

- expression_data: Matrix of expression values for all analyzed genes
  (genes × cells).

- cell_metadata: Data frame with columns: cell_id, dataset, cell_type,
  original_index, and optionally anomaly_status.

- gene_metadata: Data frame with columns: gene, pc, loading for all
  analyzed genes.

- percent_var: Named numeric vector of global percent variance explained
  for each analyzed PC.

- cell_type_variance: A data frame detailing the percent of variance a
  global PC explains within specific cell types for both query and
  reference datasets.

- anomaly_results: If `detect_anomalies` is TRUE, contains the full
  output from `detectAnomaly`.

A plot object. For heatmaps when `draw_plot = FALSE`, returns a
ComplexHeatmap object. For boxplots and barplots, returns a ggplot2
object.

## Details

This function extracts the top loading genes for each specified
principal component from the reference PCA space and performs
distributional comparisons between query and reference data. For each
gene, it performs statistical tests to identify genes that may be
causing PC-specific alignment issues between datasets. A key feature is
the calculation of cell-type-specific variance explained by global PCs,
providing a more nuanced view of how major biological axes affect
individual populations. When anomaly detection is enabled, isolation
forests are used to identify anomalous cells based on their PCA
projections.

When `anomaly_comparison = TRUE`, the statistical analysis focuses
specifically on comparing non-anomalous reference cells against
anomalous query cells. This can help identify genes that are
differentially expressed between "normal" reference cells and
potentially problematic query cells, providing insights into what makes
certain query cells anomalous.

This function visualizes the results from `calculateGeneShifts`. The
"heatmap" option displays a hierarchically clustered set of genes. The
"boxplot" option creates a two-panel plot using \`ggplot2\`: the left
panel shows horizontal expression boxplots for up to 5 PCs, while the
right panel displays their corresponding PC loadings and adjusted
p-values. The "barplot" option creates horizontal barplots showing log2
fold changes between pseudo-bulk expression profiles (query vs
reference), with genes ordered identically to the heatmap clustering.
Bars show comparisons for query non-anomaly (green), optionally all
query cells (yellow), and query anomaly cells (red) versus reference.
When anomaly detection results are available and `show_anomalies` is
TRUE, additional annotation bars or visual cues highlight anomalous
cells.

## See also

`plot.calculateGeneShiftsObject`,
[`detectAnomaly`](https://ccb-hms.github.io/scDiagnostics/reference/detectAnomaly.md)

`calculateGeneShifts`

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
