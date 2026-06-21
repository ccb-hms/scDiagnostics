# Project Query Data Onto SIR Space of Reference Data

This function projects a query `SingleCellExperiment` object onto the
SIR (supervised independent component) space of a reference
`SingleCellExperiment` object. The SVD of the reference data is computed
on conditional means per cell type, and the query data is projected
based on these reference components.

## Usage

``` r
projectSIR(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  cell_types = NULL,
  multiple_cond_means = TRUE,
  cumulative_variance_threshold = 0.7,
  n_neighbor = 1,
  assay_name = "logcounts",
  max_cells_query = 5000,
  max_cells_ref = 5000
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

  A character string specifying the column in the `colData` of
  `query_data` that identifies the cell types.

- ref_cell_type_col:

  A character string specifying the column in the `colData` of
  `reference_data` that identifies the cell types.

- cell_types:

  A character vector of cell types for which to compute conditional
  means in the reference data.

- multiple_cond_means:

  A logical value indicating whether to compute multiple conditional
  means per cell type (through PCA and clustering). Defaults to `TRUE`.

- cumulative_variance_threshold:

  A numeric value between 0 and 1 specifying the variance threshold for
  PCA when computing multiple conditional means. Defaults to `0.7`.

- n_neighbor:

  An integer specifying the number of nearest neighbors for clustering
  when computing multiple conditional means. Defaults to `1`.

- assay_name:

  A character string specifying the assay name on which to perform
  computations. Defaults to `"logcounts"`.

- max_cells_query:

  Maximum number of query cells to retain after cell type filtering. If
  NULL, no downsampling of query cells is performed. Default is 5000.

- max_cells_ref:

  Maximum number of reference cells to retain after cell type filtering.
  If NULL, no downsampling of reference cells is performed. Default is
  5000.

## Value

A list containing:

- cond_means:

  A matrix of the conditional means computed for the reference data.

- rotation_mat:

  The rotation matrix obtained from the SVD of the conditional means.

- sir_projections:

  A `data.frame` containing the SIR projections for both the reference
  and query datasets.

- percent_var:

  The percentage of variance explained by each component of the SIR
  projection.

## Details

The genes used for the projection (SVD) must be present in both the
reference and query datasets. The function first computes conditional
means for each cell type in the reference data, then performs SVD on
these conditional means to obtain the rotation matrix used for
projecting both the reference and query datasets. The query data is
centered and scaled based on the reference data.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>

## Examples

``` r
# Load data
data("reference_data")
data("query_data")

# Project the query data onto SIR space of reference
sir_output <- projectSIR(query_data = query_data,
                         reference_data = reference_data,
                         query_cell_type_col = "SingleR_annotation",
                         ref_cell_type_col = "expert_annotation")
```
