# Argument Validation for SingleCellExperiment Analysis

This function validates the input arguments for functions that analyze
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
objects. It checks that the inputs are of the correct types and formats,
and that required columns and cell types are present in the data.

## Usage

``` r
argumentCheck(
  query_data = NULL,
  reference_data = NULL,
  query_cell_type_col = NULL,
  ref_cell_type_col = NULL,
  unique_cell_type = FALSE,
  plot_function = FALSE,
  cell_names_query = NULL,
  cell_names_ref = NULL,
  pc_subset_query = NULL,
  pc_subset_ref = NULL,
  common_rotation_genes = FALSE,
  assay_name = NULL,
  max_cells_ref = NULL,
  max_cells_query = NULL
)
```

## Arguments

- query_data:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing numeric expression matrix for the query cells. If
  \`NULL\`, no check is performed.

- reference_data:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing numeric expression matrix for the reference cells.
  If \`NULL\`, no check is performed.

- query_cell_type_col:

  The column name in the `colData` of `query_data` that identifies the
  cell types. If \`NULL\`, no check is performed.

- ref_cell_type_col:

  The column name in the `colData` of `reference_data` that identifies
  the cell types. If \`NULL\`, no check is performed.

- unique_cell_type:

  If \`TRUE\`, there should only be one cell type in the provided
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  objects. Default is \`FALSE\`.

- plot_function:

  A logical value indicating whether the function is being called to
  generate a plot. Default is \`FALSE\`.

- cell_names_query:

  A character vector of cell names in query data to be analyzed. If
  \`NULL\`, no check is performed.

- cell_names_ref:

  A character vector of cell names in reference data to be analyzed. If
  \`NULL\`, no check is performed.

- pc_subset_query:

  A numeric vector specifying the principal components to be used for
  the query data. If \`NULL\`, no check is performed.

- pc_subset_ref:

  A numeric vector specifying the principal components to be used for
  the reference data. If \`NULL\`, no check is performed.

- common_rotation_genes:

  If TRUE, check the rotation matrices of the reference and query data
  and ensure they have the same genes. Default is FALSE.

- assay_name:

  Name of the assay on which to perform computations. If \`NULL\`, no
  check is performed.

- max_cells_ref:

  Maximum number of reference cells to retain. If \`NULL\`, no check is
  performed.

- max_cells_query:

  Maximum number of query cells to retain. If \`NULL\`, no check is
  performed.

## Value

None.

## Details

The function performs a series of checks to ensure that:

- \`query_data\` and \`reference_data\` are
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  objects.

- \`query_cell_type_col\` and \`ref_cell_type_col\` exist in the column
  data of their respective
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  objects.

- If \`unique_cell_type\` is \`TRUE\`, there should only be one cell
  type in the
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  objects.

- \`cell_names_query\` are valid cell names in the provided query
  dataset.

- \`cell_names_ref\` are valid cell names in the provided reference
  dataset.

- The PCA subsets specified by \`pc_subset_query\` and \`pc_subset_ref\`
  are valid.

- \`max_cells_ref\` and \`max_cells_query\` are positive integers when
  not NULL.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
