# Cell Type Selection and Validation for SingleCellExperiment Analysis

This function selects and validates cell types for functions that
analyze
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
objects. It determines which cell types to include based on availability
in datasets, applies filtering criteria, and optionally selects the top
cell types by cell count.

## Usage

``` r
selectCellTypes(
  query_data = NULL,
  reference_data = NULL,
  query_cell_type_col = NULL,
  ref_cell_type_col = NULL,
  cell_types = NULL,
  dual_only = FALSE,
  n_cell_types = NULL
)
```

## Arguments

- query_data:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing numeric expression matrix for the query cells. Can
  be `NULL` if only reference data is available.

- reference_data:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing numeric expression matrix for the reference cells.
  Can be `NULL` if only query data is available.

- query_cell_type_col:

  The column name in the `colData` of `query_data` that identifies the
  cell types. Should be `NULL` if `query_data` is `NULL`.

- ref_cell_type_col:

  The column name in the `colData` of `reference_data` that identifies
  the cell types. Should be `NULL` if `reference_data` is `NULL`.

- cell_types:

  A character vector specifying the cell types to validate. If `NULL`,
  cell types will be automatically selected based on dataset
  availability and `dual_only` setting.

- dual_only:

  A logical value indicating whether cell types must be present in both
  datasets. If `TRUE`, both `query_data` and `reference_data` must be
  provided, and only cell types present in both datasets will be
  considered. Default is `FALSE`.

- n_cell_types:

  An integer specifying the maximum number of cell types to select based
  on highest cell count. If `NULL`, all valid cell types are returned.
  Default is `NULL`.

## Value

A character vector of selected cell types that meet the specified
criteria.

## Details

The function performs the following selection and validation steps:

- Validates that at least one of `query_data` or `reference_data` is
  provided.

- When `dual_only` is TRUE, ensures both datasets are provided.

- Determines available cell types based on dataset availability and
  `dual_only` setting.

- If `cell_types` is NULL and both datasets are available, includes cell
  types based on `dual_only`.

- If `cell_types` is NULL and only one dataset is available, includes
  all cell types from that dataset.

- If `cell_types` is provided, filters to include only valid types.

- If `n_cell_types` is specified, selects the top cell types by total
  cell count.

- Returns the selected and validated cell types as character strings.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
