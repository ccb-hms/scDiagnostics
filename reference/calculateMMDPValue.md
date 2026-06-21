# Calculate Maximum Mean Discrepancy P-Values for Two-Sample Comparison

This function performs the Maximum Mean Discrepancy (MMD) test for
comparing distributions between two samples in PCA space using a custom
implementation with permutation testing for better sensitivity.

## Usage

``` r
calculateMMDPValue(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  cell_types = NULL,
  pc_subset = 1:5,
  n_permutation = 100,
  kernel_type = "gaussian",
  sigma = NULL,
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

  The column name in the `colData` of `query_data` that identifies the
  cell types.

- ref_cell_type_col:

  The column name in the `colData` of `reference_data` that identifies
  the cell types.

- cell_types:

  A character vector specifying the cell types to include in the plot.
  If NULL, all cell types are included.

- pc_subset:

  A numeric vector specifying which principal components to include in
  the plot. Default is PC1 to PC5.

- n_permutation:

  Number of permutations for p-value calculation. Default is 100.

- kernel_type:

  Type of kernel to use. Options are "gaussian" (default) or "linear".

- sigma:

  Bandwidth parameter for Gaussian kernel. If NULL, uses median
  heuristic.

- assay_name:

  Name of the assay on which to perform computations. Default is
  "logcounts".

- max_cells_query:

  Maximum number of query cells to retain after cell type filtering. If
  NULL, no downsampling of query cells is performed. Default is 5000.

- max_cells_ref:

  Maximum number of reference cells to retain after cell type filtering.
  If NULL, no downsampling of reference cells is performed. Default is
  5000.

## Value

A named vector of p-values from the MMD test for each cell type.

## Details

The function performs the following steps:

1.  Projects the data into the PCA space.

2.  Subsets the data to the specified cell types and principal
    components.

3.  Performs a custom MMD test with permutation-based p-values for each
    cell type.

## References

Gretton, A., Borgwardt, K. M., Rasch, M. J., Schölkopf, B., & Smola, A.
(2012). "A kernel two-sample test". Journal of Machine Learning
Research, 13(1), 723-773.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>

## Examples

``` r
# Load data
data("reference_data")
data("query_data")

# Calculate MMD p-values (with query data)
mmd_test <- calculateMMDPValue(reference_data = reference_data,
                              query_data = query_data,
                              ref_cell_type_col = "expert_annotation",
                              query_cell_type_col = "SingleR_annotation",
                              cell_types = c("CD4", "CD8"),
                              pc_subset = 1:5,
                              n_permutation = 30)
mmd_test
#>        CD4        CD8 
#> 0.03225806 0.03225806 
```
