# Calculate Cramer Test P-Values for Two-Sample Comparison of Multivariate ECDFs

This function performs the Cramer test for comparing multivariate
empirical cumulative distribution functions (ECDFs) between two samples.

## Usage

``` r
calculateCramerPValue(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  cell_types = NULL,
  pc_subset = 1:5,
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

A named vector of p-values from the Cramer test for each cell type.

## Details

The function performs the following steps:

1.  Projects the data into the PCA space.

2.  Subsets the data to the specified cell types and principal
    components.

3.  Performs the Cramer test for each cell type using the `cramer.test`
    function in the `cramer` package.

## References

Baringhaus, L., & Franz, C. (2004). "On a new multivariate two-sample
test". Journal of Multivariate Analysis, 88(1), 190-206.

## Examples

``` r
# Load data
data("reference_data")
data("query_data")

# Plot the PC data (with query data)
cramer_test <- calculateCramerPValue(reference_data = reference_data,
                                     query_data = query_data,
                                     ref_cell_type_col = "expert_annotation",
                                     query_cell_type_col = "SingleR_annotation",
                                     cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                                     pc_subset = 1:5)
cramer_test
#> B_and_plasma          CD4          CD8      Myeloid 
#>   0.11388611   0.08891109   0.00000000   0.92007992 
```
