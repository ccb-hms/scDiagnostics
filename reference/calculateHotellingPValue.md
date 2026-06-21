# Perform Hotelling's T-squared Test on PCA Scores for Single-cell RNA-seq Data

Computes Hotelling's T-squared test statistic and p-values for each
specified cell type based on PCA-projected data from query and reference
datasets.

## Usage

``` r
calculateHotellingPValue(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  cell_types = NULL,
  pc_subset = 1:5,
  n_permutation = 500,
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

  character. The column name in the `colData` of `query_data` that
  identifies the cell types.

- ref_cell_type_col:

  character. The column name in the `colData` of `reference_data` that
  identifies the cell types.

- cell_types:

  A character vector specifying the cell types to include in the plot.
  If NULL, all cell types are included.

- pc_subset:

  A numeric vector specifying which principal components to include in
  the plot. Default is PC1 to PC5.

- n_permutation:

  Number of permutations to perform for p-value calculation. Default is
  500.

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

A named numeric vector of p-values from Hotelling's T-squared test for
each cell type.

## Details

This function calculates Hotelling's T-squared statistic for comparing
multivariate means between reference and query datasets, projected onto
a subset of principal components (PCs). It performs a permutation test
to obtain p-values for each cell type specified.

## References

Hotelling, H. (1931). "The generalization of Student's ratio". \*Annals
of Mathematical Statistics\*. 2 (3): 360–378.
doi:10.1214/aoms/1177732979.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>

## Examples

``` r
# Load data
data("reference_data")
data("query_data")

# Get the p-values
p_values <- calculateHotellingPValue(query_data = query_data,
                                     reference_data = reference_data,
                                     query_cell_type_col = "SingleR_annotation",
                                     ref_cell_type_col = "expert_annotation",
                                     pc_subset = 1:10)
round(p_values, 5)
#> B_and_plasma          CD4          CD8      Myeloid 
#>        0.228        0.002        0.000        0.382 
```
