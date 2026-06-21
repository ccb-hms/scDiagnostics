# Compare Gene Importance Across Datasets Using Random Forest

This function identifies and compares the most important genes for
differentiating cell types between a query dataset and a reference
dataset using Random Forest.

## Usage

``` r
calculateVarImpOverlap(
  reference_data,
  query_data = NULL,
  ref_cell_type_col,
  query_cell_type_col = NULL,
  cell_types = NULL,
  n_tree = 500,
  n_top = 50,
  assay_name = "logcounts",
  max_cells_ref = 5000,
  max_cells_query = 5000
)
```

## Arguments

- reference_data:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing numeric expression matrix for the reference cells.

- query_data:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing numeric expression matrix for the query cells. If
  NULL, then the variable importance scores are only computed for the
  reference data. Default is NULL.

- ref_cell_type_col:

  A character string specifying the column name in the reference dataset
  containing cell type annotations.

- query_cell_type_col:

  A character string specifying the column name in the query dataset
  containing cell type annotations.

- cell_types:

  A character vector specifying the cell types to include in the plot.
  If NULL, all cell types are included.

- n_tree:

  An integer specifying the number of trees to grow in the Random
  Forest. Default is 500.

- n_top:

  An integer specifying the number of top genes to consider when
  comparing variable importance scores. Default is 50.

- assay_name:

  Name of the assay on which to perform computations. Defaults to
  `"logcounts"`.

- max_cells_ref:

  Maximum number of reference cells to retain after cell type filtering.
  If NULL, no downsampling of reference cells is performed. Default is
  5000.

- max_cells_query:

  Maximum number of query cells to retain after cell type filtering. If
  NULL, no downsampling of query cells is performed. Default is 5000.

## Value

A list containing three elements:

- var_imp_ref:

  A list of data frames containing variable importance scores for each
  combination of cell types in the reference dataset.

- var_imp_query:

  A list of data frames containing variable importance scores for each
  combination of cell types in the query dataset.

- var_imp_comparison:

  A named vector indicating the proportion of top genes that overlap
  between the reference and query datasets for each combination of cell
  types.

## Details

This function uses the Random Forest algorithm to calculate the
importance of genes in differentiating between cell types within both a
reference dataset and a query dataset. The function then compares the
top genes identified in both datasets to determine the overlap in their
importance scores.

## References

Breiman, L. (2001). "Random forests". \*Machine Learning\*, 45(1), 5-32.
doi:10.1023/A:1010933404324.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>

## Examples

``` r
# Load data
data("reference_data")
data("query_data")

# Compute important variables for all pairwise cell comparisons
rf_output <- calculateVarImpOverlap(reference_data = reference_data,
                                    query_data = query_data,
                                    query_cell_type_col = "SingleR_annotation",
                                    ref_cell_type_col = "expert_annotation",
                                    n_tree = 500,
                                    n_top = 50)

# Comparison table
rf_output$var_imp_comparison
#>              CD4-CD8     CD4-B_and_plasma          CD4-Myeloid 
#>                 0.76                 0.80                 0.82 
#>     CD8-B_and_plasma          CD8-Myeloid B_and_plasma-Myeloid 
#>                 0.76                 0.78                 0.80 
```
