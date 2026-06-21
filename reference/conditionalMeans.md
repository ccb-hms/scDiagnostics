# Compute Conditional Means for Cell Types

This function computes conditional means for each cell type in the
reference data. It can compute either a single conditional mean per cell
type or multiple conditional means, depending on the specified settings.
Principal component analysis (PCA) is used for dimensionality reduction
before clustering when computing multiple conditional means.

## Usage

``` r
conditionalMeans(
  reference_data,
  ref_cell_type_col,
  cell_types,
  multiple_cond_means = FALSE,
  assay_name = "logcounts",
  cumulative_variance_threshold = 0.7,
  n_neighbor = 1
)
```

## Arguments

- reference_data:

  A `SingleCellExperiment` object containing the reference data, where
  rows represent genes and columns represent cells.

- ref_cell_type_col:

  A character string specifying the column in `colData(reference_data)`
  that contains the cell type labels.

- cell_types:

  A character vector of cell types for which to compute conditional
  means.

- multiple_cond_means:

  A logical value indicating whether to compute multiple conditional
  means per cell type. Defaults to `FALSE`.

- assay_name:

  A character string specifying the name of the assay to use for the
  computation. Defaults to `"logcounts"`.

- cumulative_variance_threshold:

  A numeric value between 0 and 1 specifying the variance threshold for
  PCA when computing multiple conditional means. Defaults to `0.7`.

- n_neighbor:

  An integer specifying the number of nearest neighbors for clustering
  when computing multiple conditional means. Defaults to `1`.

## Value

A numeric matrix with the conditional means for each cell type. If
`multiple_cond_means = TRUE`, the matrix will contain multiple rows for
each cell type, representing the different conditional means computed
via clustering.

## Details

The function offers two modes of operation: - \*\*Single conditional
mean per cell type\*\*: For each cell type, it computes the mean
expression across all observations. - \*\*Multiple conditional means per
cell type\*\*: For each cell type, the function performs PCA to reduce
dimensionality, followed by clustering to compute multiple conditional
means.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
