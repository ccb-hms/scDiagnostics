# Process PCA for SingleCellExperiment Objects

This function ensures that a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object has valid PCA computed using highly variable genes when needed.
It only performs downsampling when PCA computation is required,
preserving existing valid PCA computations without modification.

## Usage

``` r
processPCA(
  sce_object,
  assay_name = "logcounts",
  n_hvgs = 2000,
  max_cells = NULL
)
```

## Arguments

- sce_object:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object to process.

- assay_name:

  Name of the assay to use for HVG selection and PCA computation. Should
  contain log-normalized expression values. Default is "logcounts".

- n_hvgs:

  Number of highly variable genes to select for PCA computation. Default
  is 2000.

- max_cells:

  Maximum number of cells to retain if downsampling is needed for PCA
  computation. If NULL, no downsampling is performed. Default is NULL.

## Value

A
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object with valid PCA in the reducedDims slot, including rotation matrix
and percentVar attributes. Will have original cell count if PCA was
valid, or at most max_cells if PCA was computed.

## Details

The function performs the following operations:

- Checks if PCA exists and is valid in the provided
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- Validates PCA integrity including rotation matrix, percentVar, gene
  consistency, and dimensions

- If PCA is valid, returns the object unchanged (no downsampling)

- If PCA is missing or invalid and dataset is large, downsamples before
  computing PCA

- Computes PCA using highly variable genes when PCA is missing or
  invalid

- Utilizes scran for HVG selection and scater for PCA computation (soft
  dependencies)

The downsampling strategy uses random sampling without replacement and
only occurs when PCA computation is necessary. This preserves expensive
pre-computed PCA results while ensuring computational efficiency for new
PCA computations.

PCA validation includes checking for:

- Presence of PCA in reducedDims

- Existence of rotation matrix and percentVar attributes

- Gene consistency between rotation matrix and current assay

- Dimension consistency between PCA coordinates and cell count

## Note

This function requires the scran and scater packages for HVG selection
and PCA computation. These packages should be installed via
BiocManager::install(c("scran", "scater")).

Objects with existing valid PCA are returned unchanged to preserve
expensive pre-computations. Only datasets requiring PCA computation are
subject to downsampling.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>

## Examples

``` r
library(SingleCellExperiment)

# Load data
data("reference_data")
data("query_data")

# Example 1: Dataset without PCA (will compute PCA)
query_no_pca <- query_data
reducedDims(query_no_pca) <- list()  # Remove existing PCA

processed_query <- processPCA(sce_object = query_no_pca, n_hvgs = 500)
#> Data missing PCA - computing...
#> Computing PCA...
#> Warning: 'fitTrendVar' is deprecated.
#> Use 'scrapper::fitVarianceTrend' instead.
#> See help("Deprecated")
#> Warning: 'combineBlocks' is deprecated.
#> See help("Deprecated")
#> Warning: 'scran::getTopHVGs' is deprecated.
#> Use 'scrapper::chooseHighlyVariableGenes' instead.
#> See help("Deprecated")
#> Using 231 highly variable genes for PCA computation
"PCA" %in% reducedDimNames(processed_query)  # Should be TRUE
#> [1] TRUE
ncol(processed_query)  # Should be 503 (unchanged)
#> [1] 503

# Example 2: Dataset with existing valid PCA (will be preserved)
processed_existing <- processPCA(sce_object = query_data, n_hvgs = 500)
#> Data already has valid PCA - returning unchanged (503 cells)
ncol(processed_existing)  # Should be 503 (unchanged, no downsampling)
#> [1] 503

# Example 3: Large dataset requiring downsampling for PCA computation
ref_no_pca <- reference_data
reducedDims(ref_no_pca) <- list()  # Remove existing PCA

processed_large <- processPCA(sce_object = ref_no_pca,
                              n_hvgs = 800,
                              max_cells = 1000)
#> Data missing PCA - computing...
#> Downsampling data from 1500 to 1000 cells before PCA computation
#> Computing PCA...
#> Warning: 'fitTrendVar' is deprecated.
#> Use 'scrapper::fitVarianceTrend' instead.
#> See help("Deprecated")
#> Warning: 'combineBlocks' is deprecated.
#> See help("Deprecated")
#> Warning: 'scran::getTopHVGs' is deprecated.
#> Use 'scrapper::chooseHighlyVariableGenes' instead.
#> See help("Deprecated")
#> Using 235 highly variable genes for PCA computation
ncol(processed_large)  # Should be 1000 (downsampled for PCA computation)
#> [1] 1000

# Example 4: Large dataset with existing PCA (no downsampling)
processed_large_existing <- processPCA(sce_object = reference_data,
                                       n_hvgs = 800,
                                       max_cells = 1000)
#> Data already has valid PCA - returning unchanged (1500 cells)
ncol(processed_large_existing)  # Should be 1500 (preserved, no downsampling)
#> [1] 1500
```
