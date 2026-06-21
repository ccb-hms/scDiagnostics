# Downsample SingleCellExperiment Objects

This internal function downsamples
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
objects while preserving reducedDims coordinate information (PCA, UMAP,
t-SNE, etc.). Optionally, it can also subset by cell types before
downsampling.

## Usage

``` r
downsampleSCE(
  sce_object,
  cell_type_col = NULL,
  cell_types = NULL,
  max_cells = 2500,
  seed = NULL
)
```

## Arguments

- sce_object:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object to potentially downsample. May contain PCA, UMAP, TSNE, or
  other reducedDims.

- cell_type_col:

  The column name in colData that contains cell type information.
  Required if cell_types is not NULL. Default is NULL.

- cell_types:

  A character vector specifying which cell types to retain. If NULL, no
  cell type filtering is performed. Default is NULL.

- max_cells:

  Maximum number of cells to retain. If the object has fewer cells, it
  is returned unchanged. If NULL, no downsampling is performed (all
  cells are kept). Default is 2500.

- seed:

  Random seed for reproducible downsampling. If NULL, no seed is set.
  Default is NULL.

## Value

A
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object with at most max_cells cells and optionally filtered by cell
types. ReducedDims coordinates are preserved through standard
subsetting. If max_cells is NULL, the original object is returned
unchanged.

## Details

The function can perform cell type filtering followed by random
downsampling without replacement when the number of cells exceeds the
specified threshold. All reducedDims coordinates are preserved through
standard sce_object subsetting operations. This function does not
preserve PCA rotation matrices or other model-specific attributes.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
