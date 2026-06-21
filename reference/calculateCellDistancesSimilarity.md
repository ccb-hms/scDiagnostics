# Function to Calculate Bhattacharyya Coefficients and Hellinger Distances

This function computes Bhattacharyya coefficients and Hellinger
distances to quantify the similarity of density distributions between
query cells and reference data for each cell type.

## Usage

``` r
calculateCellDistancesSimilarity(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  cell_types = NULL,
  cell_names_query,
  pc_subset = 1:5,
  assay_name = "logcounts",
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

- cell_names_query:

  A character vector specifying the names of the query cells for which
  to compute distance measures.

- pc_subset:

  A numeric vector specifying which principal components to include in
  the plot. Default is 1:5.

- assay_name:

  Name of the assay on which to perform computations. Default is
  "logcounts".

- max_cells_ref:

  Maximum number of reference cells to retain after cell type filtering.
  If NULL, no downsampling of reference cells is performed. Default is
  5000.

## Value

A list containing distance data for each cell type. Each entry in the
list contains:

- ref_distances:

  A vector of all pairwise distances within the reference subset for the
  cell type.

- query_to_ref_distances:

  A matrix of distances from each query cell to all reference cells for
  the cell type.

## Details

This function first computes distance data using the
`calculateCellDistances` function, which calculates pairwise distances
between cells within the reference data and between query cells and
reference cells in the PCA space. Bhattacharyya coefficients and
Hellinger distances are calculated to quantify the similarity of density
distributions between query cells and reference data for each cell type.
Bhattacharyya coefficient measures the similarity of two probability
distributions, while Hellinger distance measures the distance between
two probability distributions.

Bhattacharyya coefficients range between 0 and 1. A value closer to 1
indicates higher similarity between distributions, while a value closer
to 0 indicates lower similarity

Hellinger distances range between 0 and 1. A value closer to 0 indicates
higher similarity between distributions, while a value closer to 1
indicates lower similarity.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>

## Examples

``` r
# Load data
data("reference_data")
data("query_data")

# Plot the PC data
distance_data <- calculateCellDistances(query_data = query_data,
                                        reference_data = reference_data,
                                        query_cell_type_col = "SingleR_annotation",
                                        ref_cell_type_col = "expert_annotation",
                                        pc_subset = 1:10)

# Identify outliers for CD4
cd4_anomalies <- detectAnomaly(reference_data = reference_data,
                               query_data = query_data,
                               query_cell_type_col = "SingleR_annotation",
                               ref_cell_type_col = "expert_annotation",
                               pc_subset = 1:10,
                               n_tree = 500,
                               anomaly_threshold = 0.5)
cd4_top6_anomalies <- names(sort(cd4_anomalies$CD4$query_anomaly_scores, decreasing = TRUE)[1:6])

# Get overlap measures
overlap_measures <- calculateCellDistancesSimilarity(query_data = query_data,
                                                     reference_data = reference_data,
                                                     cell_names_query = cd4_top6_anomalies,
                                                     query_cell_type_col = "SingleR_annotation",
                                                     ref_cell_type_col = "expert_annotation",
                                                     pc_subset = 1:10)
overlap_measures
#> $bhattacharyya_coef
#>                       Cell B_and_plasma        CD4        CD8   Myeloid
#> 1 Query_AACCGCGCAAACCCAT-1    0.6063560 0.14963173 0.15973512 0.2146414
#> 2 Query_GTTACAGCAGCTGTGC-1    0.5716079 0.05978566 0.07361651 0.1379982
#> 3 Query_TGAGCATTCCAAACTG-1    0.5797450 0.21854137 0.24409183 0.2559222
#> 4 Query_TGCTACCTCGGTTCGG-1    0.3480598 0.76503890 0.96527336 0.2815120
#> 5 Query_TCAGCTCCATACTCTT-1    0.2903559 0.72368485 0.98681446 0.2965596
#> 6 Query_CTCGAAATCAGTCCCT-1    0.2975623 0.75713220 0.93154956 0.2300326
#> 
#> $hellinger_dist
#>                       Cell B_and_plasma       CD4       CD8   Myeloid
#> 1 Query_AACCGCGCAAACCCAT-1    0.6274106 0.9221541 0.9166596 0.8862046
#> 2 Query_GTTACAGCAGCTGTGC-1    0.6545167 0.9696465 0.9624882 0.9284405
#> 3 Query_TGAGCATTCCAAACTG-1    0.6482708 0.8840015 0.8694298 0.8625994
#> 4 Query_TGCTACCTCGGTTCGG-1    0.8074282 0.4847279 0.1863509 0.8476367
#> 5 Query_TCAGCTCCATACTCTT-1    0.8424037 0.5256569 0.1148283 0.8387135
#> 6 Query_CTCGAAATCAGTCCCT-1    0.8381156 0.4928162 0.2616304 0.8774779
#> 
```
