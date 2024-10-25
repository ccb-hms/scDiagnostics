# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Extract CD4 cells for subset analysis
ref_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]

# Selecting highly variable genes (can be customized by the user)
ref_top_genes <- scran::getTopHVGs(ref_data_subset, n = 500)
query_top_genes <- scran::getTopHVGs(query_data_subset, n = 500)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_top_genes, query_top_genes)
ref_data_subset <- ref_data_subset[common_genes, ]
query_data_subset <- query_data_subset[common_genes, ]

# Run PCA on reference data subset
ref_data_subset <- scater::runPCA(ref_data_subset)

# Test to ensure that the function returns a list with expected components
test_that("calculateWassersteinDistance returns correct object structure", {
    result <- calculateWassersteinDistance(
        query_data = query_data_subset,
        reference_data = ref_data_subset,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:5,
        n_resamples = 100
    )

    expect_true(is.list(result))
    expect_named(result, c("null_dist", "query_dist", "cell_type"))

    # Check components of the result
    expect_type(result$null_dist, "double")
    expect_type(result$query_dist, "double")
    expect_type(result$cell_type, "character")
})

# Test to ensure that calculateWassersteinDistance correctly identifies a specific number of unique cell types
test_that("calculateWassersteinDistance finds the right number of cell types", {
    result <- calculateWassersteinDistance(
        query_data = query_data_subset,
        reference_data = ref_data_subset,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:5,
        n_resamples = 100
    )

    expect_equal(length(result$cell_type), length(unique(ref_data_subset$expert_annotation)))
})

# Test to ensure that the function can run with default parameters
test_that("calculateWassersteinDistance works with default parameters", {
    result <- calculateWassersteinDistance(
        query_data = query_data_subset,
        reference_data = ref_data_subset,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    expect_type(result$null_dist, "double")
    expect_type(result$query_dist, "double")
})
