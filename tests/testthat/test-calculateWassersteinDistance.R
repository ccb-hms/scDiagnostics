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
    expect_named(result, c("ref_ref_dist", "ref_query_dist", "probability_superiority", "cell_types"))
    expect_s3_class(result, "calculateWassersteinDistanceObject")

    # Check components of the result
    expect_type(result$ref_ref_dist, "list")
    expect_type(result$ref_query_dist, "list")
    expect_type(result$probability_superiority, "double")
    expect_type(result$cell_types, "character")
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

    expect_equal(length(result$cell_types), length(unique(ref_data_subset$expert_annotation)))
    expect_equal(length(result$probability_superiority), length(result$cell_types))
})

# Test to ensure that the function can run with default parameters
test_that("calculateWassersteinDistance works with default parameters", {
    result <- calculateWassersteinDistance(
        query_data = query_data_subset,
        reference_data = ref_data_subset,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    expect_type(result$ref_ref_dist, "list")
    expect_type(result$ref_query_dist, "list")
    expect_type(result$probability_superiority, "double")
})

# Test parameter validation
test_that("calculateWassersteinDistance validates parameters correctly", {
    # Test invalid n_resamples
    expect_error(calculateWassersteinDistance(
        query_data = query_data_subset,
        reference_data = ref_data_subset,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        n_resamples = -10
    ), "n_resamples.*should be an integer, greater than zero")

    expect_error(calculateWassersteinDistance(
        query_data = query_data_subset,
        reference_data = ref_data_subset,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        n_resamples = "invalid"
    ), "n_resamples.*should be numeric")
})

# Test with specific cell types
test_that("calculateWassersteinDistance works with specific cell types", {
    result <- calculateWassersteinDistance(
        query_data = query_data_subset,
        reference_data = ref_data_subset,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = "CD4",
        pc_subset = 1:5,
        n_resamples = 50
    )

    expect_equal(result$cell_types, "CD4")
    expect_equal(length(result$ref_ref_dist), 1)
    expect_equal(length(result$ref_query_dist), 1)
    expect_equal(length(result$probability_superiority), 1)
})
