# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Extract CD4 cells
reference_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]

# Selecting highly variable genes (can be customized by the user)
ref_top_genes <- scran::getTopHVGs(reference_data_subset, n = 500)
query_top_genes <- scran::getTopHVGs(query_data_subset, n = 500)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_top_genes, query_top_genes)
reference_data_subset <- reference_data_subset[common_genes,]
query_data_subset <- query_data_subset[common_genes,]

# Run PCA on datasets separately
reference_data_subset <- scater::runPCA(reference_data_subset)
query_data_subset <- scater::runPCA(query_data_subset)

# Define unit tests
test_that("compareCCA function works correctly", {
    # Test correct output structure and class
    cca_comparison <- compareCCA(query_data = query_data_subset,
                                 reference_data = reference_data_subset,
                                 query_cell_type_col = "expert_annotation",
                                 ref_cell_type_col = "expert_annotation",
                                 pc_subset = 1:5)

    expect_type(cca_comparison, "list")
    expect_s3_class(cca_comparison, "compareCCAObject")
    expect_named(cca_comparison, c("coef_ref", "coef_query", "cosine_similarity", "correlations"))

    # Test with subset of principal components
    cca_comparison <- compareCCA(query_data = query_data_subset,
                                 reference_data = reference_data_subset,
                                 query_cell_type_col = "expert_annotation",
                                 ref_cell_type_col = "expert_annotation",
                                 pc_subset = 1:3)
    expect_equal(length(cca_comparison$cosine_similarity), 3)
    expect_equal(length(cca_comparison$correlations), 3)
})
