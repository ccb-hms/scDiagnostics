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

# Unit tests for comparePCA
test_that("comparePCA returns expected output structure", {
    # Run the comparePCA function
    result <- comparePCA(query_data = query_data_subset,
                         reference_data = reference_data_subset,
                         query_cell_type_col = "expert_annotation",
                         ref_cell_type_col = "expert_annotation",
                         pc_subset = 1:5,
                         n_top_vars = 50,
                         metric = "cosine",
                         correlation_method = "spearman")

    # Check the class of the output
    expect_s3_class(result, "comparePCAObject")

    # Check the dimensions of the similarity matrix
    expect_equal(dim(result$similarity_matrix), c(5, 5))

    # Check that the matrix contains numeric values
    expect_type(result$similarity_matrix, "double")
})

test_that("comparePCA handles invalid n_top_vars argument", {
    # Check if function stops with invalid n_top_vars
    expect_error(comparePCA(query_data = query_data_subset,
                            reference_data = reference_data_subset,
                            query_cell_type_col = "expert_annotation",
                            ref_cell_type_col = "expert_annotation",
                            pc_subset = 1:5,
                            n_top_vars = -10,
                            metric = "cosine",
                            correlation_method = "spearman"))
})

test_that("comparePCA handles invalid metric argument", {
    # Check if function stops with invalid metric
    expect_error(comparePCA(query_data = query_data_subset,
                            reference_data = reference_data_subset,
                            query_cell_type_col = "expert_annotation",
                            ref_cell_type_col = "expert_annotation",
                            pc_subset = 1:5,
                            n_top_vars = 50,
                            metric = "invalid_metric",
                            correlation_method = "spearman"))
})

test_that("comparePCA handles invalid correlation_method argument", {
    # Check if function stops with invalid correlation method
    expect_error(comparePCA(query_data = query_data_subset,
                            reference_data = reference_data_subset,
                            query_cell_type_col = "expert_annotation",
                            ref_cell_type_col = "expert_annotation",
                            pc_subset = 1:5,
                            n_top_vars = 50,
                            metric = "correlation",
                            correlation_method = "invalid_method"))
})

test_that("comparePCA works with correlation metric", {
    # Run the comparePCA function with correlation metric
    result <- comparePCA(query_data = query_data_subset,
                         reference_data = reference_data_subset,
                         query_cell_type_col = "expert_annotation",
                         ref_cell_type_col = "expert_annotation",
                         pc_subset = 1:5,
                         n_top_vars = 50,
                         metric = "correlation",
                         correlation_method = "spearman")

    # Check the class of the output
    expect_s3_class(result, "comparePCAObject")

    # Check the dimensions of the similarity matrix
    expect_equal(dim(result$similarity_matrix), c(5, 5))

    # Check that the matrix contains numeric values
    expect_type(result$similarity_matrix, "double")
})

test_that("comparePCA works with different correlation methods", {
    # Run the comparePCA function with pearson correlation
    result <- comparePCA(query_data = query_data_subset,
                         reference_data = reference_data_subset,
                         query_cell_type_col = "expert_annotation",
                         ref_cell_type_col = "expert_annotation",
                         pc_subset = 1:5,
                         n_top_vars = 50,
                         metric = "correlation",
                         correlation_method = "pearson")

    # Check the class of the output
    expect_s3_class(result, "comparePCAObject")

    # Check the dimensions of the similarity matrix
    expect_equal(dim(result$similarity_matrix), c(5, 5))

    # Check that the matrix contains numeric values
    expect_type(result$similarity_matrix, "double")
})

test_that("comparePCA produces expected results", {
    # Run the comparePCA function with correlation similarity
    result <- comparePCA(query_data = query_data_subset,
                         reference_data = reference_data_subset,
                         query_cell_type_col = "expert_annotation",
                         ref_cell_type_col = "expert_annotation",
                         pc_subset = 1:5,
                         n_top_vars = 50,
                         metric = "correlation",
                         correlation_method = "spearman")

    # Check if the correlation values are in the range [-1, 1]
    expect_true(all(result$similarity_matrix >= -1 & result$similarity_matrix <= 1))
})
