# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Unit tests for calculateDiscriminantSpace function
test_that("calculateDiscriminantSpace function works as expected", {

    # Test Case 1: Basic Functionality - Reference data only
    output <- calculateDiscriminantSpace(reference_data = reference_data,
                                         query_data = NULL,
                                         ref_cell_type_col = "expert_annotation",
                                         n_tree = 500,
                                         n_top = 20,
                                         eigen_threshold = 0.1,
                                         calculate_metrics = FALSE,
                                         alpha = 0.01)

    expect_type(output, "list")
    expect_s3_class(output, "calculateDiscriminantSpaceObject")
    expect_true(all(c("discriminant_eigenvalues", "discriminant_eigenvectors", "ref_proj") %in% names(output)))

    # Test Case 2: With query data but no metrics
    output_with_query <- calculateDiscriminantSpace(reference_data = reference_data,
                                                    query_data = query_data,
                                                    ref_cell_type_col = "expert_annotation",
                                                    query_cell_type_col = "SingleR_annotation",
                                                    n_tree = 500,
                                                    n_top = 20,
                                                    eigen_threshold = 0.1,
                                                    calculate_metrics = FALSE,
                                                    alpha = 0.01)

    expect_type(output_with_query, "list")
    expect_s3_class(output_with_query, "calculateDiscriminantSpaceObject")
    expect_true(all(c("discriminant_eigenvalues", "discriminant_eigenvectors", "ref_proj", "query_proj") %in% names(output_with_query)))

    # Test Case 3: With query data and metrics
    output_with_metrics <- calculateDiscriminantSpace(reference_data = reference_data,
                                                      query_data = query_data,
                                                      ref_cell_type_col = "expert_annotation",
                                                      query_cell_type_col = "SingleR_annotation",
                                                      n_tree = 100,
                                                      n_top = 10,
                                                      eigen_threshold = 0.1,
                                                      calculate_metrics = TRUE,
                                                      alpha = 0.01)

    expect_type(output_with_metrics, "list")
    expect_s3_class(output_with_metrics, "calculateDiscriminantSpaceObject")
    expect_true(all(c("discriminant_eigenvalues", "discriminant_eigenvectors", "ref_proj", "query_proj",
                      "query_mahalanobis_dist", "mahalanobis_crit", "query_cosine_similarity") %in% names(output_with_metrics)))
})

test_that("calculateDiscriminantSpace parameter validation works", {
    # Test invalid n_tree
    expect_error(calculateDiscriminantSpace(reference_data = reference_data,
                                            ref_cell_type_col = "expert_annotation",
                                            n_tree = -10),
                 "n_tree.*must be a positive integer")

    # Test invalid n_top
    expect_error(calculateDiscriminantSpace(reference_data = reference_data,
                                            ref_cell_type_col = "expert_annotation",
                                            n_top = 0),
                 "n_top.*must be a positive integer")

    # Test invalid eigen_threshold
    expect_error(calculateDiscriminantSpace(reference_data = reference_data,
                                            ref_cell_type_col = "expert_annotation",
                                            eigen_threshold = -0.1),
                 "eigen_threshold.*must be a positive number")

    # Test invalid alpha
    expect_error(calculateDiscriminantSpace(reference_data = reference_data,
                                            ref_cell_type_col = "expert_annotation",
                                            alpha = 1.5),
                 "alpha.*must be a positive number.*less than 1")
})

test_that("calculateDiscriminantSpace works with specific cell types", {
    output <- calculateDiscriminantSpace(reference_data = reference_data,
                                         query_data = query_data,
                                         ref_cell_type_col = "expert_annotation",
                                         query_cell_type_col = "SingleR_annotation",
                                         cell_types = c("CD4", "CD8"),
                                         n_tree = 100,
                                         n_top = 10,
                                         eigen_threshold = 0.1)

    expect_type(output, "list")
    expect_s3_class(output, "calculateDiscriminantSpaceObject")
    expect_true("ref_proj" %in% names(output))
    expect_true("query_proj" %in% names(output))
})
