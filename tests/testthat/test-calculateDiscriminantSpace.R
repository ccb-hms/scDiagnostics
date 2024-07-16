# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Unit tests for calculateDiscriminantSpace function
test_that("calculateDiscriminantSpace function works as expected", {
    
    # Test Case 1: Basic Functionality
    output <- calculateDiscriminantSpace(reference_data = reference_data,
                                         query_data = NULL,
                                         ref_cell_type_col = "expert_annotation",
                                         query_cell_type_col = "SingleR_annotation",
                                         n_tree = 500,
                                         n_top = 20,
                                         eigen_threshold = 0.1,
                                         calculate_metrics = FALSE,
                                         alpha = 0.01)
    
    expect_type(output, "list")
    expect_true(length(output) > 0)
    
    # Test Case 2: Projecting Reference Data Only
    output_ref_only <- calculateDiscriminantSpace(reference_data = reference_data,
                                                  query_data = NULL,
                                                  ref_cell_type_col = "expert_annotation",
                                                  query_cell_type_col = "SingleR_annotation",
                                                  n_tree = 500,
                                                  n_top = 20,
                                                  eigen_threshold = 0.1,
                                                  calculate_metrics = FALSE,
                                                  alpha = 0.01)
    
    expect_type(output_ref_only, "list")
    expect_true(length(output_ref_only) > 0)
    
    # Test Case 3: Projecting Both Reference and Query Data
    output_with_query <- calculateDiscriminantSpace(reference_data = reference_data,
                                                    query_data = query_data,
                                                    ref_cell_type_col = "expert_annotation",
                                                    query_cell_type_col = "SingleR_annotation",
                                                    n_tree = 500,
                                                    n_top = 20,
                                                    eigen_threshold = 0.1,
                                                    calculate_metrics = TRUE,
                                                    alpha = 0.01)
    
    expect_type(output_with_query, "list")
    expect_true(length(output_with_query) > 0)
    expect_true(all(c("discriminant_eigenvalues", "discriminant_eigenvectors", "ref_proj") %in% names(output_with_query[[1]])))

    # Test Case 4: Edge Cases
    # Test with different values of n_tree, eigen_threshold, alpha, etc.
    output_edge <- calculateDiscriminantSpace(reference_data = reference_data,
                                              query_data = query_data,
                                              ref_cell_type_col = "expert_annotation",
                                              query_cell_type_col = "SingleR_annotation",
                                              n_tree = 100,
                                              n_top = 10,
                                              eigen_threshold = 0.01,
                                              calculate_metrics = TRUE,
                                              alpha = 0.05)
    
    expect_type(output_edge, "list")
    expect_true(length(output_edge) > 0)
})
