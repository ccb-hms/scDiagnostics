# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("calculateSIRSpace returns correct structure", {
    # Run the calculateSIRSpace function
    result <- calculateSIRSpace(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    # Check that result is a list and contains expected components
    expect_type(result, "list")
    expect_named(result, c("cond_means", "rotation_mat", "sir_projections", "percent_var"))

    # Check that the class includes "calculateSIRSpace"
    expect_s3_class(result, "calculateSIRSpaceObject")
})

test_that("calculateSIRSpace uses the correct cumulative variance threshold", {
    # Run the function with a specific cumulative variance threshold
    result <- calculateSIRSpace(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cumulative_variance_threshold = 0.9
    )

    # Ensure the cumulative variance threshold is respected
    expect_true(sum(result$percent_var) >= 0.9)
})

test_that("calculateSIRSpace works with specific cell types", {
    # Test with a specified subset of cell types
    cell_types <- c("CD4", "CD8")

    result <- calculateSIRSpace(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = cell_types
    )

    # Check that the conditional means are computed only for the specified cell types
    expect_true(all(rownames(result$cond_means) %in% cell_types))
})

test_that("calculateSIRSpace handles multiple conditional means", {
    # Test with multiple_cond_means set to FALSE
    result <- calculateSIRSpace(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        multiple_cond_means = FALSE
    )

    # Ensure that the conditional means are calculated based on single condition
    expect_true(is.matrix(result$cond_means)) # Expect conditional means to be a matrix
})

test_that("calculateSIRSpace handles n_neighbor parameter", {
    # Test with n_neighbor = 3
    result <- calculateSIRSpace(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        n_neighbor = 3
    )

    # No direct way to test n_neighbor impact, but ensure result structure is correct
    expect_type(result$sir_projections, "list")
})

test_that("calculateSIRSpace throws error for missing arguments", {
    # Test if function throws an error for missing required arguments
    expect_error(calculateSIRSpace(query_data = query_data), "argument \"reference_data\" is missing")
    expect_error(calculateSIRSpace(reference_data = reference_data), "argument \"query_data\" is missing")
    expect_error(calculateSIRSpace(query_data = query_data, reference_data = reference_data), "argument \"query_cell_type_col\" is missing")
})
