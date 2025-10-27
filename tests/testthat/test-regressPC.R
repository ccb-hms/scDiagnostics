# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("regressPC works correctly with only query data", {
    # Perform regression with only query data
    regress_res <- regressPC(
        query_data = query_data,
        query_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    )

    # Check if the output is a list
    expect_true(is.list(regress_res))

    # Check if the list contains the correct elements
    expect_true("regression_summaries" %in% names(regress_res))
    expect_true("r_squared" %in% names(regress_res))
    expect_true("var_contributions" %in% names(regress_res))
    expect_true("total_variance_explained" %in% names(regress_res))
    expect_true("query_pca_var" %in% names(regress_res))
    expect_equal(regress_res$indep_var, "cell_type")

    # Check if R-squared values are numeric and of correct length
    expect_true(is.numeric(regress_res$r_squared))
    expect_equal(length(regress_res$r_squared), 5)

    # Check if total variance explained is a numeric value
    expect_true(is.numeric(regress_res$total_variance_explained))

    # Check if object has correct class
    expect_true("regressPCObject" %in% class(regress_res))
})

test_that("regressPC works correctly with reference and query data", {
    # Perform regression with reference and query data
    regress_res <- regressPC(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    )

    # Check if the output is a list
    expect_true(is.list(regress_res))

    # Check if the list contains the correct elements
    expect_true("regression_summaries" %in% names(regress_res))
    expect_true("r_squared" %in% names(regress_res))
    expect_true("reference_pca_var" %in% names(regress_res))
    expect_equal(regress_res$indep_var, "cell_type_dataset_interaction")

    # Check if R-squared values are numeric and of correct length
    expect_true(is.numeric(regress_res$r_squared))
    expect_equal(length(regress_res$r_squared), 5)

    # Check if object has correct class
    expect_true("regressPCObject" %in% class(regress_res))
})

test_that("regressPC works with batch information", {
    # Add batch information to query data for testing
    query_data$batch <- sample(c("batch1", "batch2"), ncol(query_data), replace = TRUE)

    # Perform regression with batch information
    regress_res <- regressPC(
        query_data = query_data,
        query_cell_type_col = "expert_annotation",
        query_batch_col = "batch",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    )

    # Check if the output is a list
    expect_true(is.list(regress_res))
    expect_equal(regress_res$indep_var, "cell_type_batch_interaction")

    # Check if object has correct class
    expect_true("regressPCObject" %in% class(regress_res))
})

test_that("regressPC handles incorrect parameters", {
    # Test for invalid query data column
    expect_error(regressPC(
        query_data = query_data,
        query_cell_type_col = "invalid_column",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    ))

    # Test for invalid reference data column
    expect_error(regressPC(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "invalid_column",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    ))

    # Test for invalid batch column
    expect_error(regressPC(
        query_data = query_data,
        query_cell_type_col = "expert_annotation",
        query_batch_col = "invalid_batch_column",
        pc_subset = 1:5
    ))
})

test_that("regressPC works with all cell types and default PC subset", {
    # Perform regression with all cell types and default PC subset
    regress_res <- regressPC(
        query_data = query_data,
        query_cell_type_col = "expert_annotation"
    )

    # Check if the output is a list
    expect_true(is.list(regress_res))

    # Check if the list contains the correct elements
    expect_true("regression_summaries" %in% names(regress_res))
    expect_true("r_squared" %in% names(regress_res))
    expect_true("var_contributions" %in% names(regress_res))
    expect_true("total_variance_explained" %in% names(regress_res))

    # Check if R-squared values are numeric and of correct length
    expect_true(is.numeric(regress_res$r_squared))
    expect_equal(length(regress_res$r_squared), 10)

    # Check if total variance explained is a numeric value
    expect_true(is.numeric(regress_res$total_variance_explained))

    # Check if object has correct class
    expect_true("regressPCObject" %in% class(regress_res))
})

