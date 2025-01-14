# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("detectAnomaly works correctly with default parameters", {
    # Store PCA anomaly data
    anomaly_output <- detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation"
    )

    # Check if the output is a list
    expect_true(is.list(anomaly_output))

    # Check if the output contains the expected elements
    expect_true("CD4" %in% names(anomaly_output))
    expect_true("Combined" %in% names(anomaly_output))

    # Check the structure of the elements for a specific cell type
    expect_true(all(c("reference_anomaly_scores", "reference_anomaly", "reference_mat_subset", "query_mat_subset", "query_anomaly_scores", "query_anomaly", "var_explained") %in% names(anomaly_output[["CD4"]])))
})

test_that("detectAnomaly works correctly without query data", {
    # Store PCA anomaly data without query data
    anomaly_output <- detectAnomaly(
        reference_data = reference_data,
        ref_cell_type_col = "expert_annotation"
    )

    # Check if the output is a list
    expect_true(is.list(anomaly_output))

    # Check if the output contains the expected elements
    expect_true("CD4" %in% names(anomaly_output))
    expect_true("Combined" %in% names(anomaly_output))

    # Check the structure of the elements for a specific cell type
    expect_true(all(c("reference_anomaly_scores", "reference_anomaly", "reference_mat_subset", "var_explained") %in% names(anomaly_output[["CD4"]])))
    expect_false("query_anomaly_scores" %in% names(anomaly_output[["CD4"]]))
})

test_that("detectAnomaly handles incorrect parameters", {
    expect_error(detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        n_tree = -10
    ), "\'n_tree\' must be a positive integer.")

    expect_error(detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        anomaly_threshold = 1.5
    ), "\'anomaly_threshold\' must be a positive number greater than 0 and less than 1.")
})

test_that("detectAnomaly works with custom cell types and pc_subset", {
    # Store PCA anomaly data with custom cell types and pc_subset
    anomaly_output <- detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        cell_types = c("CD4", "CD8"),
        pc_subset = 1:5
    )

    # Check if the output is a list
    expect_true(is.list(anomaly_output))

    # Check if the output contains the expected elements
    expect_true("CD4" %in% names(anomaly_output))
    expect_true("CD8" %in% names(anomaly_output))

    # Check the structure of the elements for a specific cell type
    expect_true(all(c("reference_anomaly_scores", "reference_anomaly", "reference_mat_subset", "query_mat_subset", "query_anomaly_scores", "query_anomaly", "var_explained") %in% names(anomaly_output[["CD4"]])))
})

