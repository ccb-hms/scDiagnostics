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

    # Check if the output is a list and has correct class
    expect_true(is.list(anomaly_output))
    expect_s3_class(anomaly_output, "detectAnomalyObject")

    # Check if the output contains the expected elements
    expect_true(length(anomaly_output) > 0)
    expect_true("Combined" %in% names(anomaly_output))

    # Check the structure of the elements for the first cell type
    first_cell_type <- names(anomaly_output)[1]
    expect_true(all(c("reference_anomaly_scores", "reference_anomaly", "reference_mat_subset") %in% names(anomaly_output[[first_cell_type]])))

    # Check query-related components if query data is provided
    if (!is.null(query_data)) {
        expect_true(all(c("query_mat_subset", "query_anomaly_scores", "query_anomaly") %in% names(anomaly_output[[first_cell_type]])))
    }
})

test_that("detectAnomaly works correctly without query data", {
    # Store PCA anomaly data without query data
    anomaly_output <- detectAnomaly(
        reference_data = reference_data,
        ref_cell_type_col = "expert_annotation"
    )

    # Check if the output is a list and has correct class
    expect_true(is.list(anomaly_output))
    expect_s3_class(anomaly_output, "detectAnomalyObject")

    # Check if the output contains the expected elements
    expect_true(length(anomaly_output) > 0)
    expect_true("Combined" %in% names(anomaly_output))

    # Check the structure of the elements for the first cell type
    first_cell_type <- names(anomaly_output)[1]
    expect_true(all(c("reference_anomaly_scores", "reference_anomaly", "reference_mat_subset") %in% names(anomaly_output[[first_cell_type]])))

    # Check that query-related components are not present
    expect_false("query_anomaly_scores" %in% names(anomaly_output[[first_cell_type]]))
    expect_false("query_anomaly" %in% names(anomaly_output[[first_cell_type]]))
})

test_that("detectAnomaly handles parameter validation correctly", {
    # Test invalid n_tree
    expect_error(detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        n_tree = -10
    ), "n_tree.*must be a positive integer")

    expect_error(detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        n_tree = 0
    ), "n_tree.*must be a positive integer")

    # Test invalid anomaly_threshold
    expect_error(detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        anomaly_threshold = 1.5
    ), "anomaly_threshold.*must be a positive number.*less than 1")

    expect_error(detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        anomaly_threshold = 0
    ), "anomaly_threshold.*must be a positive number.*greater than 0")
})

test_that("detectAnomaly works with custom parameters", {
    # Store PCA anomaly data with custom parameters
    anomaly_output <- detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        cell_types = c("CD4", "CD8"),
        pc_subset = 1:3,
        n_tree = 100,
        anomaly_threshold = 0.7
    )

    # Check if the output is a list and has correct class
    expect_true(is.list(anomaly_output))
    expect_s3_class(anomaly_output, "detectAnomalyObject")

    # Check that it contains the specified cell types plus Combined
    expect_true("CD4" %in% names(anomaly_output))
    expect_true("CD8" %in% names(anomaly_output))
})

test_that("detectAnomaly works without pc_subset (using raw data)", {
    # Test with pc_subset = NULL
    anomaly_output <- detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        pc_subset = NULL,
        n_tree = 50
    )

    # Check basic structure
    expect_true(is.list(anomaly_output))
    expect_s3_class(anomaly_output, "detectAnomalyObject")

    # When pc_subset is NULL, var_explained should not be present
    first_cell_type <- names(anomaly_output)[1]
    expect_false("var_explained" %in% names(anomaly_output[[first_cell_type]]))
})
