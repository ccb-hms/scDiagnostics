# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("calculateReconstructionError works correctly with default parameters", {
    skip_if_not_installed("scran")

    # Run reconstruction error calculation
    recon_output <- calculateReconstructionError(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        pc_subset = 1:3,
        n_hvgs = 50, # low number for speed in testing
        mad_multiplier = 2
    )

    # Check if the output is a list and has correct class
    expect_true(is.list(recon_output))
    expect_s3_class(recon_output, "calculateReconstructionErrorObject")

    # Check if the output contains the expected elements
    expect_true(length(recon_output) > 0)
    expect_true("Combined" %in% names(recon_output))

    # Check the structure of the elements for the first cell type
    first_cell_type <- names(recon_output)[1]

    # NOTE: Added reference_mat_subset and query_mat_subset to expected outputs
    expected_names <- c("reference_reconstruction_errors", "reference_anomaly", "reference_mat_subset",
                        "query_reconstruction_errors", "query_anomaly", "query_mat_subset",
                        "applied_threshold", "var_explained")
    expect_true(all(expected_names %in% names(recon_output[[first_cell_type]])))
})

test_that("calculateReconstructionError works correctly without query data", {
    skip_if_not_installed("scran")

    # Run calculation on reference data only
    recon_output <- calculateReconstructionError(
        reference_data = reference_data,
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_hvgs = 50
    )

    # Check if the output is a list and has correct class
    expect_true(is.list(recon_output))
    expect_s3_class(recon_output, "calculateReconstructionErrorObject")

    # Check the structure of the elements for the first cell type
    first_cell_type <- names(recon_output)[1]
    expect_true("reference_reconstruction_errors" %in% names(recon_output[[first_cell_type]]))
    expect_true("reference_mat_subset" %in% names(recon_output[[first_cell_type]]))

    # Check that query-related components are NOT present
    expect_false("query_reconstruction_errors" %in% names(recon_output[[first_cell_type]]))
    expect_false("query_anomaly" %in% names(recon_output[[first_cell_type]]))
    expect_false("query_mat_subset" %in% names(recon_output[[first_cell_type]]))
})

test_that("calculateReconstructionError handles parameter validation correctly", {
    # Test NULL pc_subset (strictly requires PCA)
    expect_error(calculateReconstructionError(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        pc_subset = NULL
    ), "pc_subset.*cannot be NULL")

    # Test invalid n_hvgs
    expect_error(calculateReconstructionError(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        n_hvgs = -10
    ), "n_hvgs.*must be a single positive integer")

    # Test invalid mad_multiplier
    expect_error(calculateReconstructionError(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        mad_multiplier = -2
    ), "mad_multiplier.*must be a positive numeric")
})
