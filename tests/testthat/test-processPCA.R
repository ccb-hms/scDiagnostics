# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load the datasets that come with the package
data("query_data", package = "scDiagnostics")
data("reference_data", package = "scDiagnostics")

test_that("processPCA works with single dataset without existing PCA", {
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")

    # Create test dataset without PCA
    test_data <- query_data
    reducedDims(test_data) <- list()  # Remove all existing PCA

    # Test single dataset without PCA
    result <- processPCA(sce_object = test_data, n_hvgs = 300)

    expect_s4_class(result, "SingleCellExperiment")
    expect_true("PCA" %in% reducedDimNames(result))
    expect_equal(ncol(result), 503)  # Should preserve all cells

    # Check PCA attributes exist
    pca_attrs <- attributes(reducedDim(result, "PCA"))
    expect_true("rotation" %in% names(pca_attrs))
    expect_true("percentVar" %in% names(pca_attrs))
})

test_that("processPCA preserves existing valid PCA", {
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")

    # Use dataset that already has PCA
    test_data <- query_data
    original_pca <- reducedDim(test_data, "PCA")

    # Process data that already has valid PCA
    result <- processPCA(sce_object = test_data, n_hvgs = 300)

    expect_s4_class(result, "SingleCellExperiment")
    expect_true("PCA" %in% reducedDimNames(result))
    expect_equal(ncol(result), 503)  # Should preserve all cells

    # Should return the same object unchanged
    expect_identical(reducedDim(result, "PCA"), original_pca)
})

test_that("processPCA works with larger dataset", {
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")

    # Use reference dataset (larger) without PCA
    test_data <- reference_data
    reducedDims(test_data) <- list()  # Remove existing PCA

    # Test with max_cells limit
    result <- processPCA(sce_object = test_data, n_hvgs = 500, max_cells = 800)

    expect_s4_class(result, "SingleCellExperiment")
    expect_true("PCA" %in% reducedDimNames(result))
    expect_equal(ncol(result), 800)  # Should be downsampled to max_cells
})

test_that("processPCA preserves large dataset with valid PCA", {
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")

    # Use reference dataset that already has PCA
    test_data <- reference_data

    # Process data that already has valid PCA - should NOT downsample
    result <- processPCA(sce_object = test_data, n_hvgs = 500, max_cells = 800)

    expect_s4_class(result, "SingleCellExperiment")
    expect_true("PCA" %in% reducedDimNames(result))
    expect_equal(ncol(result), 1500)  # Should preserve all 1500 cells (no downsampling)
})

test_that("processPCA input validation", {
    # Test invalid input
    expect_error(processPCA("not_an_sce"),
                 "'sce_object' must be a SingleCellExperiment object")

    expect_error(processPCA(query_data, assay_name = "nonexistent"),
                 "'sce_object' does not contain the specified assay")

    expect_error(processPCA(query_data, n_hvgs = -1),
                 "'n_hvgs' must be a single positive integer")

    expect_error(processPCA(query_data, max_cells = -1),
                 "'max_cells' must be a positive integer")
})
