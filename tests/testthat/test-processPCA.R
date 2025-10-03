# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load test data once at the top level
skip_if_not_installed("TENxPBMCData")
skip_if_not_installed("scater")
skip_if_not_installed("scran")

# Create shared test datasets
pbmc_data <- NULL
test_query <- NULL
test_ref <- NULL
test_large <- NULL

setup_test_data <- function() {
    if (is.null(pbmc_data)) {
        library(TENxPBMCData)
        pbmc_data <<- TENxPBMCData("pbmc3k")

        # Create test datasets with reasonable sizes to avoid PCA warnings
        test_query <<- pbmc_data[, 1:300]  # Larger to avoid SVD warnings
        test_query <<- scater::logNormCounts(test_query)

        test_ref <<- pbmc_data[, 301:600]
        test_ref <<- scater::logNormCounts(test_ref)

        test_large <<- pbmc_data[, 1:800]
        test_large <<- scater::logNormCounts(test_large)
    }
}

test_that("processPCA works with single dataset", {
    setup_test_data()

    # Create test dataset without PCA
    test_data <- test_query
    reducedDims(test_data) <- list()

    # Test single dataset without PCA
    result <- processPCA(sce_object = test_data, n_hvgs = 1000)

    expect_s4_class(result, "SingleCellExperiment")
    expect_true("PCA" %in% reducedDimNames(result))
    expect_equal(ncol(result), 300)
})

test_that("processPCA preserves existing PCA", {
    setup_test_data()

    # Create dataset with existing PCA
    test_data <- scater::runPCA(test_query, ncomponents = 10)
    original_cells <- ncol(test_data)

    # Test with existing PCA - should not downsample
    result <- processPCA(sce_object = test_data, max_cells = 100)

    expect_s4_class(result, "SingleCellExperiment")
    expect_true("PCA" %in% reducedDimNames(result))
    expect_equal(ncol(result), original_cells)  # Should be unchanged
})

test_that("processPCA handles downsampling correctly", {
    setup_test_data()

    # Use large dataset without PCA
    test_data <- test_large
    reducedDims(test_data) <- list()

    # Test downsampling
    result <- processPCA(sce_object = test_data, max_cells = 400, n_hvgs = 1000)

    expect_s4_class(result, "SingleCellExperiment")
    expect_true("PCA" %in% reducedDimNames(result))
    expect_equal(ncol(result), 400)  # Should be downsampled
})

test_that("processPCA input validation", {
    expect_error(
        processPCA()
    )
})

