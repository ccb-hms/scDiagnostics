# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example dataset
data("query_data")

test_that("plotGeneExpressionDimred generates plots correctly", {
    # Test PCA method - returns ggmatrix object
    p1 <- plotGeneExpressionDimred(
        sce_object = query_data,
        cell_type_col = "SingleR_annotation",
        method = "PCA",
        pc_subset = 1:5,
        feature = "VPREB3"
    )

    # Check if output is a ggmatrix object (from GGally) for PCA
    expect_s3_class(p1, "ggmatrix")

    # Test UMAP method - should return ggplot object
    # (assuming query_data has UMAP coordinates)
    if ("UMAP" %in% reducedDimNames(query_data)) {
        p2 <- plotGeneExpressionDimred(
            sce_object = query_data,
            cell_type_col = "SingleR_annotation",
            method = "UMAP",
            feature = "VPREB3"
        )
        expect_s3_class(p2, "ggplot")
    }

    # Test TSNE method - should return ggplot object
    # (assuming query_data has TSNE coordinates)
    if ("TSNE" %in% reducedDimNames(query_data)) {
        p3 <- plotGeneExpressionDimred(
            sce_object = query_data,
            cell_type_col = "SingleR_annotation",
            method = "TSNE",
            feature = "VPREB3"
        )
        expect_s3_class(p3, "ggplot")
    }
})

test_that("plotGeneExpressionDimred handles different assay names", {
    # Test with different assay name (if available)
    available_assays <- assayNames(query_data)

    if ("counts" %in% available_assays) {
        p1 <- plotGeneExpressionDimred(
            sce_object = query_data,
            cell_type_col = "SingleR_annotation",
            method = "PCA",
            pc_subset = 1:3,
            feature = "VPREB3",
            assay_name = "counts"
        )
        expect_s3_class(p1, "ggmatrix")
    }

    # Test with invalid assay name
    expect_error(
        plotGeneExpressionDimred(
            sce_object = query_data,
            cell_type_col = "SingleR_annotation",
            method = "PCA",
            pc_subset = 1:3,
            feature = "VPREB3",
            assay_name = "invalid_assay"
        )
    )
})

