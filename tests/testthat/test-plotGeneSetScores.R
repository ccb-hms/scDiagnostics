# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example dataset
data("query_data")

test_that("plotGeneSetScores generates plots correctly", {
    # Test PCA method - returns ggmatrix object
    p1 <- plotGeneSetScores(
        sce_object = query_data,
        cell_type_col = "SingleR_annotation",
        method = "PCA",
        score_col = "gene_set_scores",
        pc_subset = 1:5
    )

    # Check if output is a ggmatrix object (from GGally) for PCA
    expect_s3_class(p1, "ggmatrix")

    # Test UMAP method - should return ggplot object
    # (assuming query_data has UMAP coordinates)
    if ("UMAP" %in% reducedDimNames(query_data)) {
        p2 <- plotGeneSetScores(
            sce_object = query_data,
            cell_type_col = "SingleR_annotation",
            method = "UMAP",
            score_col = "gene_set_scores"
        )
        expect_s3_class(p2, "ggplot")
    }

    # Test TSNE method - should return ggplot object
    # (assuming query_data has TSNE coordinates)
    if ("TSNE" %in% reducedDimNames(query_data)) {
        p3 <- plotGeneSetScores(
            sce_object = query_data,
            cell_type_col = "SingleR_annotation",
            method = "TSNE",
            score_col = "gene_set_scores"
        )
        expect_s3_class(p3, "ggplot")
    }
})

test_that("plotGeneSetScores handles edge cases", {
    # Test with missing/NA values in scores (if applicable)
    # This would depend on how your function handles NA values in the score column

    # Test with single cell type
    available_cell_types <- unique(colData(query_data)[["SingleR_annotation"]])
    if (length(available_cell_types) > 0) {
        p1 <- plotGeneSetScores(
            sce_object = query_data,
            cell_type_col = "SingleR_annotation",
            method = "PCA",
            score_col = "gene_set_scores",
            pc_subset = 1:2,
            cell_types = available_cell_types[1]
        )
        expect_s3_class(p1, "ggmatrix")
    }
})

test_that("plotGeneSetScores validates score column data", {
    # Check that the score column contains numeric data
    # This assumes gene_set_scores column exists and contains numeric values
    scores <- colData(query_data)[["gene_set_scores"]]
    expect_true(is.numeric(scores))

    # Test plot generation works with actual score data
    p1 <- plotGeneSetScores(
        sce_object = query_data,
        cell_type_col = "SingleR_annotation",
        method = "PCA",
        score_col = "gene_set_scores",
        pc_subset = 1:3
    )
    expect_s3_class(p1, "ggmatrix")
})
