# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# 1. Store PCA anomaly data
anomaly_output_pca <- detectAnomaly(reference_data = reference_data,
                                    query_data = query_data,
                                    ref_cell_type_col = "expert_annotation",
                                    query_cell_type_col = "SingleR_annotation",
                                    pc_subset = 1:5,
                                    n_tree = 50, # Lowered for faster testing
                                    threshold_method = "MAD",
                                    mad_multiplier = 2)

# 2. Store HVG anomaly data
anomaly_output_hvg <- detectAnomaly(reference_data = reference_data,
                                    query_data = query_data,
                                    ref_cell_type_col = "expert_annotation",
                                    query_cell_type_col = "SingleR_annotation",
                                    pc_subset = NULL,
                                    n_hvgs = 20, # Low number for faster testing
                                    n_tree = 50,
                                    threshold_method = "MAD",
                                    mad_multiplier = 2)

test_that("plot.detectAnomalyObject generates PCA scatter plots correctly", {
    # Generate plot using the function (query data)
    p1 <- plot(anomaly_output_pca,
               cell_type = "CD4",
               pc_subset = 1:3,
               data_type = "query")

    # Check if output is a ggmatrix object (from GGally::ggpairs)
    expect_s3_class(p1, "ggmatrix")

    # Generate plot using the function (reference data)
    p2 <- plot(anomaly_output_pca,
               cell_type = "CD4",
               pc_subset = 1:3,
               data_type = "reference")

    # Check if output is a ggmatrix object
    expect_s3_class(p2, "ggmatrix")
})

test_that("plot.detectAnomalyObject generates HVG heatmaps correctly", {
    # Skip if ComplexHeatmap is not available
    skip_if_not_installed("ComplexHeatmap")

    # Generate heatmap (query data)
    h1 <- plot(anomaly_output_hvg,
               cell_type = "CD4",
               data_type = "query",
               draw_plot = FALSE) # Prevent drawing to keep tests fast

    expect_s4_class(h1, "Heatmap")

    # Generate heatmap (both data)
    h2 <- plot(anomaly_output_hvg,
               cell_type = "CD4",
               data_type = "both",
               draw_plot = FALSE)

    expect_s4_class(h2, "Heatmap")

    # Generate heatmap (reference data)
    h3 <- plot(anomaly_output_hvg,
               cell_type = "CD4",
               data_type = "reference",
               draw_plot = FALSE)

    expect_s4_class(h3, "Heatmap")
})

test_that("plot works with different facet options for PCA", {
    # Test different diagonal facet options
    p1 <- plot(anomaly_output_pca,
               cell_type = "CD4",
               pc_subset = 1:2,
               data_type = "query",
               diagonal_facet = "density",
               upper_facet = "blank")

    expect_s3_class(p1, "ggmatrix")

    # Test with ridge diagonal and contour upper
    p2 <- plot(anomaly_output_pca,
               cell_type = "CD4",
               pc_subset = 1:2,
               data_type = "query",
               diagonal_facet = "ridge",
               upper_facet = "contour")

    expect_s3_class(p2, "ggmatrix")

    # Test with boxplot diagonal and ellipse upper
    p3 <- plot(anomaly_output_pca,
               cell_type = "CD4",
               pc_subset = 1:2,
               data_type = "query",
               diagonal_facet = "boxplot",
               upper_facet = "ellipse")

    expect_s3_class(p3, "ggmatrix")
})

test_that("plot works with Combined cell type and default parameters", {
    # Test with Combined cell type (default)
    p1 <- plot(anomaly_output_pca,
               pc_subset = 1:2,
               data_type = "query")

    expect_s3_class(p1, "ggmatrix")

    # Test with explicit Combined cell type
    p2 <- plot(anomaly_output_pca,
               cell_type = "Combined",
               pc_subset = 1:2,
               data_type = "reference")

    expect_s3_class(p2, "ggmatrix")
})

test_that("plot handles error cases correctly", {
    # Test with invalid cell type
    expect_error(plot(anomaly_output_pca,
                      cell_type = "NonexistentType",
                      pc_subset = 1:2),
                 "cell_type.*not available")

    # Test with out of range pc_subset
    max_pc <- ncol(anomaly_output_pca$CD4$reference_mat_subset)
    expect_error(plot(anomaly_output_pca,
                      cell_type = "CD4",
                      pc_subset = 1:(max_pc + 5)),
                 "pc_subset.*out of range")

    # Test passing data_type = "both" to PCA object
    expect_error(plot(anomaly_output_pca,
                      cell_type = "CD4",
                      data_type = "both"),
                 "only supported for HVG heatmaps")
})

test_that("plot works with different n_tree parameter", {
    # Test with different n_tree parameter
    p1 <- plot(anomaly_output_pca,
               cell_type = "CD4",
               pc_subset = 1:2,
               data_type = "query",
               n_tree = 100)

    expect_s3_class(p1, "ggmatrix")
})
