# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Store PCA anomaly data
anomaly_output <- detectAnomaly(reference_data = reference_data,
                                query_data = query_data,
                                ref_cell_type_col = "expert_annotation",
                                query_cell_type_col = "SingleR_annotation",
                                pc_subset = 1:10,
                                n_tree = 500,
                                anomaly_threshold = 0.5)

test_that("plot.detectAnomalyObject generates plots correctly", {
    # Generate plot using the function (query data)
    p1 <- plot(anomaly_output,
               cell_type = "CD4",
               pc_subset = 1:3,
               data_type = "query")

    # Check if output is a ggmatrix object (from GGally::ggpairs)
    expect_s3_class(p1, "ggmatrix")

    # Generate plot using the function (reference data)
    p2 <- plot(anomaly_output,
               cell_type = "CD4",
               pc_subset = 1:3,
               data_type = "reference")

    # Check if output is a ggmatrix object
    expect_s3_class(p2, "ggmatrix")
})

test_that("plot works with different facet options", {
    # Test different diagonal facet options
    p1 <- plot(anomaly_output,
               cell_type = "CD4",
               pc_subset = 1:2,
               data_type = "query",
               diagonal_facet = "density",
               upper_facet = "blank")

    expect_s3_class(p1, "ggmatrix")

    # Test with ridge diagonal and contour upper
    p2 <- plot(anomaly_output,
               cell_type = "CD4",
               pc_subset = 1:2,
               data_type = "query",
               diagonal_facet = "ridge",
               upper_facet = "contour")

    expect_s3_class(p2, "ggmatrix")

    # Test with boxplot diagonal and ellipse upper
    p3 <- plot(anomaly_output,
               cell_type = "CD4",
               pc_subset = 1:2,
               data_type = "query",
               diagonal_facet = "boxplot",
               upper_facet = "ellipse")

    expect_s3_class(p3, "ggmatrix")
})

test_that("plot works with Combined cell type and default parameters", {
    # Test with Combined cell type (default)
    p1 <- plot(anomaly_output,
               pc_subset = 1:2,
               data_type = "query")

    expect_s3_class(p1, "ggmatrix")

    # Test with explicit Combined cell type
    p2 <- plot(anomaly_output,
               cell_type = "Combined",
               pc_subset = 1:2,
               data_type = "reference")

    expect_s3_class(p2, "ggmatrix")
})

test_that("plot handles error cases correctly", {
    # Test with invalid cell type
    expect_error(plot(anomaly_output,
                      cell_type = "NonexistentType",
                      pc_subset = 1:2),
                 "cell_type.*not available")

    # Test with out of range pc_subset
    max_pc <- ncol(anomaly_output$CD4$reference_mat_subset)
    expect_error(plot(anomaly_output,
                      cell_type = "CD4",
                      pc_subset = 1:(max_pc + 5)),
                 "pc_subset.*out of range")
})

test_that("plot works with different n_tree parameter", {
    # Test with different n_tree parameter
    p1 <- plot(anomaly_output,
               cell_type = "CD4",
               pc_subset = 1:2,
               data_type = "query",
               n_tree = 100)

    expect_s3_class(p1, "ggmatrix")
})
