# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Create test regression results for different scenarios
regress_res_query_only <- regressPC(query_data = query_data,
                                    query_cell_type_col = "expert_annotation",
                                    cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                                    pc_subset = 1:5)

regress_res_combined <- regressPC(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = "expert_annotation",
                                  ref_cell_type_col = "expert_annotation",
                                  cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                                  pc_subset = 1:5)

# Add batch information for batch testing
query_data_batch <- query_data
query_data_batch$batch <- sample(c("batch1", "batch2"), ncol(query_data_batch), replace = TRUE)

regress_res_batch <- regressPC(query_data = query_data_batch,
                               query_cell_type_col = "expert_annotation",
                               query_batch_col = "batch",
                               cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                               pc_subset = 1:5)

# Add multiple batch scenario for query+reference
reference_data_batch <- reference_data
query_data_multi_batch <- query_data
query_data_multi_batch$batch <- sample(c("batch1", "batch2", "batch3"), ncol(query_data_multi_batch), replace = TRUE)

regress_res_multi_batch <- regressPC(query_data = query_data_multi_batch,
                                     reference_data = reference_data_batch,
                                     query_cell_type_col = "expert_annotation",
                                     ref_cell_type_col = "expert_annotation",
                                     query_batch_col = "batch",
                                     cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                                     pc_subset = 1:5)

test_that("plot.regressPCObject generates r_squared plots correctly", {
    # Test query-only analysis
    p1 <- plot(regress_res_query_only, plot_type = "r_squared")
    expect_s3_class(p1, "ggplot")
    expect_true("GeomCol" %in% class(p1$layers[[1]]$geom))

    # Test combined analysis (should have stacked bars with components)
    p2 <- plot(regress_res_combined, plot_type = "r_squared")
    expect_s3_class(p2, "ggplot")
    expect_true("GeomCol" %in% class(p2$layers[[1]]$geom))

    # Test batch analysis (should have stacked bars with components)
    p3 <- plot(regress_res_batch, plot_type = "r_squared")
    expect_s3_class(p3, "ggplot")
    expect_true("GeomCol" %in% class(p3$layers[[1]]$geom))

    # Test multi-batch analysis
    p4 <- plot(regress_res_multi_batch, plot_type = "r_squared")
    expect_s3_class(p4, "ggplot")
    expect_true("GeomCol" %in% class(p4$layers[[1]]$geom))
})

test_that("plot.regressPCObject generates variance_contribution plots correctly", {
    # Test query-only analysis
    p1 <- plot(regress_res_query_only, plot_type = "variance_contribution")
    expect_s3_class(p1, "ggplot")
    expect_true("GeomCol" %in% class(p1$layers[[1]]$geom))

    # Test combined analysis (should have stacked bars with components)
    p2 <- plot(regress_res_combined, plot_type = "variance_contribution")
    expect_s3_class(p2, "ggplot")
    expect_true("GeomCol" %in% class(p2$layers[[1]]$geom))

    # Test batch analysis (should have stacked bars with components)
    p3 <- plot(regress_res_batch, plot_type = "variance_contribution")
    expect_s3_class(p3, "ggplot")
    expect_true("GeomCol" %in% class(p3$layers[[1]]$geom))

    # Test multi-batch analysis
    p4 <- plot(regress_res_multi_batch, plot_type = "variance_contribution")
    expect_s3_class(p4, "ggplot")
    expect_true("GeomCol" %in% class(p4$layers[[1]]$geom))
})

test_that("plot.regressPCObject generates coefficient_heatmap plots correctly", {
    # Test query-only analysis
    p1 <- plot(regress_res_query_only, plot_type = "coefficient_heatmap")
    expect_s3_class(p1, "ggplot")
    expect_true("GeomTile" %in% class(p1$layers[[1]]$geom))

    # Test combined analysis
    p2 <- plot(regress_res_combined, plot_type = "coefficient_heatmap")
    expect_s3_class(p2, "ggplot")
    expect_true("GeomTile" %in% class(p2$layers[[1]]$geom))

    # Test batch analysis
    p3 <- plot(regress_res_batch, plot_type = "coefficient_heatmap")
    expect_s3_class(p3, "ggplot")
    expect_true("GeomTile" %in% class(p3$layers[[1]]$geom))

    # Test multi-batch analysis
    p4 <- plot(regress_res_multi_batch, plot_type = "coefficient_heatmap")
    expect_s3_class(p4, "ggplot")
    expect_true("GeomTile" %in% class(p4$layers[[1]]$geom))
})

test_that("plot.regressPCObject handles default plot type correctly", {
    # Test default plot type (should be r_squared)
    p1 <- plot(regress_res_query_only)
    expect_s3_class(p1, "ggplot")

    p2 <- plot(regress_res_combined)
    expect_s3_class(p2, "ggplot")

    p3 <- plot(regress_res_batch)
    expect_s3_class(p3, "ggplot")
})

test_that("plot.regressPCObject handles invalid plot types", {
    # Test invalid plot type
    expect_error(plot(regress_res_query_only, plot_type = "invalid_type"))
    expect_error(plot(regress_res_combined, plot_type = "invalid_type"))
    expect_error(plot(regress_res_batch, plot_type = "nonexistent_plot"))
})

test_that("plot.regressPCObject respects alpha parameter for heatmaps", {
    # Test with different alpha values
    p1 <- plot(regress_res_query_only, plot_type = "coefficient_heatmap", alpha = 0.01)
    expect_s3_class(p1, "ggplot")

    p2 <- plot(regress_res_combined, plot_type = "coefficient_heatmap", alpha = 0.1)
    expect_s3_class(p2, "ggplot")

    p3 <- plot(regress_res_batch, plot_type = "coefficient_heatmap", alpha = 0.001)
    expect_s3_class(p3, "ggplot")
})

test_that("plot.regressPCObject generates correct plot elements for stacked bars", {
    # Check that component breakdown plots have legends for stacked bars
    p1 <- plot(regress_res_combined, plot_type = "r_squared")
    expect_s3_class(p1, "ggplot")

    # Extract plot data to verify stacking
    built_plot <- ggplot2::ggplot_build(p1)
    expect_true(length(built_plot$data) >= 1)

    # Test variance contribution stacked bars
    p2 <- plot(regress_res_batch, plot_type = "variance_contribution")
    expect_s3_class(p2, "ggplot")

    built_plot2 <- ggplot2::ggplot_build(p2)
    expect_true(length(built_plot2$data) >= 1)
})

test_that("plot.regressPCObject handles edge cases gracefully", {
    # Test with minimal PC subset
    regress_res_minimal <- regressPC(query_data = query_data,
                                     query_cell_type_col = "expert_annotation",
                                     cell_types = c("CD4", "CD8"),
                                     pc_subset = 1:2)

    p1 <- plot(regress_res_minimal, plot_type = "r_squared")
    expect_s3_class(p1, "ggplot")

    p2 <- plot(regress_res_minimal, plot_type = "variance_contribution")
    expect_s3_class(p2, "ggplot")

    p3 <- plot(regress_res_minimal, plot_type = "coefficient_heatmap")
    expect_s3_class(p3, "ggplot")
})

test_that("plot.regressPCObject component decomposition works correctly", {
    # Test that objects with component decomposition create stacked bars
    if (!is.null(regress_res_combined[["r_squared_components"]]) &&
        length(regress_res_combined[["r_squared_components"]]) > 1) {

        p1 <- plot(regress_res_combined, plot_type = "r_squared")
        expect_s3_class(p1, "ggplot")

        # Check that the plot has fill aesthetic (indicating stacked bars)
        expect_true("fill" %in% names(p1$mapping))
    }

    # Test batch interaction component decomposition
    if (!is.null(regress_res_batch[["r_squared_components"]]) &&
        length(regress_res_batch[["r_squared_components"]]) > 1) {

        p2 <- plot(regress_res_batch, plot_type = "variance_contribution")
        expect_s3_class(p2, "ggplot")

        # Check that the plot has fill aesthetic (indicating stacked bars)
        expect_true("fill" %in% names(p2$mapping))
    }
})

test_that("plot.regressPCObject handles different model types in heatmaps", {
    # Test cell type only model
    p1 <- plot(regress_res_query_only, plot_type = "coefficient_heatmap")
    expect_s3_class(p1, "ggplot")

    # Test cell type * dataset interaction
    p2 <- plot(regress_res_combined, plot_type = "coefficient_heatmap")
    expect_s3_class(p2, "ggplot")

    # Test cell type * batch interaction
    p3 <- plot(regress_res_batch, plot_type = "coefficient_heatmap")
    expect_s3_class(p3, "ggplot")

    # Test multi-batch scenario
    p4 <- plot(regress_res_multi_batch, plot_type = "coefficient_heatmap")
    expect_s3_class(p4, "ggplot")
})

test_that("plot.regressPCObject produces plots with correct titles and labels", {
    # Test titles for different plot types
    p1 <- plot(regress_res_query_only, plot_type = "r_squared")
    expect_true(grepl("R-squared Values", p1$labels$title))

    p2 <- plot(regress_res_combined, plot_type = "variance_contribution")
    expect_true(grepl("Variance Contribution", p2$labels$title))

    p3 <- plot(regress_res_batch, plot_type = "coefficient_heatmap")
    expect_true(grepl("Regression Coefficients", p3$labels$title))

    # Test that y-axis labels are correct
    expect_true(grepl("R", as.character(p1$labels$y)) || is.expression(p1$labels$y))
    expect_true(grepl("Variance", p2$labels$y))
})

test_that("plot.regressPCObject handles missing or NULL components gracefully", {
    # Create a minimal regression object without component decomposition
    minimal_regress <- list(
        r_squared = c(PC1 = 0.5, PC2 = 0.3),
        var_contributions = c(PC1 = 2.5, PC2 = 1.5),
        total_variance_explained = 4.0,
        query_pca_var = c(5.0, 5.0, 4.0, 4.0, 3.0),
        indep_var = "cell_type",
        reference_cell_type = "CD4",
        regression_summaries = list(
            PC1 = list(coefficients = data.frame(coef = 0.1, se = 0.05, t = 2.0, p.value = 0.05, p.adjusted = 0.1)),
            PC2 = list(coefficients = data.frame(coef = 0.2, se = 0.06, t = 3.0, p.value = 0.01, p.adjusted = 0.02))
        )
    )
    class(minimal_regress) <- c("list", "regressPCObject")

    # Test that plots work even without component decomposition
    p1 <- plot(minimal_regress, plot_type = "r_squared")
    expect_s3_class(p1, "ggplot")

    p2 <- plot(minimal_regress, plot_type = "variance_contribution")
    expect_s3_class(p2, "ggplot")
})

test_that("plot.regressPCObject additional arguments are passed correctly", {
    # Test that additional arguments don't break the function
    p1 <- plot(regress_res_query_only, plot_type = "r_squared", extra_arg = "test")
    expect_s3_class(p1, "ggplot")

    p2 <- plot(regress_res_combined, plot_type = "coefficient_heatmap", alpha = 0.05, another_arg = TRUE)
    expect_s3_class(p2, "ggplot")
})

