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

test_that("plot.regressPCObject generates plots correctly for query-only analysis", {
    # Test r_squared plot
    p1 <- plot(regress_res_query_only, plot_type = "r_squared")
    expect_s3_class(p1, "ggplot")

    # Test heatmap plot
    p2 <- plot(regress_res_query_only, plot_type = "heatmap")
    expect_s3_class(p2, "ggplot")

    # Test variance_contribution plot
    p3 <- plot(regress_res_query_only, plot_type = "variance_contribution")
    expect_s3_class(p3, "ggplot")

    # Test default plot type
    p4 <- plot(regress_res_query_only)
    expect_s3_class(p4, "ggplot")
})

test_that("plot.regressPCObject generates plots correctly for combined analysis", {
    # Test r_squared plot
    p1 <- plot(regress_res_combined, plot_type = "r_squared")
    expect_s3_class(p1, "ggplot")

    # Test heatmap plot
    p2 <- plot(regress_res_combined, plot_type = "heatmap")
    expect_s3_class(p2, "ggplot")

    # Test interaction_effects plot
    p3 <- plot(regress_res_combined, plot_type = "interaction_effects")
    expect_s3_class(p3, "ggplot")

    # Test default plot type
    p4 <- plot(regress_res_combined)
    expect_s3_class(p4, "ggplot")
})

test_that("plot.regressPCObject generates plots correctly for batch analysis", {
    # Test r_squared plot
    p1 <- plot(regress_res_batch, plot_type = "r_squared")
    expect_s3_class(p1, "ggplot")

    # Test heatmap plot
    p2 <- plot(regress_res_batch, plot_type = "heatmap")
    expect_s3_class(p2, "ggplot")

    # Test interaction_effects plot
    p3 <- plot(regress_res_batch, plot_type = "interaction_effects")
    expect_s3_class(p3, "ggplot")
})

test_that("plot.regressPCObject handles invalid plot types", {
    # Test invalid plot type for query-only analysis
    expect_error(plot(regress_res_query_only, plot_type = "invalid_type"))

    # Test invalid plot type for combined analysis
    expect_error(plot(regress_res_combined, plot_type = "invalid_type"))

    # Test unavailable plot type for query-only (interaction_effects)
    expect_error(plot(regress_res_query_only, plot_type = "interaction_effects"))

    # Test unavailable plot type for combined (variance_contribution)
    expect_error(plot(regress_res_combined, plot_type = "variance_contribution"))
})

test_that("plot.regressPCObject respects alpha parameter", {
    # Test with different alpha values
    p1 <- plot(regress_res_query_only, plot_type = "heatmap", alpha = 0.01)
    expect_s3_class(p1, "ggplot")

    p2 <- plot(regress_res_combined, plot_type = "interaction_effects", alpha = 0.1)
    expect_s3_class(p2, "ggplot")
})

test_that("plot.regressPCObject handles edge cases", {
    # Test with minimal PC subset
    regress_res_minimal <- regressPC(query_data = query_data,
                                     query_cell_type_col = "expert_annotation",
                                     cell_types = c("CD4", "CD8"),
                                     pc_subset = 1:2)

    p1 <- plot(regress_res_minimal, plot_type = "r_squared")
    expect_s3_class(p1, "ggplot")

    p2 <- plot(regress_res_minimal, plot_type = "heatmap")
    expect_s3_class(p2, "ggplot")
})
