# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")


test_that("plot method works correctly", {

    # Use package data to create a proper calculateGeneShiftsObject
    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20
    )

    # Get available cell types for testing
    available_cell_types <- unique(result$cell_metadata$cell_type)
    test_cell_type <- available_cell_types[1]

    # Test input validation
    expect_error(
        plot(result),
        "cell_type must be specified and be a single character string"
    )

    expect_error(
        plot(result, cell_type = c("CD4", "CD8")),
        "cell_type must be specified and be a single character string"
    )

    expect_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot", n_genes = NULL),
        "For 'boxplot' plot type, 'n_genes' must be specified"
    )

    expect_error(
        plot(result, cell_type = test_cell_type, plot_type = "heatmap",
             n_genes = NULL, significance_threshold = NULL),
        "For 'heatmap' plot type, at least one of 'n_genes' or 'significance_threshold' must be specified"
    )

    expect_error(
        plot(result, cell_type = test_cell_type, n_genes = -5),
        "If not NULL, 'n_genes' must be a positive integer"
    )

    expect_error(
        plot(result, cell_type = test_cell_type, n_genes = 0),
        "If not NULL, 'n_genes' must be a positive integer"
    )

    expect_error(
        plot(result, cell_type = test_cell_type, significance_threshold = -0.1),
        "If not NULL, 'significance_threshold' must be a numeric value between 0 and 1"
    )

    expect_error(
        plot(result, cell_type = test_cell_type, significance_threshold = 1.5),
        "If not NULL, 'significance_threshold' must be a numeric value between 0 and 1"
    )

    # Test nonexistent cell type
    expect_error(
        plot(result, cell_type = "NonexistentType"),
        "Cell type 'NonexistentType' not found in results"
    )

    # Test no available PCs
    expect_error(
        plot(result, cell_type = test_cell_type, pc_subset = 10:15),
        "None of the requested PCs were found in the results object"
    )

    # Test boxplot PC limit warning
    expect_warning(
        {
            # Create result with more PCs to trigger the warning
            result_many_pcs <- calculateGeneShifts(
                query_data = query_data,
                reference_data = reference_data,
                query_cell_type_col = "SingleR_annotation",
                ref_cell_type_col = "expert_annotation",
                pc_subset = 1:8,  # Create with 8 PCs
                n_top_loadings = 20
            )
            plot(result_many_pcs,
                 cell_type = test_cell_type,
                 plot_type = "boxplot",
                 pc_subset = 1:8)  # Request all 8 PCs
        },
        "Boxplot type can display at most 5 PCs. Using the first 5"
    )
})

test_that("plot method argument matching works", {

    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20
    )

    test_cell_type <- unique(result$cell_metadata$cell_type)[1]

    # Test plot_type matching
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot")
    )

    # Test plot_by matching
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_by = "p_adjusted")
    )

    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_by = "top_loading")
    )
})

test_that("plot method returns appropriate objects for boxplot", {

    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20
    )

    test_cell_type <- unique(result$cell_metadata$cell_type)[1]

    # Test boxplot returns ggplot object
    plot_result <- plot(result, cell_type = test_cell_type, plot_type = "boxplot")
    expect_s3_class(plot_result, "ggplot")
})

test_that("plot method object structure validation works", {

    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20
    )

    test_cell_type <- unique(result$cell_metadata$cell_type)[1]

    # Test missing required elements
    incomplete_data <- result
    incomplete_data$expression_data <- NULL

    expect_error(
        plot(incomplete_data, cell_type = test_cell_type),
        "Input object 'x' is missing required elements"
    )
})

test_that("plot method works with different parameter combinations", {

    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20
    )

    test_cell_type <- unique(result$cell_metadata$cell_type)[1]

    # Test boxplot with different plot_by options
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot",
             plot_by = "top_loading")
    )

    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot",
             plot_by = "p_adjusted")
    )

    # Test different n_genes values
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot", n_genes = 5)
    )

    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot", n_genes = 15)
    )
})

test_that("plot method handles different pc_subset configurations", {

    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:5,
        n_top_loadings = 20
    )

    test_cell_type <- unique(result$cell_metadata$cell_type)[1]

    # Test single PC
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot", pc_subset = 1)
    )

    # Test multiple PCs
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot", pc_subset = 1:3)
    )

    # Test different PC combinations
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot", pc_subset = c(1, 3, 5))
    )
})

test_that("plot method handles significance_threshold correctly", {

    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20
    )

    test_cell_type <- unique(result$cell_metadata$cell_type)[1]

    # Test different significance thresholds
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot",
             significance_threshold = 0.01)
    )

    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot",
             significance_threshold = 0.1)
    )

    # Test with NULL significance_threshold
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot",
             significance_threshold = NULL)
    )
})

# Conditional tests for heatmap functionality (only if ComplexHeatmap is available)
test_that("plot method heatmap functionality works when ComplexHeatmap is available", {

    skip_if_not_installed("ComplexHeatmap")

    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20
    )

    test_cell_type <- unique(result$cell_metadata$cell_type)[1]

    # Test heatmap with n_genes only
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "heatmap", n_genes = 5)
    )

    # Test heatmap with significance_threshold only
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "heatmap",
             n_genes = NULL, significance_threshold = 0.05)
    )

    # Test heatmap with both n_genes and significance_threshold
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "heatmap",
             n_genes = 5, significance_threshold = 0.05)
    )
})

test_that("plot method handles empty gene selection gracefully", {

    result <- calculateGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20
    )

    test_cell_type <- unique(result$cell_metadata$cell_type)[1]

    # Test boxplot with very small n_genes (should still work)
    expect_no_error(
        plot(result, cell_type = test_cell_type, plot_type = "boxplot", n_genes = 1)
    )

    # This would only test heatmap empty selection if ComplexHeatmap is available
    skip_if_not_installed("ComplexHeatmap")

    # Test with very restrictive significance threshold for heatmap
    expect_message(
        heatmap_result <- plot(result, cell_type = test_cell_type, plot_type = "heatmap",
                               n_genes = NULL, significance_threshold = 0.0001),
        "No data available to generate a heatmap"
    )
    expect_null(heatmap_result)
})

