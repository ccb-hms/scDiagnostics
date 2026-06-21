# Load necessary libraries
library(testthat)
library(scDiagnostics)
library(ggplot2)

# Load example datasets
data("reference_data")
data("query_data")

# Create a helper object to test the plotting function
# Skipping at this level if scran is missing to avoid test failure during object creation
if (requireNamespace("scran", quietly = TRUE)) {
    recon_output <- calculateReconstructionError(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        pc_subset = 1:3,
        n_hvgs = 50,
        mad_multiplier = 2
    )

    recon_output_ref_only <- calculateReconstructionError(
        reference_data = reference_data,
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_hvgs = 50
    )
}

test_that("plot.calculateReconstructionErrorObject generates valid ggplot objects", {
    skip_if_not_installed("scran")

    # Test Violin plot (Default)
    p_violin <- plot(recon_output, plot_type = "violin")
    expect_s3_class(p_violin, "ggplot")

    # Test Boxplot
    p_boxplot <- plot(recon_output, plot_type = "boxplot")
    expect_s3_class(p_boxplot, "ggplot")

    # Test Ridge plot
    p_ridge <- plot(recon_output, plot_type = "ridge")
    expect_s3_class(p_ridge, "ggplot")
})

test_that("plot.calculateReconstructionErrorObject generates valid Heatmap objects", {
    skip_if_not_installed("scran")
    skip_if_not_installed("ComplexHeatmap")
    skip_if_not_installed("circlize")

    # Test Heatmap (Both datasets)
    h_both <- plot(recon_output, plot_type = "heatmap", data_type = "both", draw_plot = FALSE)
    expect_s4_class(h_both, "Heatmap")

    # Test Heatmap (Query only)
    h_query <- plot(recon_output, plot_type = "heatmap", data_type = "query", draw_plot = FALSE)
    expect_s4_class(h_query, "Heatmap")

    # Test Heatmap (Reference only)
    h_ref <- plot(recon_output_ref_only, plot_type = "heatmap", data_type = "reference", draw_plot = FALSE)
    expect_s4_class(h_ref, "Heatmap")
})

test_that("plot handles data_type argument correctly for distribution plots", {
    skip_if_not_installed("scran")

    # Test plotting query only
    p_query <- plot(recon_output, data_type = "query", plot_type = "boxplot")
    expect_s3_class(p_query, "ggplot")

    # Test plotting reference only
    p_ref <- plot(recon_output, data_type = "reference", plot_type = "boxplot")
    expect_s3_class(p_ref, "ggplot")

    # Test plotting both (default)
    p_both <- plot(recon_output, data_type = "both", plot_type = "boxplot")
    expect_s3_class(p_both, "ggplot")
})

test_that("plot handles specific cell types and defaults", {
    skip_if_not_installed("scran")

    # Should default to Combined if available
    p_default <- plot(recon_output)
    expect_s3_class(p_default, "ggplot")

    # Test explicit cell type
    first_cell_type <- names(recon_output)[names(recon_output) != "Combined"][1]
    p_explicit <- plot(recon_output, cell_type = first_cell_type)
    expect_s3_class(p_explicit, "ggplot")
})

test_that("plot handles error cases correctly", {
    skip_if_not_installed("scran")

    # Test with invalid cell type
    expect_error(
        plot(recon_output, cell_type = "NonexistentType"),
        "is not available in the provided object"
    )

    # Test asking for query data when none exists
    expect_error(
        plot(recon_output_ref_only, data_type = "query"),
        "no query data available"
    )
})
