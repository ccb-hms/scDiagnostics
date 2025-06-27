# Load necessary libraries
library(testthat)
library(scDiagnostics)
library(ggplot2)

# Load example datasets and prepare data as before
data("reference_data")
data("query_data")

# Generate Wasserstein distance data
wasserstein_data <- calculateWassersteinDistance(
    query_data = query_data,
    reference_data = reference_data,
    query_cell_type_col = "expert_annotation",
    ref_cell_type_col = "expert_annotation",
    pc_subset = 1:5,
    n_resamples = 100
)

test_that("plot.calculateWassersteinDistanceObject returns a ggplot object", {
    plot_obj <- plot(wasserstein_data)
    expect_s3_class(plot_obj, "ggplot")
})

test_that("plot creates ridge plots with correct structure", {
    plot_obj <- plot(wasserstein_data)

    # Check that it has ridge layers
    expect_true(length(plot_obj$layers) > 0)

    # Check correct labels
    expect_equal(plot_obj$labels$title, "Comparison of Wasserstein Distance Distributions by Cell Type")
    expect_equal(plot_obj$labels$x, "Wasserstein Distance")
})

test_that("plot works with specific cell types", {
    available_cell_types <- wasserstein_data$cell_types[1:min(2, length(wasserstein_data$cell_types))]
    plot_obj <- plot(wasserstein_data, cell_types = available_cell_types)
    expect_s3_class(plot_obj, "ggplot")
})

test_that("plot includes probability annotations", {
    plot_obj <- plot(wasserstein_data)

    # Check for text annotations (probability of superiority)
    text_layers <- sapply(plot_obj$layers, function(layer) inherits(layer$geom, "GeomText"))
    expect_true(any(text_layers))
})
