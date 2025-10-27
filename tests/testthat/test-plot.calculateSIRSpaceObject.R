# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Compute important variables for all pairwise cell comparisons
sir_output <- calculateSIRSpace(reference_data = reference_data,
                                query_data = query_data,
                                query_cell_type_col = "SingleR_annotation",
                                ref_cell_type_col = "expert_annotation",
                                multiple_cond_means = TRUE,
                                cumulative_variance_threshold = 0.9,
                                n_neighbor = 1)

test_that("plot.calculateSIRSpaceObject generates scores plots correctly", {
    # Generate scores plot using the function (default plot_type)
    p1 <- plot(sir_output, plot_type = "scores", sir_subset = 1:3)

    # Check if output is a ggmatrix object (from GGally::ggpairs)
    expect_s3_class(p1, "ggmatrix")
})

test_that("plot.calculateSIRSpaceObject generates loadings plots correctly", {
    # Generate loadings plot using the function
    p1 <- plot(sir_output, plot_type = "loadings", sir_subset = 1:3, n_top = 5)

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
})

test_that("plot.calculateSIRSpaceObject works with different facet options", {
    # Test different facet combinations for scores plot
    p1 <- plot(sir_output,
               plot_type = "scores",
               sir_subset = 1:2,
               lower_facet = "scatter",
               diagonal_facet = "density",
               upper_facet = "blank")

    expect_s3_class(p1, "ggmatrix")

    # Test with different options
    p2 <- plot(sir_output,
               plot_type = "scores",
               sir_subset = 1:2,
               lower_facet = "contour",
               diagonal_facet = "boxplot",
               upper_facet = "ellipse")

    expect_s3_class(p2, "ggmatrix")
})

test_that("plot.calculateSIRSpaceObject works with specific cell types", {
    # Test with specific cell types
    available_cell_types <- unique(sir_output$sir_projections$cell_type)
    selected_cell_types <- available_cell_types[1:min(2, length(available_cell_types))]

    p1 <- plot(sir_output,
               plot_type = "scores",
               cell_types = selected_cell_types,
               sir_subset = 1:2)

    expect_s3_class(p1, "ggmatrix")
})

test_that("plot.calculateSIRSpaceObject parameter validation works", {
    # Test invalid n_top for loadings
    expect_error(plot(sir_output, plot_type = "loadings", n_top = -1),
                 "n_top must be a positive integer")

    # Test invalid sir_subset (assuming there are fewer SIR components available)
    max_sir <- ncol(sir_output$rotation_mat)
    expect_error(plot(sir_output, sir_subset = 1:(max_sir + 5)),
                 "sir_subset contains values outside the range")
})
