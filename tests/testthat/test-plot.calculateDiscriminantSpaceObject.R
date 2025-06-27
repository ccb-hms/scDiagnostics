# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Compute important variables for all pairwise cell comparisons
disc_output <- calculateDiscriminantSpace(reference_data = reference_data,
                                          query_data = query_data,
                                          query_cell_type_col = "expert_annotation",
                                          ref_cell_type_col = "expert_annotation",
                                          eigen_threshold  = 1e-1,
                                          n_tree = 500,
                                          n_top = 50,
                                          calculate_metrics = FALSE,
                                          alpha = 0.01)

test_that("plot.calculateDiscriminantSpaceObject generates plots correctly", {
    # Generate plot using the function (returns ggmatrix, not ggplot)
    p1 <- plot(disc_output)

    # Check if output is a ggmatrix object (from GGally::ggpairs)
    expect_s3_class(p1, "ggmatrix")
})

test_that("plot works with different facet options", {
    # Test different facet combinations
    p1 <- plot(disc_output,
               lower_facet = "scatter",
               diagonal_facet = "density",
               upper_facet = "blank")

    expect_s3_class(p1, "ggmatrix")

    # Test with different options
    p2 <- plot(disc_output,
               lower_facet = "contour",
               diagonal_facet = "boxplot",
               upper_facet = "ellipse")

    expect_s3_class(p2, "ggmatrix")
})

test_that("plot works with specific cell types and dv_subset", {
    # Test with specific cell types
    available_cell_types <- unique(c(disc_output$ref_proj$cell_type,
                                     disc_output$query_proj$cell_type))
    selected_cell_types <- available_cell_types[1:min(2, length(available_cell_types))]

    p1 <- plot(disc_output, cell_types = selected_cell_types)
    expect_s3_class(p1, "ggmatrix")

    # Test with dv_subset
    max_dv <- ncol(disc_output$discriminant_eigenvectors)
    dv_subset <- 1:min(2, max_dv)

    p2 <- plot(disc_output, dv_subset = dv_subset)
    expect_s3_class(p2, "ggmatrix")
})

test_that("plot handles error cases correctly", {
    # Create a disc_output without query data
    disc_output_no_query <- calculateDiscriminantSpace(reference_data = reference_data,
                                                       query_data = NULL,
                                                       ref_cell_type_col = "expert_annotation",
                                                       n_tree = 100,
                                                       n_top = 10)

    # Should error when no query data is available
    expect_error(plot(disc_output_no_query), "There is no query data to plot")
})
