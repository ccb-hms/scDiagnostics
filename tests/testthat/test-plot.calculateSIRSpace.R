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

test_that("plot.calculateSIRSpace generates plots correctly", {
    # Generate plot using the function
    p1 <- plot(sir_output, plot_type = "scatterplot", sir_subset = 1:6)
    p2 <- plot(sir_output, plot_type = "boxplot", sir_subset = 1:6)
    p3 <- plot(sir_output, plot_type = "varplot", sir_subset = 1:6)

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    expect_s3_class(p3, "ggplot")
})

