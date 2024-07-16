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

test_that("plot.calculateDiscriminantSpace generates plots correctly", {
    # Generate plot using the function
    p1 <- plot(disc_output, plot_type = "scatterplot")
    p2 <- plot(disc_output, cell_types = "CD4-CD8", plot_type = "boxplot")
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})
