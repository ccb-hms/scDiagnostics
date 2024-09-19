# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Project the query data onto PCA space of reference
nn_output <- calculateNearestNeighborProbabilities(query_data = query_data, 
                                                   reference_data = reference_data,
                                                   query_cell_type_col = "SingleR_annotation", 
                                                   ref_cell_type_col = "expert_annotation",
                                                   n_neighbor = 20, 
                                                   pc_subset = 1:5)

test_that("plot.calculateNearestNeighborProbabilities generates plots correctly", {
    # Generate plot using the function (no query data)
    p1 <- plot(nn_output)

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
})
