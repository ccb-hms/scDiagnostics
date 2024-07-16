# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load data
data("reference_data")
data("query_data")

# Compute pairwise correlations
cor_matrix_avg <- calculateAveragePairwiseCorrelation(query_data = query_data, 
                                                      reference_data = reference_data, 
                                                      query_cell_type_col = "expert_annotation", 
                                                      ref_cell_type_col = "expert_annotation", 
                                                      cell_types = c("CD4", "CD8", "B_and_plasma"), 
                                                      pc_subset = 1:5,
                                                      correlation_method = "spearman")

test_that("plot.calculateAveragePairwiseCorrelation generates plots correctly", {
    # Generate plot using the function
    p1 <- plot(cor_matrix_avg)
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
})
