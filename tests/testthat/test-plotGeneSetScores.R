# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example dataset
data("query_data")

test_that("plotGeneSetScores generates plots correctly", {
    # Generate plot using the function
    p1 <- plotGeneSetScores(
        se_object = query_data,
        method = "PCA",
        score_col = "gene_set_scores",
        pc_subset = 1:5
    )
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    
})

test_that("plotGeneSetScores handles invalid input gracefully", {
    # Test with invalid method
    expect_error(plotGeneSetScores(
        se_object = query_data,
        method = "InvalidMethod",
        score_col = "gene_set_scores",
        pc_subset = 1:5
    ))
    
    # Test with non-existent score_col
    expect_error(plotGeneSetScores(
        se_object = query_data,
        method = "PCA",
        score_col = "invalid_score_column",
        pc_subset = 1:5
    ))
})
