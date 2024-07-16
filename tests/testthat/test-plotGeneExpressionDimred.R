# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example dataset
data("query_data")

test_that("plotGeneExpressionDimred generates plots correctly", {
    # Generate plot using the function with PCA method
    p1 <- plotGeneExpressionDimred(
        se_object = query_data,
        method = "PCA",
        pc_subset = 1:5,
        feature = "VPREB3"
    )
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    
})

test_that("plotGeneExpressionDimred handles invalid input gracefully", {
    # Test with invalid method
    expect_error(plotGeneExpressionDimred(
        se_object = query_data,
        method = "InvalidMethod",
        pc_subset = 1:5,
        feature = "VPREB3"
    ))
    
    # Test with non-existent feature
    expect_error(plotGeneExpressionDimred(
        se_object = query_data,
        method = "PCA",
        pc_subset = 1:5,
        feature = "InvalidFeature"
    ))
})
