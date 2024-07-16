# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("plotPairwiseDistancesDensity generates density plot correctly", {
    # Generate plot using the function
    p1 <- plotPairwiseDistancesDensity(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_type_query = "CD8",
        cell_type_ref = "CD8",
        pc_subset = 1:5,
        distance_metric = "euclidean",
        correlation_method = "pearson"
    )
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    
})

test_that("plotPairwiseDistancesDensity handles invalid input gracefully", {
    # Test with non-existent column names in query_data or reference_data
    expect_error(plotPairwiseDistancesDensity(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "invalid_column",
        ref_cell_type_col = "expert_annotation",
        cell_type_query = "CD8",
        cell_type_ref = "CD8",
        pc_subset = 1:5,
        distance_metric = "euclidean",
        correlation_method = "pearson"
    ))
    
    # Test with non-existent cell types
    expect_error(plotPairwiseDistancesDensity(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_type_query = "invalid_cell_type",
        cell_type_ref = "CD8",
        pc_subset = 1:5,
        distance_metric = "euclidean",
        correlation_method = "pearson"
    ))
    
    # Test with invalid distance_metric
    expect_error(plotPairwiseDistancesDensity(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_type_query = "CD8",
        cell_type_ref = "CD8",
        pc_subset = 1:5,
        distance_metric = "invalid_metric",
        correlation_method = "pearson"
    ))
    
    # Test with invalid correlation_method
    expect_error(plotPairwiseDistancesDensity(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_type_query = "CD8",
        cell_type_ref = "CD8",
        pc_subset = 1:5,
        distance_metric = "correlation",
        correlation_method = "invalid_method"
    ))
})
