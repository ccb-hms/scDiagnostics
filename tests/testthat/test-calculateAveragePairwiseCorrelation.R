# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Test cases
test_that("calculateAveragePairwiseCorrelation computes average pairwise correlations", {
    
    # Test with typical input
    cor_matrix_avg <- calculateAveragePairwiseCorrelation(query_data = query_data, 
                                                          reference_data = reference_data,
                                                          query_cell_type_col = "SingleR_annotation",
                                                          ref_cell_type_col = "expert_annotation",
                                                          cell_types = c("CD4", "CD8", "B_and_plasma"),
                                                          pc_subset = 1:10,
                                                          correlation_method = "pearson")
    
    # Check if the result is a matrix
    expect_true(is.matrix(cor_matrix_avg))
    
    # Check dimensions of the output matrix
    expect_equal(dim(cor_matrix_avg), c(length(c("Query-CD4", "Query-CD8", "Query-B_and_plasma")), 
                                        length(c("Ref-CD4", "Ref-CD8", "Ref-B_and_plasma"))))
    
    # Check specific values if needed
    
})

test_that("calculateAveragePairwiseCorrelation handles NULL cell_types gracefully", {
    
    # Test with NULL cell_types
    cor_matrix_avg <- calculateAveragePairwiseCorrelation(query_data = query_data, 
                                                          reference_data = reference_data,
                                                          query_cell_type_col = "SingleR_annotation",
                                                          ref_cell_type_col = "expert_annotation",
                                                          pc_subset = 1:10,
                                                          correlation_method = "spearman")
    
    # Check if the result is a matrix
    expect_true(is.matrix(cor_matrix_avg))
    
    # Check dimensions of the output matrix
    expect_equal(dim(cor_matrix_avg)[1], dim(cor_matrix_avg)[2])  # Square matrix
    
})

test_that("calculateAveragePairwiseCorrelation handles single cell type gracefully", {
    
    # Test with a single cell type
    cor_matrix_avg <- calculateAveragePairwiseCorrelation(query_data = query_data, 
                                                          reference_data = reference_data,
                                                          query_cell_type_col = "SingleR_annotation",
                                                          ref_cell_type_col = "expert_annotation",
                                                          cell_types = "CD4",
                                                          pc_subset = 1:10,
                                                          correlation_method = "pearson")
    
    # Check if the result is a matrix
    expect_true(is.matrix(cor_matrix_avg))
    
    # Check dimensions of the output matrix
    expect_equal(dim(cor_matrix_avg), c(1, 1))
    
})

test_that("calculateAveragePairwiseCorrelation handles invalid correlation methods", {
    
    # Test with an invalid correlation method
    expect_error(calculateAveragePairwiseCorrelation(query_data = query_data, 
                                                     reference_data = reference_data,
                                                     query_cell_type_col = "SingleR_annotation",
                                                     ref_cell_type_col = "expert_annotation",
                                                     cell_types = c("CD4", "CD8", "B_and_plasma"),
                                                     pc_subset = 1:10,
                                                     correlation_method = "invalid_method"))
    
})
