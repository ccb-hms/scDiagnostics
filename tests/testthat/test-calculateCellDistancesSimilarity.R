# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Test cases
test_that("calculateCellDistancesSimilarity computes similarity measures correctly", {
    
    # Test with typical input
    overlap_measures <- calculateCellDistancesSimilarity(query_data = query_data, 
                                                         reference_data = reference_data,
                                                         query_cell_type_col = "SingleR_annotation",
                                                         ref_cell_type_col = "expert_annotation",
                                                         cell_names = c("AAACGGGCACATGACT-1", "ATTGGTGGTAACGCGA-1", "ATTGGACCATACTACG-1"),
                                                         pc_subset = 1:10)
    
    # Check if the result is a list
    expect_true(is.list(overlap_measures))
    
    # Check if the list contains entries for Bhattacharyya coefficients and Hellinger distances
    expect_named(overlap_measures, c("bhattacharyya_coef", "hellinger_dist"))

})

test_that("calculateCellDistancesSimilarity handles NULL cell_types gracefully", {
    
    # Test with NULL cell_types
    overlap_measures <- calculateCellDistancesSimilarity(query_data = query_data, 
                                                         reference_data = reference_data,
                                                         query_cell_type_col = "SingleR_annotation",
                                                         ref_cell_type_col = "expert_annotation",
                                                         cell_names = c("AAACGGGCACATGACT-1", "ATTGGTGGTAACGCGA-1", "ATTGGACCATACTACG-1"),
                                                         pc_subset = 1:10)
    
    # Check if the result is a list
    expect_true(is.list(overlap_measures))
    
    # Check if the list contains entries for Bhattacharyya coefficients and Hellinger distances
    expect_named(overlap_measures, c("bhattacharyya_coef", "hellinger_dist"))
    
    # Check structure of Bhattacharyya coefficients and Hellinger distances entries
    expect_equal(ncol(overlap_measures$bhattacharyya_coef), length(unique(c(reference_data$expert_annotation, query_data$SingleR_annotation))))
    expect_equal(ncol(overlap_measures$hellinger_dist), length(unique(c(reference_data$expert_annotation, query_data$SingleR_annotation))))
    
})

test_that("calculateCellDistancesSimilarity handles single cell names gracefully", {
    
    # Test with a single cell name
    overlap_measures <- calculateCellDistancesSimilarity(query_data = query_data, 
                                                         reference_data = reference_data,
                                                         query_cell_type_col = "SingleR_annotation",
                                                         ref_cell_type_col = "expert_annotation",
                                                         cell_names = "AAACGGGCACATGACT-1",
                                                         pc_subset = 1:10)
    
    # Check if the result is a list
    expect_true(is.list(overlap_measures))
    
    # Check if the list contains entries for Bhattacharyya coefficients and Hellinger distances
    expect_named(overlap_measures, c("bhattacharyya_coef", "hellinger_dist"))

})
