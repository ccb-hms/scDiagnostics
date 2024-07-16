# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Test cases
test_that("calculateCellDistances computes distances correctly", {
    
    # Test with typical input
    distance_data <- calculateCellDistances(query_data = query_data, 
                                            reference_data = reference_data,
                                            query_cell_type_col = "SingleR_annotation",
                                            ref_cell_type_col = "expert_annotation",
                                            cell_types = c("CD4", "CD8", "B_and_plasma"),
                                            pc_subset = 1:5)
    
    # Check if the result is a list
    expect_true(is.list(distance_data))
    
    # Check if the list contains entries for each cell type
    expect_equal(length(distance_data), length(c("CD4", "CD8", "B_and_plasma")))
    
    # Check structure of each entry
    expect_named(distance_data[["CD4"]], c("ref_distances", "query_to_ref_distances"))
    expect_named(distance_data[["CD8"]], c("ref_distances", "query_to_ref_distances"))
    expect_named(distance_data[["B_and_plasma"]], c("ref_distances", "query_to_ref_distances"))
    
    # Check dimensions and types of data within each entry if needed
    
})

test_that("calculateCellDistances handles NULL cell_types gracefully", {
    
    # Test with NULL cell_types
    distance_data <- calculateCellDistances(query_data = query_data, 
                                            reference_data = reference_data,
                                            query_cell_type_col = "SingleR_annotation",
                                            ref_cell_type_col = "expert_annotation",
                                            pc_subset = 1:5)
    
    # Check if the result is a list
    expect_true(is.list(distance_data))
    
    # Check if the list contains entries for each cell type present in the data
    expect_equal(length(distance_data), length(na.omit(unique(c(reference_data$expert_annotation, 
                                                                query_data$SingleR_annotation)))))
    
})

test_that("calculateCellDistances handles single cell type gracefully", {
    
    # Test with a single cell type
    distance_data <- calculateCellDistances(query_data = query_data, 
                                            reference_data = reference_data,
                                            query_cell_type_col = "SingleR_annotation",
                                            ref_cell_type_col = "expert_annotation",
                                            cell_types = "CD4",
                                            pc_subset = 1:5)
    
    # Check if the result is a list
    expect_true(is.list(distance_data))
    
    # Check if the list contains an entry for the specified cell type
    expect_named(distance_data, "CD4")
    
    # Check structure of the CD4 entry
    expect_named(distance_data[["CD4"]], c("ref_distances", "query_to_ref_distances"))
    
})

