# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

## Unit tests for calculateCramerPValue function
test_that("calculateCramerPValue handles NULL query_data gracefully", {
    
    # Test with NULL query_data
    expect_error(calculateCramerPValue(
        reference_data = reference_data,
        query_data = NULL,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        cell_types = c("CD4", "CD8"),
        pc_subset = 1:5
    ))
    
})

test_that("calculateCramerPValue handles NULL query_cell_type_col gracefully", {
    
    # Test with NULL query_cell_type_col
    expect_error(calculateCramerPValue(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = NULL,
        cell_types = c("CD4", "CD8"),
        pc_subset = 1:5
    ))
    
})


test_that("calculateCramerPValue handles NULL cell_types gracefully", {
    
    # Test with NULL cell_types
    p_values <- calculateCramerPValue(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        cell_types = NULL,
        pc_subset = 1:5
    )
    
    # Check if the result is a named vector of p-values
    expect_type(p_values, "double")
    expect_true(length(p_values) > 0)  # Ensure some p-values are returned
    
    # Check the range of p-values
    expect_true(all(p_values >= 0 & p_values <= 1))
    
})

test_that("calculateCramerPValue returns correct p-values for specified cell types", {
    
    # Test with specific cell_types
    p_values <- calculateCramerPValue(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        cell_types = c("CD4", "CD8"),
        pc_subset = 1:5
    )
    
    # Check if the result is a named vector of p-values
    expect_type(p_values, "double")
    expect_true(all(names(p_values) %in% c("CD4", "CD8")))
    
    # Check the range of p-values
    expect_true(all(p_values >= 0 & p_values <= 1))
    
    # Add more specific checks as needed based on the expected behavior
    
})
