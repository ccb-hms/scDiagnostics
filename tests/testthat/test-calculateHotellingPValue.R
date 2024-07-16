# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Unit tests for calculateHotellingPValue function
test_that("calculateHotellingPValue function tests", {
    
    # Test case 1: Check if function returns a named numeric vector
    p_values <- calculateHotellingPValue(query_data = query_data,
                                         reference_data = reference_data,
                                         query_cell_type_col = "SingleR_annotation",
                                         ref_cell_type_col = "expert_annotation",
                                         pc_subset = 1:5,
                                         n_permutation = 100)
    
    expect_type(p_values, "double")
    expect_true(length(p_values) > 0)

    # Test case 2: Check if p-values are between 0 and 1
    expect_true(all(p_values >= 0 & p_values <= 1))
    
})