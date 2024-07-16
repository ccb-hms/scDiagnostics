# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Test cases
test_that("boxplotPCA correctly handles typical inputs", {
    
    # Test with typical input
    plot <- boxplotPCA(query_data = query_data, 
                       reference_data = reference_data,
                       query_cell_type_col = "SingleR_annotation",
                       ref_cell_type_col = "expert_annotation",
                       cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                       pc_subset = 1:5)
    
    # Check if the result is a ggplot object
    expect_true("ggplot" %in% class(plot))
    
})

test_that("boxplotPCA handles NULL cell_types gracefully", {
    
    # Test with NULL cell_types
    plot <- boxplotPCA(query_data = query_data, 
                       reference_data = reference_data,
                       query_cell_type_col = "SingleR_annotation",
                       ref_cell_type_col = "expert_annotation",
                       pc_subset = 1:5)
    
    # Check if the result is a ggplot object
    expect_true("ggplot" %in% class(plot))
    
})

test_that("boxplotPCA handles invalid cell_types gracefully", {
    
    # Test with invalid cell_types
    expect_error(boxplotPCA(query_data = query_data, 
                            reference_data = reference_data,
                            query_cell_type_col = "SingleR_annotation",
                            ref_cell_type_col = "expert_annotation",
                            cell_types = c("InvalidCellType"),
                            pc_subset = 1:5))
    
})

test_that("boxplotPCA correctly handles single principal component subset", {
    
    # Test with single principal component subset
    plot <- boxplotPCA(query_data = query_data, 
                       reference_data = reference_data,
                       query_cell_type_col = "SingleR_annotation",
                       ref_cell_type_col = "expert_annotation",
                       cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                       pc_subset = 1)
    
    # Check if the result is a ggplot object
    expect_true("ggplot" %in% class(plot))

})
