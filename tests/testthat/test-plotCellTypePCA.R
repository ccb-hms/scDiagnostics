# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("plotCellTypePCA generates plots correctly", {
    # Generate plot using the function
    p1 <- plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    )
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    
})

test_that("plotCellTypePCA handles invalid input gracefully", {
    # Test with invalid column names for cell types
    expect_error(plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "InvalidColumn",
        ref_cell_type_col = "InvalidColumn",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    ))
    
    # Test with non-existent cell types
    expect_error(plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = "NonExistentCellType",
        pc_subset = 1:5
    ))
    
    # Test with invalid principal component subset
    expect_error(plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = c(10, 20, 30)
    ))
})
