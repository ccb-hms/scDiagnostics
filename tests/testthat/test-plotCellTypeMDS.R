# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("plotCellTypeMDS generates plots correctly", {
    # Generate plot using the function
    p1 <- plotCellTypeMDS(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid")[1:4]
    )
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")

})

test_that("plotCellTypeMDS handles invalid input gracefully", {
    # Test with invalid cell type column names
    expect_error(plotCellTypeMDS(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "InvalidColumn",
        ref_cell_type_col = "InvalidColumn",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid")[1:4]
    ))
    
    # Test with non-existent cell types
    expect_error(plotCellTypeMDS(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = "NonExistentCellType"
    ))
})
