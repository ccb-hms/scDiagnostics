# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("plotMarkerExpression generates density plots correctly", {
    # Generate plot using the function
    p1 <- plotMarkerExpression(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        gene_name = "VPREB3",
        cell_type = c("B_and_plasma")
    )
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    
})

test_that("plotMarkerExpression handles invalid input gracefully", {
    # Test with non-existent column names in query_data or reference_data
    expect_error(plotMarkerExpression(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "invalid_column",
        query_cell_type_col = "SingleR_annotation",
        gene_name = "VPREB3",
        cell_type = c("B_and_plasma")
    ))
    
    # Test with non-existent gene_name
    expect_error(plotMarkerExpression(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        gene_name = "invalid_gene_name",
        cell_type = c("B_and_plasma")
    ))
})

