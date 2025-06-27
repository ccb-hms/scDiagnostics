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

    # Check if output is a ggmatrix object (from GGally::ggpairs)
    expect_s3_class(p1, "ggmatrix")
})

test_that("plotCellTypePCA works with different facet options", {
    # Test different facet combinations
    p1 <- plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8"),
        pc_subset = 1:3,
        lower_facet = "scatter",
        diagonal_facet = "density",
        upper_facet = "blank"
    )

    expect_s3_class(p1, "ggmatrix")

    # Test with different options
    p2 <- plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8"),
        pc_subset = 1:3,
        lower_facet = "contour",
        diagonal_facet = "boxplot",
        upper_facet = "ellipse"
    )

    expect_s3_class(p2, "ggmatrix")
})

test_that("plotCellTypePCA works with default parameters", {
    # Test with minimal parameters (should use defaults)
    p1 <- plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    expect_s3_class(p1, "ggmatrix")
})

test_that("plotCellTypePCA handles invalid input gracefully", {
    # Test with invalid column names for cell types
    expect_error(plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "InvalidColumn",
        ref_cell_type_col = "InvalidColumn",
        cell_types = c("CD4", "CD8"),
        pc_subset = 1:3
    ))

    # Test with invalid principal component subset (too high)
    expect_error(plotCellTypePCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8"),
        pc_subset = c(100, 200, 300)
    ))
})
