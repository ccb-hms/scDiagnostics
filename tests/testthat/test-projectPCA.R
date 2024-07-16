# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("projectPCA projects query data onto PCA space of reference data", {
    # Perform projection of query data onto reference PCA space
    pca_output <- projectPCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:10
    )

    # Check if output has expected number of rows
    expect_equal(nrow(pca_output), ncol(reference_data) + ncol(query_data))
    
    # Check if output has expected number of columns
    expect_equal(ncol(pca_output), length(1:10) + 2)  # +2 for dataset and cell_type columns
    
    # Check if dataset column has correct values
    expect_true(all(pca_output$dataset %in% c("Reference", "Query")))
    
    # Check if cell_type column has correct values
    expect_true(all(is.na(pca_output$cell_type) | pca_output$cell_type %in% unique(reference_data$expert_annotation)))
})

test_that("projectPCA handles incorrect parameters", {
    # Test for invalid query data column
    expect_error(projectPCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "invalid_column",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:10
    ))
    
    # Test for invalid reference data column
    expect_error(projectPCA(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "invalid_column",
        pc_subset = 1:10
    ))
})
