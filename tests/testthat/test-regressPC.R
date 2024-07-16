# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("regressPC works correctly with only reference data", {
    # Perform regression with only reference data
    regress_res <- regressPC(
        reference_data = reference_data,
        ref_cell_type_col = "expert_annotation",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    )
    
    # Check if the output is a list
    expect_true(is.list(regress_res))
    
    # Check if the list contains the correct elements
    expect_true("regression_summaries" %in% names(regress_res))
    expect_true("r_squared" %in% names(regress_res))
    expect_true("var_contributions" %in% names(regress_res))
    expect_true("total_variance_explained" %in% names(regress_res))
    
    # Check if R-squared values are numeric and of correct length
    expect_true(is.numeric(regress_res$r_squared))
    expect_equal(length(regress_res$r_squared), 5)
    
    # Check if total variance explained is a numeric value
    expect_true(is.numeric(regress_res$total_variance_explained))
})

test_that("regressPC works correctly with reference and query data", {
    # Perform regression with reference and query data
    regress_res <- regressPC(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "SingleR_annotation",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    )
    
    # Check if the output is a list
    expect_true(is.list(regress_res))
    
    # Check if the list contains the correct elements
    for (cell_type in c("CD4", "CD8", "B_and_plasma", "Myeloid")) {
        expect_true(cell_type %in% names(regress_res))
    }
    
    expect_true("indep_var" %in% names(regress_res))
    expect_equal(regress_res$indep_var, "dataset")
    
    # Check if each cell type has regression summaries
    for (cell_type in c("CD4", "CD8", "B_and_plasma", "Myeloid")) {
        expect_true(is.list(regress_res[[cell_type]]))
        expect_equal(length(regress_res[[cell_type]]), 5)
    }
})

test_that("regressPC handles incorrect parameters", {
    # Test for invalid reference data column
    expect_error(regressPC(
        reference_data = reference_data,
        ref_cell_type_col = "invalid_column",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    ))
    
    # Test for invalid query data column
    expect_error(regressPC(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = "expert_annotation",
        query_cell_type_col = "invalid_column",
        cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
        pc_subset = 1:5
    ))
    
    # Test for invalid cell types
    expect_error(regressPC(
        reference_data = reference_data,
        ref_cell_type_col = "expert_annotation",
        cell_types = c("Invalid_Cell_Type"),
        pc_subset = 1:5
    ))
})

test_that("regressPC works with all cell types and default PC subset", {
    # Perform regression with all cell types and default PC subset
    regress_res <- regressPC(
        reference_data = reference_data,
        ref_cell_type_col = "expert_annotation"
    )
    
    # Check if the output is a list
    expect_true(is.list(regress_res))
    
    # Check if the list contains the correct elements
    expect_true("regression_summaries" %in% names(regress_res))
    expect_true("r_squared" %in% names(regress_res))
    expect_true("var_contributions" %in% names(regress_res))
    expect_true("total_variance_explained" %in% names(regress_res))
    
    # Check if R-squared values are numeric and of correct length
    expect_true(is.numeric(regress_res$r_squared))
    expect_equal(length(regress_res$r_squared), 10)
    
    # Check if total variance explained is a numeric value
    expect_true(is.numeric(regress_res$total_variance_explained))
})

