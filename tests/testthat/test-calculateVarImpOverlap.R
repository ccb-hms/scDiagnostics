# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Tests for the calculateVarImpOverlap function
test_that("calculateVarImpOverlap handles incorrect input arguments", {
    expect_error(calculateVarImpOverlap(reference_data = NULL, ref_cell_type_col = "expert_annotation"))
    expect_error(calculateVarImpOverlap(reference_data = reference_data, ref_cell_type_col = NULL))
    expect_error(calculateVarImpOverlap(reference_data = reference_data, ref_cell_type_col = "expert_annotation", n_tree = -1))
    expect_error(calculateVarImpOverlap(reference_data = reference_data, ref_cell_type_col = "expert_annotation", n_top = -1))
})

test_that("calculateVarImpOverlap calculates variable importance scores correctly for reference data", {
    result <- calculateVarImpOverlap(reference_data = reference_data, ref_cell_type_col = "expert_annotation", n_tree = 50, n_top = 10)
    expect_type(result, "list")
    expect_true("var_imp_ref" %in% names(result))
    expect_type(result$var_imp_ref, "list")
    for (cell_comb in names(result$var_imp_ref)) {
        expect_true(all(c("Gene", "RF_Importance") %in% colnames(result$var_imp_ref[[cell_comb]])))
    }
})

test_that("calculateVarImpOverlap calculates variable importance scores correctly for both reference and query data", {
    result <- calculateVarImpOverlap(reference_data = reference_data, query_data = query_data, ref_cell_type_col = "expert_annotation", 
                                     query_cell_type_col = "expert_annotation", n_tree = 50, n_top = 10)
    expect_type(result, "list")
    expect_true(all(c("var_imp_ref", "var_imp_query", "var_imp_comparison") %in% names(result)))
    expect_type(result$var_imp_ref, "list")
    expect_type(result$var_imp_query, "list")
    expect_type(result$var_imp_comparison, "double")
    for (cell_comb in names(result$var_imp_ref)) {
        expect_true(all(c("Gene", "RF_Importance") %in% colnames(result$var_imp_ref[[cell_comb]])))
        expect_true(all(c("Gene", "RF_Importance") %in% colnames(result$var_imp_query[[cell_comb]])))
    }
    expect_equal(length(result$var_imp_comparison), length(result$var_imp_ref))
})

test_that("calculateVarImpOverlap correctly compares top genes between reference and query data", {
    result <- calculateVarImpOverlap(reference_data = reference_data, query_data = query_data, ref_cell_type_col = "expert_annotation", 
                                     query_cell_type_col = "expert_annotation", n_tree = 50, n_top = 10)
    for (cell_comb in names(result$var_imp_ref)) {
        top_genes_ref <- result$var_imp_ref[[cell_comb]]$Gene[1:10]
        top_genes_query <- result$var_imp_query[[cell_comb]]$Gene[1:10]
        overlap <- length(intersect(top_genes_ref, top_genes_query)) / 10
        expect_equal(as.numeric(result$var_imp_comparison[cell_comb]), overlap)
    }
})
