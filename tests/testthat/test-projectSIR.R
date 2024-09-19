# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

test_that("projectSIR returns expected output structure", {
    # Test that the function returns a list with the correct components
    result <- projectSIR(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    expect_type(result, "list")
    expect_named(result, c("cond_means", "rotation_mat", "sir_projections", "percent_var"))
})

test_that("projectSIR handles genes not found in query data", {
    # Modify query_data to exclude certain genes
    modified_query_data <- query_data[-1, ] # Remove the first gene

    expect_error(
        projectSIR(
            query_data = modified_query_data,
            reference_data = reference_data,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation"
        ),
        "Genes in reference SVD are not found in query data."
    )
})

test_that("projectSIR returns projections with correct number of components", {
    result <- projectSIR(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    expect_equal(ncol(result$sir_projections), length(result$percent_var) + 2)
})

