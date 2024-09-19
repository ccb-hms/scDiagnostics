# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Function to run all tests
test_that("calculateNearestNeighborProbabilities function tests", {

    # Test case 1: Check if function returns a list
    nn_output <- calculateNearestNeighborProbabilities(query_data = query_data,
                                                       reference_data = reference_data,
                                                       query_cell_type_col = "SingleR_annotation",
                                                       ref_cell_type_col = "expert_annotation")
    expect_type(nn_output, "list")

    # Test case 2: Check if function returns correct class
    expect_true(inherits(nn_output, "calculateNearestNeighborProbabilitiesObject"))

    # Test case 3: Check for non-SingleCellExperiment input
    expect_error(calculateNearestNeighborProbabilities(query_data = matrix(runif(1000), nrow = 100),
                                                       reference_data = reference_data,
                                                       query_cell_type_col = "SingleR_annotation",
                                                       ref_cell_type_col = "expert_annotation"))
    expect_error(calculateNearestNeighborProbabilities(query_data = query_data,
                                                       reference_data = matrix(runif(1000), nrow = 100),
                                                       query_cell_type_col = "SingleR_annotation",
                                                       ref_cell_type_col = "expert_annotation"))

    # Test case 4: Check for incorrect column names
    expect_error(calculateNearestNeighborProbabilities(query_data = query_data,
                                                       reference_data = reference_data,
                                                       query_cell_type_col = "Nonexistent_column",
                                                       ref_cell_type_col = "expert_annotation"))
    expect_error(calculateNearestNeighborProbabilities(query_data = query_data,
                                                       reference_data = reference_data,
                                                       query_cell_type_col = "SingleR_annotation",
                                                       ref_cell_type_col = "Nonexistent_column"))

    # Test case 5: Check for invalid n_neighbor parameter
    expect_error(calculateNearestNeighborProbabilities(query_data = query_data,
                                                       reference_data = reference_data,
                                                       query_cell_type_col = "SingleR_annotation",
                                                       ref_cell_type_col = "expert_annotation",
                                                       n_neighbor = -5),
                 "'n_neighbor' must be a positive integer.")
    expect_error(calculateNearestNeighborProbabilities(query_data = query_data,
                                                       reference_data = reference_data,
                                                       query_cell_type_col = "SingleR_annotation",
                                                       ref_cell_type_col = "expert_annotation",
                                                       n_neighbor = "twenty"))

    # Test case 6: Check for empty cell types
    nn_output <- calculateNearestNeighborProbabilities(query_data = query_data,
                                                       reference_data = reference_data,
                                                       query_cell_type_col = "SingleR_annotation",
                                                       ref_cell_type_col = "expert_annotation",
                                                       cell_types = character(0))
    expect_true(length(nn_output) == 0)

    # Test case 7: Check for correct computation of probabilities
    nn_output <- calculateNearestNeighborProbabilities(query_data = query_data,
                                                       reference_data = reference_data,
                                                       query_cell_type_col = "SingleR_annotation",
                                                       ref_cell_type_col = "expert_annotation",
                                                       cell_types = "CD4",
                                                       pc_subset = 1:5,
                                                       n_neighbor = 5)
    expect_true(all(names(nn_output) == "CD4"))
    expect_type(nn_output$CD4$query_prob, "double")

})
