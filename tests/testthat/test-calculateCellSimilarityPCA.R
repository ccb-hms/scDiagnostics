# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Test cases
test_that("calculateCellSimilarityPCA computes cosine similarity correctly", {

    # Test with typical input
    cosine_similarities <- calculateCellSimilarityPCA(sce_object = reference_data,
                                                      cell_names = colnames(reference_data)[1:5],
                                                      pc_subset = 1:5,
                                                      n_top_vars = 50)

    # Check if the result is a data frame
    expect_type(cosine_similarities, "double")

    # Check dimensions of the result
    expect_equal(nrow(cosine_similarities), length(colnames(reference_data)[1:5]))
    expect_equal(ncol(cosine_similarities), 5)  # Check against pc_subset

    # Check if column names correspond to PCs
    expect_true(all(grepl("^PC", colnames(cosine_similarities))))

    # Check range of cosine similarity values
    expect_true(all(cosine_similarities >= -1 & cosine_similarities <= 1))

})

test_that("calculateCellSimilarityPCA handles NULL n_top_vars gracefully.", {

    expect_error(cosine_similarities <-
                     calculateCellSimilarityPCA(sce_object = reference_data,
                                                cell_names = rownames(reference_data)[1:5],
                                                pc_subset = 1:5,
                                                n_top_vars = NULL))

})

