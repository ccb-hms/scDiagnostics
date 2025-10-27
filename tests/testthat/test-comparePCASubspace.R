# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Extract CD4 cells
reference_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]

# Selecting highly variable genes (can be customized by the user)
ref_top_genes <- scran::getTopHVGs(reference_data_subset, n = 500)
query_top_genes <- scran::getTopHVGs(query_data_subset, n = 500)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_top_genes, query_top_genes)
reference_data_subset <- reference_data_subset[common_genes,]
query_data_subset <- query_data_subset[common_genes,]

# Run PCA on datasets separately
reference_data_subset <- scater::runPCA(reference_data_subset)
query_data_subset <- scater::runPCA(query_data_subset)

# Unit tests for comparePCASubspace
test_that("comparePCASubspace returns expected output structure", {
    # Run the comparePCASubspace function
    subspace_comparison <- comparePCASubspace(query_data = query_data_subset,
                                              reference_data = reference_data_subset,
                                              query_cell_type_col = "expert_annotation",
                                              ref_cell_type_col = "expert_annotation",
                                              pc_subset = 1:5,
                                              n_top_vars = 50)

    # Check the class of the output
    expect_s3_class(subspace_comparison, "comparePCASubspaceObject")

    # Check the structure of the output list
    expect_true(is.list(subspace_comparison))
    expect_named(subspace_comparison, c("cosine_similarity", "cosine_id", "var_explained_ref", "var_explained_query", "var_explained_avg", "weighted_cosine_similarity"))
})

test_that("comparePCASubspace handles invalid n_top_vars argument", {
    # Check if function stops with invalid n_top_vars
    expect_error(comparePCASubspace(query_data = query_data_subset,
                                    reference_data = reference_data_subset,
                                    query_cell_type_col = "expert_annotation",
                                    ref_cell_type_col = "expert_annotation",
                                    pc_subset = 1:5,
                                    n_top_vars = -10))
})

test_that("comparePCASubspace handles missing common genes", {
    # Create a subset without common genes
    ref_subset_no_common <- reference_data[setdiff(rownames(reference_data), common_genes), ]
    query_subset_no_common <- query_data[setdiff(rownames(query_data), common_genes), ]

    # Check if function handles no common genes
    expect_error(comparePCASubspace(query_data = query_subset_no_common,
                                    reference_data = ref_subset_no_common,
                                    query_cell_type_col = "expert_annotation",
                                    ref_cell_type_col = "expert_annotation",
                                    pc_subset = 1:5,
                                    n_top_vars = 50))
})

test_that("comparePCASubspace calculates cosine similarities correctly", {
    # Run the comparePCASubspace function
    subspace_comparison <- comparePCASubspace(query_data = query_data_subset,
                                              reference_data = reference_data_subset,
                                              query_cell_type_col = "expert_annotation",
                                              ref_cell_type_col = "expert_annotation",
                                              pc_subset = 1:5,
                                              n_top_vars = 50)

    # Check that cosine_similarity contains numeric values
    expect_type(subspace_comparison$cosine_similarity, "double")
    # Check that cosine_similarity values are within range [0, 1]
    expect_true(all(subspace_comparison$cosine_similarity >= 0 & subspace_comparison$cosine_similarity <= 1))
})

test_that("comparePCASubspace returns correct dimensions for cosine_id", {
    # Run the comparePCASubspace function
    subspace_comparison <- comparePCASubspace(query_data = query_data_subset,
                                              reference_data = reference_data_subset,
                                              query_cell_type_col = "expert_annotation",
                                              ref_cell_type_col = "expert_annotation",
                                              pc_subset = 1:5,
                                              n_top_vars = 50)

    # Check the dimensions of the cosine_id matrix
    expect_equal(dim(subspace_comparison$cosine_id), c(5, 2))
})

test_that("comparePCASubspace returns variance explained components correctly", {
    # Run the comparePCASubspace function
    subspace_comparison <- comparePCASubspace(query_data = query_data_subset,
                                              reference_data = reference_data_subset,
                                              query_cell_type_col = "expert_annotation",
                                              ref_cell_type_col = "expert_annotation",
                                              pc_subset = 1:5,
                                              n_top_vars = 50)

    # Check that var_explained_ref is a numeric vector
    expect_type(subspace_comparison$var_explained_ref, "double")
    expect_equal(length(subspace_comparison$var_explained_ref), length(1:5))

    # Check that var_explained_query is a numeric vector
    expect_type(subspace_comparison$var_explained_query, "double")
    expect_equal(length(subspace_comparison$var_explained_query), length(1:5))

    # Check that var_explained_avg is a numeric vector
    expect_type(subspace_comparison$var_explained_avg, "double")
    expect_equal(length(subspace_comparison$var_explained_avg), length(1:5))
})

test_that("comparePCASubspace returns a valid weighted cosine similarity score", {
    # Run the comparePCASubspace function
    subspace_comparison <- comparePCASubspace(query_data = query_data_subset,
                                              reference_data = reference_data_subset,
                                              query_cell_type_col = "expert_annotation",
                                              ref_cell_type_col = "expert_annotation",
                                              pc_subset = 1:5,
                                              n_top_vars = 50)

    # Check that weighted_cosine_similarity is a numeric value
    expect_type(subspace_comparison$weighted_cosine_similarity, "double")
    # Check that weighted_cosine_similarity is within range [0, 1]
    expect_true(subspace_comparison$weighted_cosine_similarity >= 0 & subspace_comparison$weighted_cosine_similarity <= 1)
})
