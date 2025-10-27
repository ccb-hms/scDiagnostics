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

# Compare PCA subspaces
subspace_comparison <- comparePCASubspace(reference_data = reference_data_subset,
                                          query_data = query_data_subset,
                                          query_cell_type_col = "expert_annotation",
                                          ref_cell_type_col = "expert_annotation",
                                          n_top_vars = 50,
                                          pc_subset = 1:5)

test_that("plot.comparePCASubspaceObject generates plots correctly", {
    # Generate plot using the function
    p1 <- plot(subspace_comparison)

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
})

test_that("plot.comparePCASubspaceObject creates expected plot elements", {
    # Generate plot using the function
    p1 <- plot(subspace_comparison)

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")

    # Check that the plot has the expected layers
    expect_true(length(p1$layers) > 0)

    # Check that the plot has appropriate labels
    expect_true(!is.null(p1$labels$title))
    expect_true(!is.null(p1$labels$subtitle))
    expect_true(!is.null(p1$labels$x))
    expect_true(!is.null(p1$labels$y))
})
