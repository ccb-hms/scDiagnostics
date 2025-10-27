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

# Call the PCA comparison function
similarity_result <- comparePCA(reference_data = reference_data_subset,
                                query_data = query_data_subset,
                                query_cell_type_col = "expert_annotation",
                                ref_cell_type_col = "expert_annotation",
                                pc_subset = 1:5,
                                n_top_vars = 50,
                                metric = c("cosine", "correlation")[1],
                                correlation_method = c("spearman", "pearson")[1])

test_that("plot.comparePCAObject generates plots correctly", {
    # Generate plot using the function
    p1 <- plot(similarity_result)

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
})

test_that("plot.comparePCAObject works with different parameters", {
    # Test with show_values = FALSE
    p1 <- plot(similarity_result, show_values = FALSE)
    expect_s3_class(p1, "ggplot")

    # Test with show_significance = FALSE
    p2 <- plot(similarity_result, show_significance = FALSE)
    expect_s3_class(p2, "ggplot")

    # Test with custom color limits
    p3 <- plot(similarity_result, color_limits = c(-1, 1))
    expect_s3_class(p3, "ggplot")
})

test_that("plot.comparePCAObject works with correlation metric", {
    # Create result with correlation metric
    similarity_result_corr <- comparePCA(reference_data = reference_data_subset,
                                         query_data = query_data_subset,
                                         query_cell_type_col = "expert_annotation",
                                         ref_cell_type_col = "expert_annotation",
                                         pc_subset = 1:5,
                                         n_top_vars = 50,
                                         metric = "correlation",
                                         correlation_method = "spearman")

    # Generate plot
    p1 <- plot(similarity_result_corr)
    expect_s3_class(p1, "ggplot")
})
