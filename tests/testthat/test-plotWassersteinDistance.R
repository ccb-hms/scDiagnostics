# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Extract CD4 cells for subset analysis
ref_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]

# Selecting highly variable genes (can be customized by the user)
ref_top_genes <- scran::getTopHVGs(ref_data_subset, n = 500)
query_top_genes <- scran::getTopHVGs(query_data_subset, n = 500)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_top_genes, query_top_genes)
ref_data_subset <- ref_data_subset[common_genes, ]
query_data_subset <- query_data_subset[common_genes, ]

# Run PCA on reference data subset
ref_data_subset <- scater::runPCA(ref_data_subset)

test_that("plotWassersteinDistance computes and plots Wasserstein distances correctly", {
    
    # Generate plot using the function
    wasserstein_plot <- plotWassersteinDistance(
        query_data = query_data_subset,
        reference_data = ref_data_subset,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:5,
        n_resamples = 50,
        alpha = 0.05
    )
    
    # Check if output is a ggplot object
    expect_s3_class(wasserstein_plot, "ggplot")

    # Check if there are two vertical lines in the plot
    expect_equal(length(wasserstein_plot$layers[[2]]$data$xintercept), 2)
})

test_that("plotWassersteinDistance handles invalid input gracefully", {
    # Test with non-existent column names in query data
    expect_error(plotWassersteinDistance(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "invalid_column",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:5,
        n_resamples = 100,
        alpha = 0.05
    ))
    
    # Test with non-existent column names in reference data
    expect_error(plotWassersteinDistance(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "invalid_column",
        pc_subset = 1:5,
        n_resamples = 100,
        alpha = 0.05
    ))
    
    # Test with invalid alpha value
    expect_error(plotWassersteinDistance(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "expert_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:5,
        n_resamples = 100,
        alpha = 1.5
    ))
})
