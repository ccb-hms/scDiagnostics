# Load necessary libraries
library(testthat)
library(scDiagnostics)
library(ggplot2)

# Load example datasets and prepare data as before
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

# Generate Wasserstein distance data
wasserstein_data <- calculateWassersteinDistance(
    query_data = query_data_subset,
    reference_data = ref_data_subset,
    query_cell_type_col = "expert_annotation",
    ref_cell_type_col = "expert_annotation",
    pc_subset = 1:5,
    n_resamples = 100
)

# Test to ensure that the plot function returns a ggplot object
test_that("plot method for calculateWassersteinDistanceObject returns a ggplot object", {
    plot_obj <- plot.calculateWassersteinDistanceObject(wasserstein_data)
    expect_true(inherits(plot_obj, "ggplot"))
})

# Test to ensure plot structure contains expected layers and aesthetics
test_that("plot structure and aesthetics are correct", {
    plot_obj <- plot.calculateWassersteinDistanceObject(wasserstein_data)

    # Check that the plot contains a density layer
    expect_true(any(sapply(plot_obj$layers, function(layer) inherits(layer$geom, "GeomDensity"))))

    # Check for vertical lines
    expect_true(any(sapply(plot_obj$layers, function(layer) inherits(layer$geom, "GeomVline"))))

    # Verify labels are correct
    expect_equal(plot_obj$labels$title, paste0("Density of Wasserstein Distances For Reference Distribution of ", wasserstein_data$cell_type))
    expect_equal(plot_obj$labels$x, "Wasserstein Distances")
    expect_equal(plot_obj$labels$y, "Density")
})

# Test to ensure that alpha level alters the plot appropriately
test_that("alpha level impacts the plot threshold line", {
    alpha_test <- 0.01
    plot_obj <- plot.calculateWassersteinDistanceObject(wasserstein_data, alpha = alpha_test)

    # Extract the x-intercepts of the significance threshold line
    vline_data <- ggplot_build(plot_obj)$data[[2]]

    # Check that first vline matches quantile of null_dist based on alpha level
    expect_equal(vline_data$xintercept[1], as.numeric(quantile(wasserstein_data$null_dist, 1 - alpha_test)), tolerance = 0.01)
})
