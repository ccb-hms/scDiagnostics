# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load data
data("reference_data")
data("query_data")

# Store PCA anomaly data and plots
anomaly_output <- detectAnomaly(reference_data = reference_data,
                                query_data = query_data,
                                ref_cell_type_col = "expert_annotation",
                                query_cell_type_col = "SingleR_annotation",
                                pc_subset = 1:10,
                                n_tree = 500,
                                anomaly_threshold = 0.5)
top6_anomalies <- names(sort(anomaly_output$Combined$reference_anomaly_scores, decreasing = TRUE)[1:6])

# Compute cosine similarity between anomalies and top PCs
cosine_similarities <- calculateCellSimilarityPCA(sce_object = reference_data,
                                                  cell_names = top6_anomalies,
                                                  pc_subset = 1:25,
                                                  n_top_vars = 50)

test_that("plot.calculateCellSimilarityPCA generates plots correctly", {
    # Generate plot using the function
    p1 <- plot(cosine_similarities, pc_subset = 15:25)

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
})
