# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load data
data("reference_data")
data("query_data")

# Plot the PC data
distance_data <- calculateCellDistances(query_data = query_data,
                                        reference_data = reference_data,
                                        query_cell_type_col = "SingleR_annotation",
                                        ref_cell_type_col = "expert_annotation",
                                        pc_subset = 1:10)

# Identify outliers for CD4
cd4_anomalies <- detectAnomaly(reference_data = reference_data,
                               query_data = query_data,
                               query_cell_type_col = "SingleR_annotation",
                               ref_cell_type_col = "expert_annotation",
                               pc_subset = 1:10,
                               n_tree = 500,
                               anomaly_threshold = 0.5)$CD4
cd4_top6_anomalies <- names(sort(cd4_anomalies$query_anomaly_scores, decreasing = TRUE)[1:6])

test_that("plot.calculateCellDistances generates plots correctly", {
    # Generate plot using the function
    p1 <- plot(distance_data, ref_cell_type = "CD4", cell_names = cd4_top6_anomalies)
    p2 <- plot(distance_data, ref_cell_type = "CD8", cell_names = cd4_top6_anomalies)

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})
