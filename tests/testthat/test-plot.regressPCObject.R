# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Plot the PC data (no query data)
regress_res_ref <- regressPC(reference_data = reference_data,
                             ref_cell_type_col = "expert_annotation", 
                             cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                             pc_subset = 1:25)

# Plot the PC data (with query data)
regress_res_query <- regressPC(reference_data = reference_data,
                               query_data = query_data,
                               ref_cell_type_col = "expert_annotation", 
                               query_cell_type_col = "expert_annotation",
                               cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
                               pc_subset = 1:25)

test_that("plot.regressPC generates plots correctly", {
    # Generate plot using the function (no query data)
    p1 <- plot(regress_res_ref, plot_type = "r_squared")
    p2 <- plot(regress_res_ref, plot_type = "p-value", alpha = 0.05)
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    
    # Generate plot using the function (with query data)
    p3 <- plot(regress_res_query, plot_type = "r_squared")
    p4 <- plot(regress_res_query, plot_type = "p-value")
    
    # Check if output is a ggplot object
    expect_s3_class(p3, "ggplot")
    expect_s3_class(p4, "ggplot")
})
