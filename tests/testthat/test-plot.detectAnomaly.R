# Load necessary libraries
library(testthat)
library(scDiagnostics) 

# Load example datasets
data("reference_data")
data("query_data")

# Store PCA anomaly data
anomaly_output <- detectAnomaly(reference_data = reference_data, 
                                query_data = query_data, 
                                ref_cell_type_col = "expert_annotation", 
                                query_cell_type_col = "SingleR_annotation",
                                pc_subset = 1:10,
                                n_tree = 500,
                                anomaly_treshold = 0.5) 


test_that("plot.detectAnomaly generates plots correctly", {
    # Generate plot using the function (query data)
    p1 <- plot(anomaly_output, 
               cell_type = "CD4", 
               pc_subset = 1:5, 
               data_type = "query")
    
    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")
    
    # Generate plot using the function (reference data)
    p2 <- plot(anomaly_output, 
               cell_type = "CD4", 
               pc_subset = 1:5, 
               data_type = "reference")
    
    # Check if output is a ggplot object
    expect_s3_class(p2, "ggplot")
})
