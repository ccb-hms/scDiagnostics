# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("qc_data")

# Remove cell types with very few cells
qc_data_subset <- qc_data[, !(qc_data$SingleR_annotation
                              %in% c("Chondrocytes", "DC",
                                     "Neurons","Platelets"))]

test_that("plotQCvsAnnotation generates scatter plot correctly", {
    # Generate plot using the function
    p1 <- plotQCvsAnnotation(
        sce_object = qc_data_subset,
        cell_type_col = "SingleR_annotation",
        cell_types = NULL,
        qc_col = "total",
        score_col = "annotation_scores"
    )

    # Check if output is a ggplot object
    expect_s3_class(p1, "ggplot")

})

test_that("plotQCvsAnnotation handles invalid input gracefully", {
    # Test with non-existent column names in sce_object
    expect_error(plotQCvsAnnotation(
        sce_object = qc_data_subset,
        cell_type_col = "invalid_column",
        cell_types = NULL,
        qc_col = "total",
        score_col = "annotation_scores"
    ))

    # Test with non-existent qc_col
    expect_error(plotQCvsAnnotation(
        sce_object = qc_data_subset,
        cell_type_col = "SingleR_annotation",
        cell_types = NULL,
        qc_col = "invalid_column",
        score_col = "annotation_scores"
    ))

    # Test with non-existent score_col
    expect_error(plotQCvsAnnotation(
        sce_object = qc_data_subset,
        cell_type_col = "SingleR_annotation",
        cell_types = NULL,
        qc_col = "total",
        score_col = "invalid_column"
    ))
})
