# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("query_data")

test_that("histQCvsAnnotation works correctly with default parameters", {
    # Generate histograms
    histograms <- histQCvsAnnotation(
        sce_object = query_data,
        cell_type_col = "SingleR_annotation",
        qc_col = "percent_mito",
        score_col = "annotation_scores"
    )

    # Check if the output is a ggplot object
    expect_true(inherits(histograms, "ggplot"))
})

test_that("histQCvsAnnotation works correctly with specific cell types", {
    # Generate histograms for specific cell types
    histograms <- histQCvsAnnotation(
        sce_object = query_data,
        cell_type_col = "SingleR_annotation",
        cell_types = c("CD4", "CD8"),
        qc_col = "percent_mito",
        score_col = "annotation_scores"
    )

    # Check if the output is a ggplot object
    expect_true(inherits(histograms, "ggplot"))
})

test_that("histQCvsAnnotation handles incorrect parameters", {
    expect_error(histQCvsAnnotation(
        sce_object = query_data,
        cell_type_col = "SingleR_annotation",
        qc_col = "invalid_column",
        score_col = "annotation_scores"
    ), "qc_col: 'invalid_column' is not a valid column name in sce_object.")

    expect_error(histQCvsAnnotation(
        sce_object = query_data,
        cell_type_col = "SingleR_annotation",
        qc_col = "percent_mito",
        score_col = "invalid_column"
    ), "score_col: 'invalid_column' is not a valid column name in sce_object.")
})

test_that("histQCvsAnnotation works with all cell types", {
    # Generate histograms for all cell types
    histograms <- histQCvsAnnotation(
        sce_object = query_data,
        cell_type_col = "SingleR_annotation",
        cell_types = NULL,
        qc_col = "percent_mito",
        score_col = "annotation_scores"
    )

    # Check if the output is a ggplot object
    expect_true(inherits(histograms, "ggplot"))
})
