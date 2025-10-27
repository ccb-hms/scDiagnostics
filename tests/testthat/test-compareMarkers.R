# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Test cases
test_that("compareMarkers correctly handles typical inputs", {

    # Basic functionality test
    result <- compareMarkers(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = "expert_annotation",
                             ref_cell_type_col = "expert_annotation")

    # Check that result is a list with correct class
    expect_true(is.list(result))
    expect_true("compareMarkersObject" %in% class(result))

    # Check all expected components are present
    expected_components <- c("marker_overlap", "expression_consistency", "quality_scores",
                             "markers_query", "markers_ref", "common_cell_types",
                             "n_cells_query", "n_cells_ref", "anomaly_filter_used",
                             "selected_cell_types", "anomaly_output")
    expect_true(all(expected_components %in% names(result)))

    # Check that marker_overlap is numeric and properly named
    expect_true(is.numeric(result$marker_overlap))
    expect_true(all(result$marker_overlap >= 0 & result$marker_overlap <= 1))
    expect_true(!is.null(names(result$marker_overlap)))

    # Check that expression_consistency is numeric and properly named
    expect_true(is.numeric(result$expression_consistency))
    expect_true(all(result$expression_consistency >= 0 & result$expression_consistency <= 1))
    expect_true(!is.null(names(result$expression_consistency)))

    # Check quality_scores are valid categories
    expect_true(all(result$quality_scores %in% c("Good", "Moderate", "Poor")))
    expect_true(!is.null(names(result$quality_scores)))

    # Check that markers are lists
    expect_true(is.list(result$markers_query))
    expect_true(is.list(result$markers_ref))

    # Check cell counts are numeric and positive
    expect_true(is.numeric(result$n_cells_query))
    expect_true(is.numeric(result$n_cells_ref))
    expect_true(all(result$n_cells_query > 0))
    expect_true(all(result$n_cells_ref > 0))

    # Check anomaly_filter_used defaults to "none"
    expect_equal(result$anomaly_filter_used, "none")
    expect_null(result$anomaly_output)
})

test_that("compareMarkers works with different parameter settings", {

    # Test with different n_markers
    result <- compareMarkers(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = "expert_annotation",
                             ref_cell_type_col = "expert_annotation",
                             n_markers = 20)

    expect_true("compareMarkersObject" %in% class(result))
    expect_true(length(result$common_cell_types) > 0)

    # Test with different min_cells
    result2 <- compareMarkers(query_data = query_data,
                              reference_data = reference_data,
                              query_cell_type_col = "expert_annotation",
                              ref_cell_type_col = "expert_annotation",
                              min_cells = 5)

    expect_true("compareMarkersObject" %in% class(result2))

    # Test with specific cell_types
    # First get available cell types
    query_types <- unique(colData(query_data)[["expert_annotation"]])
    ref_types <- unique(colData(reference_data)[["expert_annotation"]])
    common_types <- intersect(query_types, ref_types)

    if (length(common_types) > 1) {
        selected_types <- common_types[1:min(2, length(common_types))]
        result3 <- compareMarkers(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = "expert_annotation",
                                  ref_cell_type_col = "expert_annotation",
                                  cell_types = selected_types)

        expect_true("compareMarkersObject" %in% class(result3))
        expect_true(all(result3$selected_cell_types %in% selected_types))
    }
})

test_that("compareMarkers handles edge cases appropriately", {

    # Test error when no common cell types exist
    # Create a modified query with different cell type names
    query_modified <- query_data
    colData(query_modified)[["different_annotation"]] <- paste0("type_",
                                                                colData(query_data)[["expert_annotation"]])

    expect_error(
        compareMarkers(query_data = query_modified,
                       reference_data = reference_data,
                       query_cell_type_col = "different_annotation",
                       ref_cell_type_col = "expert_annotation"),
        "No common cell types"
    )

    # Test with very high min_cells requirement
    expect_error(
        compareMarkers(query_data = query_data,
                       reference_data = reference_data,
                       query_cell_type_col = "expert_annotation",
                       ref_cell_type_col = "expert_annotation",
                       min_cells = 10000),
        "No common cell types"
    )
})

test_that("compareMarkers anomaly filtering works correctly", {

    # Test that anomaly_filter argument is properly handled
    result_none <- compareMarkers(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = "expert_annotation",
                                  ref_cell_type_col = "expert_annotation",
                                  anomaly_filter = "none")

    expect_equal(result_none$anomaly_filter_used, "none")
    expect_null(result_none$anomaly_output)

    # Note: Testing actual anomaly filtering would require the detectAnomaly function
    # to work properly, so we mainly test that the parameter is handled correctly
})

test_that("compareMarkers produces consistent marker gene results", {

    result <- compareMarkers(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = "expert_annotation",
                             ref_cell_type_col = "expert_annotation")

    # Check that each cell type has corresponding markers
    for (cell_type in result$common_cell_types) {
        expect_true(cell_type %in% names(result$markers_query))
        expect_true(cell_type %in% names(result$markers_ref))

        # Check marker data frame structure
        query_markers <- result$markers_query[[cell_type]]
        ref_markers <- result$markers_ref[[cell_type]]

        if (nrow(query_markers) > 0) {
            expect_true(all(c("gene", "pval", "adj_pval", "logFC") %in% colnames(query_markers)))
            expect_true(all(query_markers$logFC > 0))
            expect_true(all(query_markers$adj_pval < 0.05))
        }

        if (nrow(ref_markers) > 0) {
            expect_true(all(c("gene", "pval", "adj_pval", "logFC") %in% colnames(ref_markers)))
            expect_true(all(ref_markers$logFC > 0))
            expect_true(all(ref_markers$adj_pval < 0.05))
        }
    }
})

test_that("compareMarkers handles different assay names", {

    # Test with default assay
    result <- compareMarkers(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = "expert_annotation",
                             ref_cell_type_col = "expert_annotation",
                             assay_name = "logcounts")

    expect_true("compareMarkersObject" %in% class(result))

    # Test error with non-existent assay
    expect_error(
        compareMarkers(query_data = query_data,
                       reference_data = reference_data,
                       query_cell_type_col = "expert_annotation",
                       ref_cell_type_col = "expert_annotation",
                       assay_name = "nonexistent_assay")
    )
})

test_that("compareMarkers argument validation works", {

    # Test invalid anomaly_filter argument
    expect_error(
        compareMarkers(query_data = query_data,
                       reference_data = reference_data,
                       query_cell_type_col = "expert_annotation",
                       ref_cell_type_col = "expert_annotation",
                       anomaly_filter = "invalid_option")
    )

    # Test with non-existent cell type column
    expect_error(
        compareMarkers(query_data = query_data,
                       reference_data = reference_data,
                       query_cell_type_col = "nonexistent_column",
                       ref_cell_type_col = "expert_annotation")
    )
})
