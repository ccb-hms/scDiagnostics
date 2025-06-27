# Load necessary libraries
library(testthat)
library(scDiagnostics)
library(scater)
library(SingleR)

# Load example datasets
data("reference_data")
data("query_data")

# Function to run all tests
test_that("calculateGraphIntegration", {

    # Create modified reference data (remove Myeloid)
    ref_data_mod <- reference_data[, reference_data$expert_annotation != "Myeloid"]
    ref_data_mod <- runPCA(ref_data_mod, ncomponents = 50)

    # Create SingleR annotations
    SingleR_annotation <- SingleR(query_data, ref_data_mod,
                                  labels = ref_data_mod$expert_annotation)
    query_data_mod <- query_data
    query_data_mod$SingleR_annotation <- SingleR_annotation$labels

    # Test 1: Basic functionality with default parameters
    result <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    # Check return object structure
    expect_s3_class(result, "calculateGraphIntegrationObject")
    expected_names <- c("high_query_prop_analysis", "cross_type_mixing",
                        "local_annotation_inconsistencies", "local_inconsistency_summary",
                        "community_composition", "annotation_consistency",
                        "overall_metrics", "graph_info", "parameters", "cell_info")
    expect_equal(names(result), expected_names)

    # Check data frame structures
    expect_true(is.data.frame(result$high_query_prop_analysis))
    expect_true(is.data.frame(result$cross_type_mixing))
    expect_true(is.data.frame(result$local_annotation_inconsistencies))
    expect_true(is.data.frame(result$community_composition))
    expect_true(is.data.frame(result$annotation_consistency))
    expect_true(is.data.frame(result$cell_info))

    # Check overall metrics structure
    expect_true(is.list(result$overall_metrics))
    expect_true(all(c("total_communities", "modularity", "mean_query_isolation_rate") %in%
                        names(result$overall_metrics)))

    # Check parameters are stored correctly
    expect_equal(result$parameters$pc_subset, 1:10)
    expect_equal(result$parameters$k_neighbors, 30)
    expect_equal(result$parameters$resolution, 0.1)

    # Test 2: Custom parameters
    result_custom <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:5,
        k_neighbors = 15,
        resolution = 0.2,
        high_query_prop_threshold = 0.8,
        cross_type_threshold = 0.15,
        local_consistency_threshold = 0.5,
        local_confidence_threshold = 0.3
    )

    expect_equal(result_custom$parameters$pc_subset, 1:5)
    expect_equal(result_custom$parameters$k_neighbors, 15)
    expect_equal(result_custom$parameters$resolution, 0.2)
    expect_equal(result_custom$parameters$high_query_prop_threshold, 0.8)

    # Test 3: Specific cell types
    specific_types <- c("CD4", "CD8", "B_and_plasma")
    result_specific <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = specific_types
    )

    expect_true(all(result_specific$parameters$cell_types_analyzed %in% specific_types))

    # Test 4: Different assay name
    # Add logcounts assay if it doesn't exist
    if (!"logcounts" %in% assayNames(query_data_mod)) {
        assay(query_data_mod, "logcounts") <- assay(query_data_mod, "counts")
    }
    if (!"logcounts" %in% assayNames(ref_data_mod)) {
        assay(ref_data_mod, "logcounts") <- assay(ref_data_mod, "counts")
    }

    result_logcounts <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        assay_name = "logcounts"
    )

    expect_s3_class(result_logcounts, "calculateGraphIntegrationObject")

    # Test 5: Error handling - invalid parameters

    # Invalid k_neighbors
    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            k_neighbors = -5
        ),
        "k_neighbors.*must be a positive integer"
    )

    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            k_neighbors = 1.5
        ),
        "k_neighbors.*must be a positive integer"
    )

    # Invalid resolution
    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            resolution = -0.1
        ),
        "resolution.*must be a positive number"
    )

    # Invalid high_query_prop_threshold
    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            high_query_prop_threshold = 0.3
        ),
        "high_query_prop_threshold.*must be between 0.5 and 1"
    )

    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            high_query_prop_threshold = 1.2
        ),
        "high_query_prop_threshold.*must be between 0.5 and 1"
    )

    # Invalid cross_type_threshold
    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            cross_type_threshold = 0.6
        ),
        "cross_type_threshold.*must be between 0 and 0.5"
    )

    # Invalid local_consistency_threshold
    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            local_consistency_threshold = 1.5
        ),
        "local_consistency_threshold.*must be between 0 and 1"
    )

    # Invalid local_confidence_threshold
    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            local_confidence_threshold = -0.1
        ),
        "local_confidence_threshold.*must be between 0 and 1"
    )

    # Test 6: Edge case - very small dataset
    # Create small subset
    small_query <- query_data_mod[, 1:50]
    small_ref <- ref_data_mod[, 1:50]
    small_ref <- runPCA(small_ref, ncomponents = 10)

    # Should handle small datasets gracefully
    result_small <- calculateGraphIntegration(
        query_data = small_query,
        reference_data = small_ref,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        min_cells_per_celltype = 5,
        min_cells_per_community = 3
    )

    expect_s3_class(result_small, "calculateGraphIntegrationObject")

    # Test 7: Edge case - k_neighbors larger than dataset
    # Should automatically adjust k_neighbors
    expect_warning(
        result_large_k <- calculateGraphIntegration(
            query_data = small_query,
            reference_data = small_ref,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            k_neighbors = 1000,
            min_cells_per_celltype = 5
        ),
        "k_neighbors reduced"
    )

    expect_s3_class(result_large_k, "calculateGraphIntegrationObject")

    # Test 8: Invalid cell type columns
    expect_error(
        calculateGraphIntegration(
            query_data = query_data_mod,
            reference_data = ref_data_mod,
            query_cell_type_col = "nonexistent_column",
            ref_cell_type_col = "expert_annotation"
        ),
        # This should trigger the argumentCheck function
        class = "error"
    )

    # Test 9: No common cell types after filtering
    # Create data with no overlapping cell types that meet minimum requirements
    query_small_types <- query_data_mod[, 1:10]  # Very small subset
    ref_small_types <- ref_data_mod[, 1:10]      # Very small subset

    expect_error(
        calculateGraphIntegration(
            query_data = query_small_types,
            reference_data = ref_small_types,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            min_cells_per_celltype = 50  # Very high threshold
        ),
        "No cell types meet the minimum cell count requirement"
    )

    # Test 10: Check numerical results make sense
    expect_true(result$overall_metrics$total_communities >= 0)
    expect_true(result$overall_metrics$modularity >= -1 && result$overall_metrics$modularity <= 1)
    expect_true(result$overall_metrics$mean_query_isolation_rate >= 0 &&
                    result$overall_metrics$mean_query_isolation_rate <= 1)
    expect_true(result$overall_metrics$mean_local_inconsistency_rate >= 0 &&
                    result$overall_metrics$mean_local_inconsistency_rate <= 1)

    # Test 11: Check community composition consistency
    if (nrow(result$community_composition) > 0) {
        expect_true(all(result$community_composition$total_cells >=
                            (result$community_composition$n_reference + result$community_composition$n_query)))
        expect_true(all(result$community_composition$query_proportion >= 0 &
                            result$community_composition$query_proportion <= 1))
    }

    # Test 12: Check annotation consistency metrics
    if (nrow(result$annotation_consistency) > 0) {
        expect_true(all(result$annotation_consistency$query_isolation_rate >= 0 &
                            result$annotation_consistency$query_isolation_rate <= 1))
        expect_true(all(result$annotation_consistency$true_cross_mixing_rate >= 0 &
                            result$annotation_consistency$true_cross_mixing_rate <= 1))
    }

    # Test 13: Check cell_info consistency
    expect_equal(nrow(result$cell_info),
                 sum(result$annotation_consistency$total_query_cells) +
                     sum(table(ref_data_mod$expert_annotation)[result$parameters$cell_types_analyzed]))

    expect_true(all(c("Query", "Reference") %in% result$cell_info$dataset))
    expect_true(all(result$cell_info$cell_type %in% result$parameters$cell_types_analyzed))

    # Test 14: Graph info structure
    expect_true(is.list(result$graph_info))
    expect_true(all(c("layout", "edges", "n_nodes") %in% names(result$graph_info)))
    expect_true(is.matrix(result$graph_info$layout))
    expect_equal(ncol(result$graph_info$layout), 2)  # 2D layout
    expect_equal(result$graph_info$n_nodes, nrow(result$cell_info))

    # Test 15: Local inconsistencies structure
    if (nrow(result$local_annotation_inconsistencies) > 0) {
        expect_true(all(result$local_annotation_inconsistencies$confidence_score >= 0))
        expect_true(all(result$local_annotation_inconsistencies$current_support_prop >= 0 &
                            result$local_annotation_inconsistencies$current_support_prop <= 1))
        expect_true(all(result$local_annotation_inconsistencies$suggested_support_prop >= 0 &
                            result$local_annotation_inconsistencies$suggested_support_prop <= 1))
    }

    # Test 16: Check that different thresholds produce different results
    result_strict <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        high_query_prop_threshold = 0.95,
        local_consistency_threshold = 0.8
    )

    result_lenient <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        high_query_prop_threshold = 0.7,
        local_consistency_threshold = 0.4
    )

    # Strict thresholds should generally detect fewer issues than lenient ones
    expect_true(nrow(result_strict$high_query_prop_analysis) <= nrow(result_lenient$high_query_prop_analysis))

    # Test 17: Different PC subsets
    result_fewer_pcs <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:5
    )

    result_more_pcs <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:20
    )

    expect_equal(result_fewer_pcs$parameters$pc_subset, 1:5)
    expect_equal(result_more_pcs$parameters$pc_subset, 1:20)

    # Both should produce valid results
    expect_s3_class(result_fewer_pcs, "calculateGraphIntegrationObject")
    expect_s3_class(result_more_pcs, "calculateGraphIntegrationObject")
})
