# Load necessary libraries
library(testthat)
library(scDiagnostics)
library(ggplot2)
library(scater)
library(SingleR)

# Load example datasets
data("reference_data")
data("query_data")

test_that("plot.calculateGraphIntegrationObject", {

    # Create modified reference data (remove Myeloid)
    ref_data_mod <- reference_data[, reference_data$expert_annotation != "Myeloid"]
    ref_data_mod <- runPCA(ref_data_mod, ncomponents = 50)

    # Create SingleR annotations
    SingleR_annotation <- SingleR(query_data, ref_data_mod,
                                  labels = ref_data_mod$expert_annotation)
    query_data_mod <- query_data
    query_data_mod$SingleR_annotation <- SingleR_annotation$labels

    # Create diagnostic results for testing
    diag_result <- calculateGraphIntegration(
        query_data = query_data_mod,
        reference_data = ref_data_mod,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    # Test 1: Default plot (community_network)
    p1 <- plot(diag_result)
    expect_s3_class(p1, "ggplot")
    expect_true(grepl("Community Network", p1$labels$title))

    # Test 2: Community network plot with cell_type coloring
    p2 <- plot(diag_result, plot_type = "community_network", color_by = "cell_type")
    expect_s3_class(p2, "ggplot")
    expect_equal(p2$labels$title, "Community Network: Cell Type Composition")

    # Test 3: Community network plot with community_type coloring
    p3 <- plot(diag_result, plot_type = "community_network", color_by = "community_type")
    expect_s3_class(p3, "ggplot")
    expect_equal(p3$labels$title, "Community Network: Inter-Community Connection Strength")

    # Test 4: Cell network plot with cell_type coloring
    p4 <- plot(diag_result, plot_type = "cell_network", color_by = "cell_type")
    expect_s3_class(p4, "ggplot")
    expect_equal(p4$labels$title, "Cell Network: Cell Type Distribution")

    # Test 5: Cell network plot with community_type coloring
    p5 <- plot(diag_result, plot_type = "cell_network", color_by = "community_type")
    expect_s3_class(p5, "ggplot")
    expect_equal(p5$labels$title, "Cell Network: Annotation Consistency Diagnostics")

    # Test 6: Community data plot
    p6 <- plot(diag_result, plot_type = "community_data")
    expect_s3_class(p6, "ggplot")
    expect_equal(p6$labels$title, "Community Data Overview")
    expect_equal(p6$labels$x, "Query Proportion")
    expect_equal(p6$labels$y, "Total Cells in Community")

    # Test 7: Summary plot
    p7 <- plot(diag_result, plot_type = "summary")
    expect_s3_class(p7, "ggplot")
    expect_equal(p7$labels$title, "Total Annotation Issues by Cell Type")
    expect_true("GeomCol" %in% class(p7$layers[[1]]$geom))

    # Test 8: Local issues plot
    p8 <- plot(diag_result, plot_type = "local_issues")
    expect_s3_class(p8, "ggplot")
    expect_true(grepl("Local", p8$labels$title))

    # Test 9: Annotation issues plot
    p9 <- plot(diag_result, plot_type = "annotation_issues")
    expect_s3_class(p9, "ggplot")
    expect_true(grepl("Annotation", p9$labels$title))

    # Test 10: Invalid plot_type
    expect_error(
        plot(diag_result, plot_type = "invalid_type"),
        "'arg' should be one of"
    )

    # Test 11: Invalid color_by
    expect_error(
        plot(diag_result, plot_type = "community_network", color_by = "invalid_color"),
        "'arg' should be one of"
    )

    # Test 12: Custom parameters
    p12 <- plot(diag_result, plot_type = "cell_network",
                max_nodes = 500, point_size = 1.5)
    expect_s3_class(p12, "ggplot")

    # Test 13: exclude_reference_only parameter
    p13 <- plot(diag_result, plot_type = "community_data",
                exclude_reference_only = TRUE)
    expect_s3_class(p13, "ggplot")

    p14 <- plot(diag_result, plot_type = "cell_network",
                exclude_reference_only = TRUE)
    expect_s3_class(p14, "ggplot")

    # Test 14: Test color schemes are applied correctly
    p_colors <- plot(diag_result, plot_type = "community_data")
    expect_true("ScaleDiscrete" %in% class(p_colors$scales$get_scales("colour")))

    # Test 15: Test with large dataset (triggering sampling)
    p_large <- plot(diag_result, plot_type = "cell_network", max_nodes = 100)
    expect_s3_class(p_large, "ggplot")

    # Test 16: Test subtitle generation
    p_subtitle <- plot(diag_result, plot_type = "cell_network")
    expect_true(!is.null(p_subtitle$labels$subtitle))

    # Test 17: Test coordinate systems
    p_coord <- plot(diag_result, plot_type = "summary")
    expect_true("CoordCartesian" %in% class(p_coord$coordinates))

    # Test 18: Test themes
    p_void <- plot(diag_result, plot_type = "cell_network")
    expect_true("theme" %in% names(p_void))

    # Test 19: Test scale transformations
    p_scales <- plot(diag_result, plot_type = "community_data")
    size_scale <- p_scales$scales$get_scales("size")
    expect_true("ScaleContinuous" %in% class(size_scale))

    # Test 20: Test with exclude_reference_only filtering
    p_with_ref <- plot(diag_result, plot_type = "community_network",
                       exclude_reference_only = FALSE)
    p_without_ref <- plot(diag_result, plot_type = "community_network",
                          exclude_reference_only = TRUE)
    expect_s3_class(p_with_ref, "ggplot")
    expect_s3_class(p_without_ref, "ggplot")

    # Test 21: Test legend positioning
    p_legend <- plot(diag_result, plot_type = "community_data")
    expect_equal(p_legend$theme$legend.position, "right")

    p_legend_bottom <- plot(diag_result, plot_type = "summary")
    expect_equal(p_legend_bottom$theme$legend.position, "bottom")

    # Test 22: Test coordinate flipping in annotation_issues plot
    if (nrow(diag_result$high_query_prop_analysis) > 0 ||
        nrow(diag_result$cross_type_mixing) > 0 ||
        nrow(diag_result$local_annotation_inconsistencies) > 0) {

        p_flip <- plot(diag_result, plot_type = "annotation_issues")
        expect_true("CoordFlip" %in% class(p_flip$coordinates))
    }

    # Test 23: Test no local issues scenario
    diag_no_local <- diag_result
    diag_no_local$local_annotation_inconsistencies <-
        diag_no_local$local_annotation_inconsistencies[0, ]

    p_no_local <- plot(diag_no_local, plot_type = "local_issues")
    expect_s3_class(p_no_local, "ggplot")

    # Test 24: Test no annotation issues scenario
    diag_no_issues <- diag_result
    diag_no_issues$high_query_prop_analysis <- diag_no_issues$high_query_prop_analysis[0, ]
    diag_no_issues$cross_type_mixing <- diag_no_issues$cross_type_mixing[0, ]
    diag_no_issues$local_annotation_inconsistencies <-
        diag_no_issues$local_annotation_inconsistencies[0, ]

    p_no_issues <- plot(diag_no_issues, plot_type = "annotation_issues")
    expect_s3_class(p_no_issues, "ggplot")

    # Test 25: Test plot building without errors
    expect_silent(ggplot_build(plot(diag_result, plot_type = "community_data")))
    expect_silent(ggplot_build(plot(diag_result, plot_type = "summary")))

    # Test 26: Test different point sizes
    p_small <- plot(diag_result, plot_type = "cell_network", point_size = 0.5)
    p_large <- plot(diag_result, plot_type = "cell_network", point_size = 2.0)
    expect_s3_class(p_small, "ggplot")
    expect_s3_class(p_large, "ggplot")

    # Test 27: Test all plot types work
    plot_types <- c("community_network", "cell_network", "community_data",
                    "summary", "local_issues", "annotation_issues")

    for (plot_type in plot_types) {
        p_test <- plot(diag_result, plot_type = plot_type)
        expect_s3_class(p_test, "ggplot")
    }

    # Test 28: Test both color_by options for network plots
    color_options <- c("cell_type", "community_type")

    for (color_opt in color_options) {
        p_comm <- plot(diag_result, plot_type = "community_network", color_by = color_opt)
        p_cell <- plot(diag_result, plot_type = "cell_network", color_by = color_opt)
        expect_s3_class(p_comm, "ggplot")
        expect_s3_class(p_cell, "ggplot")
    }

    # Test 29: Test additional arguments
    p_extra <- plot(diag_result, plot_type = "summary",
                    environment = environment())
    expect_s3_class(p_extra, "ggplot")
})
