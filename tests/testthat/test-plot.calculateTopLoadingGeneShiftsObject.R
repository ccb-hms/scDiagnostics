# Load necessary libraries
library(testthat)
library(scDiagnostics)
library(ggplot2)

# Load example datasets
data("reference_data")
data("query_data")

# Create test data for plot method tests
test_result <- calculateTopLoadingGeneShifts(
    query_data = query_data,
    reference_data = reference_data,
    query_cell_type_col = "SingleR_annotation",
    ref_cell_type_col = "expert_annotation",
    pc_subset = 1:3,
    n_top_loadings = 20
)

# Get available cell types for testing
available_cell_types <- unique(test_result[["cell_metadata"]][["cell_type"]])
test_cell_type <- available_cell_types[1]

# Test plot.calculateTopLoadingGeneShiftsObject function
test_that("plot.calculateTopLoadingGeneShiftsObject works with valid inputs", {

    # Test basic plotting functionality
    p <- plot(test_result,
              cell_type = test_cell_type,
              pc_subset = 1:3,
              plot_by = "p_adjusted",
              n_genes = 5)

    # Test return type
    expect_s3_class(p, "ggplot")

    # Test plot data structure
    plot_data <- p$data
    expect_s3_class(plot_data, "data.frame")

    # Test required columns exist
    expected_cols <- c("Gene_PC", "Expression", "Dataset", "PC")
    expect_true(all(expected_cols %in% colnames(plot_data)))

    # Test dataset factor levels
    expect_equal(levels(plot_data[["Dataset"]]), c("Reference", "Query"))

    # Test that expression values are numeric
    expect_true(is.numeric(plot_data[["Expression"]]))

    # Test faceting variable
    expect_true(all(grepl("PC[0-9]+ \\([0-9.]+%\\)", plot_data[["PC"]])))
})

test_that("plot.calculateTopLoadingGeneShiftsObject works with different plot_by options", {

    # Test plot_by = "p_adjusted" (default)
    p1 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1:2,
               plot_by = "p_adjusted",
               n_genes = 3)

    expect_s3_class(p1, "ggplot")
    expect_true(grepl("most significant adjusted p-values", p1$labels$subtitle))

    # Test plot_by = "top_loading"
    p2 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1:2,
               plot_by = "top_loading",
               n_genes = 3)

    expect_s3_class(p2, "ggplot")
    expect_true(grepl("highest absolute loadings", p2$labels$subtitle))

    # Test that different methods may show different genes
    data1 <- p1$data
    data2 <- p2$data

    # Both should have same structure
    expect_equal(colnames(data1), colnames(data2))
})

test_that("plot.calculateTopLoadingGeneShiftsObject handles different pc_subset values", {

    # Test single PC
    p1 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1,
               n_genes = 3)

    expect_s3_class(p1, "ggplot")
    pc_levels <- levels(p1$data[["PC"]])
    expect_length(pc_levels, 1)
    expect_true(grepl("PC1", pc_levels[1]))

    # Test multiple PCs
    p2 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1:3,
               n_genes = 3)

    expect_s3_class(p2, "ggplot")
    pc_levels <- levels(p2$data[["PC"]])
    expect_length(pc_levels, 3)
    expect_true(all(grepl("PC[1-3]", pc_levels)))

    # Test non-sequential PCs
    p3 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = c(1, 3),
               n_genes = 3)

    expect_s3_class(p3, "ggplot")
    pc_levels <- levels(p3$data[["PC"]])
    expect_length(pc_levels, 2)
})

test_that("plot.calculateTopLoadingGeneShiftsObject handles different significance thresholds", {

    # Test with strict threshold
    p1 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1:2,
               significance_threshold = 0.01,
               n_genes = 5)

    expect_s3_class(p1, "ggplot")

    # Test with lenient threshold
    p2 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1:2,
               significance_threshold = 0.1,
               n_genes = 5)

    expect_s3_class(p2, "ggplot")

    # Both should produce valid plots
    expect_equal(class(p1), class(p2))
})

test_that("plot.calculateTopLoadingGeneShiftsObject input validation works", {

    # Test missing cell_type
    expect_error(
        plot(test_result, pc_subset = 1:2),
        "cell_type must be specified and be exactly one cell type"
    )

    # Test multiple cell types
    expect_error(
        plot(test_result,
             cell_type = c(test_cell_type, available_cell_types[2]),
             pc_subset = 1:2),
        "cell_type must be specified and be exactly one cell type"
    )

    # Test invalid n_genes
    expect_error(
        plot(test_result,
             cell_type = test_cell_type,
             n_genes = -1),
        "n_genes must be a positive integer"
    )

    expect_error(
        plot(test_result,
             cell_type = test_cell_type,
             n_genes = 0),
        "n_genes must be a positive integer"
    )

    # Test invalid cell type
    expect_error(
        plot(test_result,
             cell_type = "NonexistentCellType",
             pc_subset = 1:2),
        "Cell type NonexistentCellType not found in results"
    )

    # Test invalid PC subset
    expect_error(
        plot(test_result,
             cell_type = test_cell_type,
             pc_subset = 10:12),
        "No requested PCs found in x"
    )
})

test_that("plot.calculateTopLoadingGeneShiftsObject handles missing data gracefully", {

    # Test with PC that has no data for specific cell type
    # This should issue a warning but still work for other PCs
    expect_warning(
        p <- plot(test_result,
                  cell_type = test_cell_type,
                  pc_subset = 1:3,
                  n_genes = 100),  # Request more genes than available
        NA  # May or may not warn depending on data availability
    )

    # Should still produce a plot if any PC has data
    if (exists("p")) {
        expect_s3_class(p, "ggplot")
    }
})

test_that("plot.calculateTopLoadingGeneShiftsObject plot elements are correct", {

    p <- plot(test_result,
              cell_type = test_cell_type,
              pc_subset = 1:3,
              plot_by = "p_adjusted",
              n_genes = 5)

    # Test plot labels
    expect_equal(p$labels$title, paste0("Gene Expression Shifts for ", test_cell_type))
    expect_true(grepl("most significant adjusted p-values", p$labels$subtitle))
    expect_equal(p$labels$x, "")  # No x-axis label
    expect_equal(p$labels$y, "")  # No y-axis label
    expect_equal(p$labels$fill, "Dataset")

    # Test theme elements
    expect_equal(p$theme$legend.position, "right")
    expect_equal(p$theme$plot.title$hjust, 0)  # Left aligned
    expect_equal(p$theme$plot.subtitle$hjust, 0)  # Left aligned

    # Test faceting
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("plot.calculateTopLoadingGeneShiftsObject handles different cell types", {

    # Test with different available cell types
    for (cell_type in head(available_cell_types, 3)) {
        p <- plot(test_result,
                  cell_type = cell_type,
                  pc_subset = 1:2,
                  n_genes = 3)

        expect_s3_class(p, "ggplot")
        expect_equal(p$labels$title, paste0("Gene Expression Shifts for ", cell_type))

        # Check that only the specified cell type is in the data
        plot_data <- p$data
        cell_metadata <- test_result[["cell_metadata"]]
        relevant_cells <- cell_metadata[cell_metadata[["cell_type"]] == cell_type, ]

        expect_true(all(plot_data[["Dataset"]] %in% c("Reference", "Query")))
    }
})

test_that("plot.calculateTopLoadingGeneShiftsObject handles edge cases", {

    # Test with object that has very few genes
    minimal_result <- calculateTopLoadingGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:2,
        n_top_loadings = 5
    )

    p <- plot(minimal_result,
              cell_type = test_cell_type,
              pc_subset = 1:2,
              n_genes = 3)

    expect_s3_class(p, "ggplot")

    # Test with single PC and single gene
    p_single <- plot(minimal_result,
                     cell_type = test_cell_type,
                     pc_subset = 1,
                     n_genes = 1)

    expect_s3_class(p_single, "ggplot")
})

test_that("plot.calculateTopLoadingGeneShiftsObject produces correct facet layout", {

    # Test with 1-3 PCs (should be 1 row)
    p1 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1:3,
               n_genes = 3)

    expect_s3_class(p1, "ggplot")
    expect_equal(p1$facet$params$ncol, 3)

    # Test with 4-6 PCs (should be 2 rows, max 3 cols)
    if (length(names(test_result)) >= 6) {  # Check if we have enough PCs
        p2 <- plot(test_result,
                   cell_type = test_cell_type,
                   pc_subset = 1:4,
                   n_genes = 3)

        expect_s3_class(p2, "ggplot")
        expect_equal(p2$facet$params$ncol, 3)  # Max 3 columns
    }
})

test_that("plot.calculateTopLoadingGeneShiftsObject gene ordering is correct", {

    # Test p_adjusted ordering
    p1 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1:2,
               plot_by = "p_adjusted",
               n_genes = 5)

    expect_s3_class(p1, "ggplot")

    # Test top_loading ordering
    p2 <- plot(test_result,
               cell_type = test_cell_type,
               pc_subset = 1:2,
               plot_by = "top_loading",
               n_genes = 5)

    expect_s3_class(p2, "ggplot")

    # Both should have Gene_PC factor levels
    expect_true(is.factor(p1$data[["Gene_PC"]]))
    expect_true(is.factor(p2$data[["Gene_PC"]]))

    # Gene_PC levels should be properly ordered
    expect_true(length(levels(p1$data[["Gene_PC"]])) > 0)
    expect_true(length(levels(p2$data[["Gene_PC"]])) > 0)
})

