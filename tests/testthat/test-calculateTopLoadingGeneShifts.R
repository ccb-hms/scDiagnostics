# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Test calculateTopLoadingGeneShifts function
test_that("calculateTopLoadingGeneShifts works with valid inputs", {

    # Test with default parameters
    result <- calculateTopLoadingGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation"
    )

    # Test return structure
    expect_type(result, "list")
    expect_s3_class(result, "calculateTopLoadingGeneShiftsObject")

    # Test required elements exist
    expected_elements <- c("expression_data", "cell_metadata", "gene_metadata", "percent_var")
    expect_true(all(expected_elements %in% names(result)))

    # Test PC elements exist (default is 1:5)
    pc_names <- paste0("PC", 1:5)
    expect_true(all(pc_names %in% names(result)))

    # Test expression_data structure
    expect_true(nrow(result[["expression_data"]]) > 0)
    expect_true(ncol(result[["expression_data"]]) > 0)

    # Test cell_metadata structure
    expect_s3_class(result[["cell_metadata"]], "data.frame")
    expected_cell_cols <- c("cell_id", "dataset", "cell_type", "original_index")
    expect_true(all(expected_cell_cols %in% colnames(result[["cell_metadata"]])))
    expect_true(all(result[["cell_metadata"]][["dataset"]] %in% c("Reference", "Query")))

    # Test gene_metadata structure
    expect_s3_class(result[["gene_metadata"]], "data.frame")
    expected_gene_cols <- c("gene", "pc", "loading")
    expect_true(all(expected_gene_cols %in% colnames(result[["gene_metadata"]])))

    # Test percent_var structure
    expect_type(result[["percent_var"]], "double")
    expect_named(result[["percent_var"]], paste0("PC", 1:5))
    expect_true(all(result[["percent_var"]] >= 0))
    expect_true(all(result[["percent_var"]] <= 100))
})

test_that("calculateTopLoadingGeneShifts works with custom parameters", {

    # Test with custom pc_subset
    result <- calculateTopLoadingGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:3,
        n_top_loadings = 20,
        p_value_threshold = 0.01,
        adjust_method = "bonferroni"
    )

    # Test correct PC subset
    expect_true(all(paste0("PC", 1:3) %in% names(result)))
    expect_false("PC4" %in% names(result))
    expect_false("PC5" %in% names(result))

    # Test percent_var has correct names
    expect_named(result[["percent_var"]], paste0("PC", 1:3))

    # Test gene metadata has correct PCs
    expect_true(all(result[["gene_metadata"]][["pc"]] %in% 1:3))
})

test_that("calculateTopLoadingGeneShifts works with specific cell types", {

    # Get available cell types
    common_cell_types <- intersect(
        unique(query_data[["SingleR_annotation"]]),
        unique(reference_data[["expert_annotation"]])
    )

    # Test with specific cell types
    selected_cell_types <- head(common_cell_types, 2)
    result <- calculateTopLoadingGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        cell_types = selected_cell_types,
        pc_subset = 1:2
    )

    # Test that only selected cell types are present
    unique_cell_types <- unique(result[["cell_metadata"]][["cell_type"]])
    expect_true(all(unique_cell_types %in% selected_cell_types))

    # Test PC results only contain selected cell types
    for (pc_name in paste0("PC", 1:2)) {
        if (nrow(result[[pc_name]]) > 0) {
            expect_true(all(result[[pc_name]][["cell_type"]] %in% selected_cell_types))
        }
    }
})

test_that("calculateTopLoadingGeneShifts handles edge cases", {

    # Test with single PC
    result <- calculateTopLoadingGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1,
        n_top_loadings = 10
    )

    expect_true("PC1" %in% names(result))
    expect_false("PC2" %in% names(result))
    expect_length(result[["percent_var"]], 1)

    # Test with small n_top_loadings
    result_small <- calculateTopLoadingGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:2,
        n_top_loadings = 5
    )

    expect_true(nrow(result_small[["gene_metadata"]]) <= 10)  # 5 genes * 2 PCs
})

test_that("calculateTopLoadingGeneShifts input validation works", {

    # Test invalid n_top_loadings
    expect_error(
        calculateTopLoadingGeneShifts(
            query_data = query_data,
            reference_data = reference_data,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            n_top_loadings = -5
        ),
        "n_top_loadings must be a positive integer"
    )

    expect_error(
        calculateTopLoadingGeneShifts(
            query_data = query_data,
            reference_data = reference_data,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            n_top_loadings = 0
        ),
        "n_top_loadings must be a positive integer"
    )

    # Test invalid p_value_threshold
    expect_error(
        calculateTopLoadingGeneShifts(
            query_data = query_data,
            reference_data = reference_data,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            p_value_threshold = -0.1
        ),
        "p_value_threshold must be between 0 and 1"
    )

    expect_error(
        calculateTopLoadingGeneShifts(
            query_data = query_data,
            reference_data = reference_data,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            p_value_threshold = 1.5
        ),
        "p_value_threshold must be between 0 and 1"
    )

    # Test invalid cell type columns
    expect_error(
        calculateTopLoadingGeneShifts(
            query_data = query_data,
            reference_data = reference_data,
            query_cell_type_col = "nonexistent_column",
            ref_cell_type_col = "expert_annotation"
        )
    )

    expect_error(
        calculateTopLoadingGeneShifts(
            query_data = query_data,
            reference_data = reference_data,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "nonexistent_column"
        )
    )
})

test_that("calculateTopLoadingGeneShifts handles different adjust methods", {

    methods <- c("fdr", "bonferroni", "holm", "none")

    for (method in methods) {
        result <- calculateTopLoadingGeneShifts(
            query_data = query_data,
            reference_data = reference_data,
            query_cell_type_col = "SingleR_annotation",
            ref_cell_type_col = "expert_annotation",
            pc_subset = 1:2,
            n_top_loadings = 10,
            adjust_method = method
        )

        expect_s3_class(result, "calculateTopLoadingGeneShiftsObject")

        # Test that p_adjusted values are reasonable
        for (pc_name in paste0("PC", 1:2)) {
            if (nrow(result[[pc_name]]) > 0) {
                expect_true(all(result[[pc_name]][["p_adjusted"]] >= 0))
                expect_true(all(result[[pc_name]][["p_adjusted"]] <= 1))

                # For most methods, p_adjusted should be >= p_value
                if (method != "none") {
                    expect_true(all(result[[pc_name]][["p_adjusted"]] >=
                                        result[[pc_name]][["p_value"]]))
                }
            }
        }
    }
})


test_that("processGenesSimple helper function works correctly", {

    # Create test data
    test_genes <- c("GENE1", "GENE2")
    test_loadings <- c(0.5, -0.3)

    # Create test expression matrices
    query_matrix <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
    ref_matrix <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 2, ncol = 3)
    rownames(query_matrix) <- test_genes
    rownames(ref_matrix) <- test_genes

    result <- processGenesSimple(
        top_genes = test_genes,
        top_loadings = test_loadings,
        query_expr_matrix = query_matrix,
        ref_expr_matrix = ref_matrix,
        cell_type = "TestCellType"
    )

    expect_type(result, "list")
    expect_length(result, 2)

    # Test first gene result
    gene1_result <- result[[1]]
    expect_s3_class(gene1_result, "data.frame")
    expect_equal(gene1_result[["gene"]], "GENE1")
    expect_equal(gene1_result[["loading"]], 0.5)
    expect_equal(gene1_result[["cell_type"]], "TestCellType")
    expect_true(is.numeric(gene1_result[["p_value"]]))
    expect_true(is.numeric(gene1_result[["mean_query"]]))
    expect_true(is.numeric(gene1_result[["mean_reference"]]))
})

test_that("calculateTopLoadingGeneShifts preserves gene order by significance", {

    result <- calculateTopLoadingGeneShifts(
        query_data = query_data,
        reference_data = reference_data,
        query_cell_type_col = "SingleR_annotation",
        ref_cell_type_col = "expert_annotation",
        pc_subset = 1:2,
        n_top_loadings = 20
    )

    # Test that results are sorted by p_adjusted (ascending)
    for (pc_name in paste0("PC", 1:2)) {
        if (nrow(result[[pc_name]]) > 1) {
            p_values <- result[[pc_name]][["p_adjusted"]]
            expect_true(all(p_values == sort(p_values)))
        }
    }
})

