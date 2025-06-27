# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Load example datasets
data("reference_data")
data("query_data")

# Test cases
test_that("calculateMMDPValue correctly handles typical inputs", {

    # Basic functionality test
    result <- calculateMMDPValue(reference_data = reference_data,
                                 query_data = query_data,
                                 ref_cell_type_col = "expert_annotation",
                                 query_cell_type_col = "expert_annotation",
                                 n_permutation = 5)

    # Check that result is a named numeric vector
    expect_true(is.numeric(result))
    expect_true(!is.null(names(result)))

    # Check that p-values are in valid range [0, 1] or NA
    valid_pvals <- result[!is.na(result)]
    expect_true(all(valid_pvals >= 0 & valid_pvals <= 1))

    # Check that cell type names are present
    expect_true(length(names(result)) > 0)
})

test_that("calculateMMDPValue works with different parameters", {

    # Test with specific cell types
    query_types <- unique(colData(query_data)[["expert_annotation"]])
    ref_types <- unique(colData(reference_data)[["expert_annotation"]])
    common_types <- intersect(query_types, ref_types)

    if (length(common_types) > 1) {
        selected_types <- common_types[1:min(2, length(common_types))]

        result <- calculateMMDPValue(reference_data = reference_data,
                                     query_data = query_data,
                                     ref_cell_type_col = "expert_annotation",
                                     query_cell_type_col = "expert_annotation",
                                     cell_types = selected_types,
                                     n_permutation = 5)

        expect_true(is.numeric(result))
        expect_true(all(names(result) %in% selected_types))
    }

    # Test with different PC subset
    result2 <- calculateMMDPValue(reference_data = reference_data,
                                  query_data = query_data,
                                  ref_cell_type_col = "expert_annotation",
                                  query_cell_type_col = "expert_annotation",
                                  pc_subset = 1:3,
                                  n_permutation = 5)

    expect_true(is.numeric(result2))

    # Test with fewer permutations
    result3 <- calculateMMDPValue(reference_data = reference_data,
                                  query_data = query_data,
                                  ref_cell_type_col = "expert_annotation",
                                  query_cell_type_col = "expert_annotation",
                                  n_permutation = 5)

    expect_true(is.numeric(result3))
})

test_that("calculateMMDPValue handles edge cases", {

    # Test with very small permutation number
    result <- calculateMMDPValue(reference_data = reference_data,
                                 query_data = query_data,
                                 ref_cell_type_col = "expert_annotation",
                                 query_cell_type_col = "expert_annotation",
                                 n_permutation = 5)

    expect_true(is.numeric(result))

    # Test error with invalid column names
    expect_error(
        calculateMMDPValue(reference_data = reference_data,
                           query_data = query_data,
                           ref_cell_type_col = "nonexistent_column",
                           query_cell_type_col = "expert_annotation",
                           n_permutation = 5)
    )

    # Test error with invalid assay name
    expect_error(
        calculateMMDPValue(reference_data = reference_data,
                           query_data = query_data,
                           ref_cell_type_col = "expert_annotation",
                           query_cell_type_col = "expert_annotation",
                           assay_name = "nonexistent_assay",
                           n_permutation = 5)
    )
})

