# Load necessary libraries
library(testthat)
library(scDiagnostics)

# Unit tests for calculateHVGOverlap function
test_that("calculateHVGOverlap function tests", {
    
    # Test case 1: Check if function returns a numeric value
    reference_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
    query_genes <- c("Gene3", "Gene4", "Gene5", "Gene6", "Gene7")
    overlap_coefficient <- calculateHVGOverlap(reference_genes, query_genes)
    expect_type(overlap_coefficient, "double")
    
    # Test case 2: Check if function returns a value between 0 and 1
    expect_true(overlap_coefficient >= 0 && overlap_coefficient <= 1)
    
    # Test case 3: Check if function handles complete overlap correctly
    reference_genes <- c("Gene1", "Gene2", "Gene3")
    query_genes <- c("Gene1", "Gene2", "Gene3")
    overlap_coefficient <- calculateHVGOverlap(reference_genes, query_genes)
    expect_equal(overlap_coefficient, 1)
    
    # Test case 4: Check if function handles no overlap correctly
    reference_genes <- c("Gene1", "Gene2", "Gene3")
    query_genes <- c("Gene4", "Gene5", "Gene6")
    overlap_coefficient <- calculateHVGOverlap(reference_genes, query_genes)
    expect_equal(overlap_coefficient, 0)
    
    # Test case 5: Check if function handles partial overlap correctly
    reference_genes <- c("Gene1", "Gene2", "Gene3")
    query_genes <- c("Gene2", "Gene3", "Gene4")
    overlap_coefficient <- calculateHVGOverlap(reference_genes, query_genes)

    # Test case 6: Check if function handles different sizes of input vectors
    reference_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
    query_genes <- c("Gene3", "Gene4")
    overlap_coefficient <- calculateHVGOverlap(reference_genes, query_genes)
    expect_equal(overlap_coefficient, 2 / 2)
    
    # Test case 7: Check for correct error handling of non-character input
    expect_error(calculateHVGOverlap(reference_genes = 1:5, query_genes = c("Gene1", "Gene2")),
                 "reference_genes must be a character vector.")
    expect_error(calculateHVGOverlap(reference_genes = c("Gene1", "Gene2"), query_genes = 1:5),
                 "query_genes must be a character vector.")
    
    # Test case 8: Check for correct error handling of empty input vectors
    expect_error(calculateHVGOverlap(reference_genes = character(0), query_genes = c("Gene1", "Gene2")),
                 "Input vectors must not be empty.")
    expect_error(calculateHVGOverlap(reference_genes = c("Gene1", "Gene2"), query_genes = character(0)),
                 "Input vectors must not be empty.")
})
