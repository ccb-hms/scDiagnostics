test_that("detectAnomaly handles a reference cell type absent from the query", {
    data("reference_data", package = "scDiagnostics", envir = environment())
    data("query_data", package = "scDiagnostics", envir = environment())

    # Drop one cell type from the query entirely, leaving it in the reference.
    dropped <- "CD4"
    keep <- query_data[["SingleR_annotation"]] != dropped
    skip_if_not(any(!keep), "test cell type not present in query_data")
    query_missing <- query_data[, keep]

    expect_no_error(
        res <- detectAnomaly(
            reference_data = reference_data,
            query_data = query_missing,
            ref_cell_type_col = "expert_annotation",
            query_cell_type_col = "SingleR_annotation",
            pc_subset = 1:5,
            n_tree = 50
        )
    )

    # The reference-only cell type is still reported, with reference scores and
    # an empty query score vector rather than an error.
    expect_true(dropped %in% names(res))
    expect_gt(length(res[[dropped]][["reference_anomaly_scores"]]), 0)
    expect_length(res[[dropped]][["query_anomaly_scores"]], 0)

    # Cell types present in both are unaffected.
    shared <- setdiff(names(res), c(dropped, "Combined"))
    expect_gt(length(shared), 0)
    for (ct in shared) {
        expect_gt(length(res[[ct]][["query_anomaly_scores"]]), 0)
    }
})

test_that("detectAnomaly skips cell types with fewer than two reference cells", {
    data("reference_data", package = "scDiagnostics", envir = environment())
    data("query_data", package = "scDiagnostics", envir = environment())

    types <- reference_data[["expert_annotation"]]
    target <- names(sort(table(types)))[1]
    # Keep a single cell of `target` in the reference.
    idx <- which(types == target)
    ref_thin <- reference_data[, c(setdiff(seq_along(types), idx), idx[1])]

    expect_warning(
        res <- detectAnomaly(
            reference_data = ref_thin,
            query_data = query_data,
            ref_cell_type_col = "expert_annotation",
            query_cell_type_col = "SingleR_annotation",
            pc_subset = 1:5,
            n_tree = 50
        ),
        "isolation forest"
    )
    expect_false(target %in% names(res))
})
