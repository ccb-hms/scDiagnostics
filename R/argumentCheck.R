#' @title Argument Validation for SingleCellExperiment Analysis
#'
#' @description
#' This function validates the input arguments for functions that analyze \code{\linkS4class{SingleCellExperiment}}
#' objects. It checks that the inputs are of the correct types and formats, and that required columns and cell types
#' are present in the data.
#'
#' @details
#' The function performs a series of checks to ensure that:
#' \itemize{
#'  \item `query_data` and `reference_data` are \code{\linkS4class{SingleCellExperiment}} objects.
#'  \item `query_cell_type_col` and `ref_cell_type_col` exist in the column data of their respective \code{\linkS4class{SingleCellExperiment}} objects.
#'  \item The specified `cell_types` are available in the provided datasets.
#'  \item If `unique_cell_type` is `TRUE`, there should only be one cell type in the \code{\linkS4class{SingleCellExperiment}} objects.
#'  \item If `plot_function` is `TRUE`, the number of unique `cell_types` does not exceed 10.
#'  \item `cell_names_query` are valid cell names in the provided query dataset.
#'  \item `cell_names_ref` are valid cell names in the provided reference dataset.
#'  \item The PCA subsets specified by `pc_subset_query` and `pc_subset_ref` are valid.
#' }
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' If `NULL`, no check is performed.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' If `NULL`, no check is performed.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data}
#' that identifies the cell types. If `NULL`, no check is performed.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data}
#' that identifies the cell types. If `NULL`, no check is performed.
#' @param cell_types A character vector specifying the cell types to include in the plot. If `NULL`, no check is performed.
#' @param unique_cell_type If `TRUE`, there should only be one cell type in the provided \code{\linkS4class{SingleCellExperiment}} objects.
#' Default is `FALSE`.
#' @param plot_function A logical value indicating whether the function is being called to generate a plot. Default is `FALSE`.
#' @param cell_names_query A character vector of cell names in query data to be analyzed. If `NULL`, no check is performed.
#' @param cell_names_ref A character vector of cell names in reference data to be analyzed. If `NULL`, no check is performed.
#' @param pc_subset_query A numeric vector specifying the principal components to be used for the query data. If `NULL`, no check is performed.
#' @param pc_subset_ref A numeric vector specifying the principal components to be used for the reference data. If `NULL`, no check is performed.
#' @param common_rotation_genes If TRUE, check the rotation matrices of the reference and query data and ensure they have the same genes.
#' Default is FALSE.
#' @param assay_name Name of the assay on which to perform computations. If `NULL`, no check is performed.
#'
#' @keywords internal
#'
#' @return None.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to check standard arguments for functions in the package
argumentCheck <- function(query_data = NULL,
                          reference_data = NULL,
                          query_cell_type_col = NULL,
                          ref_cell_type_col = NULL,
                          cell_types = NULL,
                          unique_cell_type = FALSE,
                          plot_function = FALSE,
                          cell_names_query = NULL,
                          cell_names_ref = NULL,
                          pc_subset_query = NULL,
                          pc_subset_ref = NULL,
                          common_rotation_genes = FALSE,
                          assay_name = NULL) {

    # Check if query_data is a SingleCellExperiment object
    if (!is.null(query_data)) {

        if (!is(query_data, "SingleCellExperiment")) {
            stop("'query_data' must be a SingleCellExperiment object.")
        }

        if(!is.null(assay_name) && !(assay_name %in% SummarizedExperiment::assayNames(query_data)))
            stop("'query_data' does not contain the specified assay.")
    }

    # Check if reference_data is a SingleCellExperiment object
    if (!is.null(reference_data)) {

        if (!is(reference_data, "SingleCellExperiment")) {
            stop("'reference_data' must be a SingleCellExperiment object.")
        }

        if(!is.null(assay_name) && !(assay_name %in% SummarizedExperiment::assayNames(reference_data)))
            stop("'reference_data' does not contain the specified assay.")
    }

    # Check if query_cell_type_col is a character string of length 1 and exists in query_data
    if (!is.null(query_cell_type_col)) {

        if (!is.null(query_data)) {
            if (!is.character(query_cell_type_col) ||
                length(query_cell_type_col) != 1) {
                stop("'query_cell_type_col' must be a character string of length 1")
            }

            if (!(query_cell_type_col %in% names(colData(query_data)))) {
                stop("'query_cell_type_col' is not an existing column of 'query_data'.")
            }
        }
    }

    # Check if ref_cell_type_col is a character string of length 1 and exists in reference_data
    if (!is.null(ref_cell_type_col)) {

        if (!is.null(reference_data)) {
            if (!is.character(ref_cell_type_col) ||
                length(ref_cell_type_col) != 1) {
                stop("'ref_cell_type_col' must be a character string of length 1")
            }

            if (!(ref_cell_type_col %in% names(colData(reference_data)))) {
                stop("'ref_cell_type_col' is not an existing column of 'reference_data'.")
            }
        }
    }

    # Check if cell_types are available in the SingleCellExperiment object(s)
    if (!is.null(cell_types)) {

        if (!is.null(query_data)) {
            if (!all(cell_types %in%
                     unique(query_data[[query_cell_type_col]]))) {
                stop("'cell_types' contains one or more cell types that are not available in 'query_data'.")
            }
        }

        if (!is.null(reference_data)) {
            if (!all(cell_types %in%
                     unique(reference_data[[ref_cell_type_col]]))) {
                stop("'cell_types' contains one or more cell types that are not available in 'reference_data'.")
            }
        }
    }

    # Check that the SingleCellExperiment object(s) have a unique cell type
    if (isTRUE(unique_cell_type)) {

        if (!is.null(query_data)) {
            if (length(unique(query_data[[query_cell_type_col]])) > 1) {
                stop("This function should be used when there is only one cell type in 'query_data'.")
            }
        }

        if (!is.null(reference_data)) {

            if (length(unique(reference_data[[ref_cell_type_col]])) > 1) {
                stop("This function should be used when there is only one cell type in 'reference_data'.")
            }
        }

        if (!is.null(reference_data) && !is.null(query_data)) {
            if (unique(query_data[[query_cell_type_col]]) != unique(reference_data[[ref_cell_type_col]])) {
                stop("The cell type of the query data does not match the cell type of the reference data.")
            }
        }
    }

    # Check the number of cell types for plot function
    if (plot_function == TRUE) {

        if (length(unique(cell_types)) > 10) {
            stop("The maximum number of cell types for plotting is 10.")
        }
    }

    # Check cell_names contain valid cell names in query_data
    if (!is.null(cell_names_query)) {

        if (!all(cell_names_query %in% colnames(query_data))) {
            stop("'cell_names' contains one or more cells that are not available in 'query_data'.")
        }
    }

    # Check cell_names contain valid cell names in reference_data
    if (!is.null(cell_names_ref)) {

        if (!all(cell_names_ref %in% colnames(reference_data))) {
            stop("'cell_names' contains one or more cells that are not available in 'reference_data'.")
        }
    }

    # Check PC subset for query_data
    if (!is.null(pc_subset_query)) {

        # Check if "PCA" is present in query's reduced dimensions
        if (!"PCA" %in% names(reducedDims(query_data))) {
            stop("'query_data' must have pre-computed PCA in 'reducedDims'. ",
                 "Use processPCA() to compute PCA: ",
                 "query_data <- processPCA(query_data = query_data)")
        }

        # Check input if PC subset is valid
        if (!all(pc_subset_query %in% seq_len(ncol(reducedDim(query_data, "PCA"))))) {
            stop("'pc_subset' is out of range for 'query_data'.")
        }
    }

    # Check PC subset for reference_data
    if (!is.null(pc_subset_ref)) {

        # Check if "PCA" is present in reference's reduced dimensions
        if (!"PCA" %in% names(reducedDims(reference_data))) {
            stop("'reference_data' must have pre-computed PCA in 'reducedDims'. ",
                 "Use processPCA() to compute PCA: ",
                 "reference_data <- processPCA(reference_data = reference_data)")
        }

        # Check input if PC subset is valid
        if (!all(pc_subset_ref %in% seq_len(ncol(reducedDim(reference_data, "PCA"))))) {
            stop("'pc_subset' is out of range for 'reference_data'.")
        }
    }

    # Check if the rotation matrices have the same genes in the same order
    if (common_rotation_genes == TRUE) {

        # Check if both datasets have PCA before comparing rotation matrices
        if (!"PCA" %in% names(reducedDims(query_data))) {
            stop("'query_data' must have pre-computed PCA in 'reducedDims' for rotation matrix comparison. ",
                 "Use processPCA() to compute PCA for both datasets: ",
                 "result <- processPCA(query_data = query_data, reference_data = reference_data)")
        }

        if (!"PCA" %in% names(reducedDims(reference_data))) {
            stop("'reference_data' must have pre-computed PCA in 'reducedDims' for rotation matrix comparison. ",
                 "Use processPCA() to compute PCA for both datasets: ",
                 "result <- processPCA(query_data = query_data, reference_data = reference_data)")
        }

        # Check if the rotation matrices have the same number of genes
        if (ncol(attributes(reducedDim(query_data, "PCA"))[["rotation"]]) !=
            ncol(attributes(reducedDim(reference_data, "PCA"))[["rotation"]])) {
            stop("The number of genes in the rotation matrices differ.")
        }

        # Check if genes in both rotation matrices are the same
        if (!all(rownames(attributes(reducedDim(query_data, "PCA"))[["rotation"]]) %in%
                 rownames(attributes(reducedDim(reference_data, "PCA"))[["rotation"]]))) {
            stop("The genes in the rotation matrices differ.")
        }
    }
}
