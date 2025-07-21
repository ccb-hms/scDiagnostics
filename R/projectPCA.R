#' @title Project Query Data Onto PCA Space of Reference Data
#'
#' @description
#' This function projects a query singleCellExperiment object onto the PCA space of a reference
#' singleCellExperiment object. The PCA analysis on the reference data is assumed to be pre-computed
#' and stored within the object.
#'
#' @details
#' This function assumes that the "PCA" element exists within the \code{reducedDims} of the reference data
#' (obtained using \code{reducedDim(reference_data)}) and that the genes used for PCA are present in both
#' the reference and query data. It performs centering and scaling of the query data based on the reference
#' data before projection.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix
#' for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix
#' for the reference cells.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data}
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data}
#' that identifies the cell types.
#' @param pc_subset A numeric vector specifying the subset of principal components (PCs) to compare. Default is 1:10.
#' @param assay_name Name of the assay on which to perform computations. Defaults to \code{"logcounts"}.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A \code{data.frame} containing the projected data in rows (reference and query data combined).
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Project the query data onto PCA space of reference
#' pca_output <- projectPCA(query_data = query_data,
#'                          reference_data = reference_data,
#'                          query_cell_type_col = "SingleR_annotation",
#'                          ref_cell_type_col = "expert_annotation",
#'                          pc_subset = 1:10)
#'
# Function to project query data onto PCA space of reference data
projectPCA <- function(query_data,
                       reference_data,
                       query_cell_type_col,
                       ref_cell_type_col,
                       pc_subset = 1:10,
                       assay_name = "logcounts",
                       max_cells = 2500){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Downsample query and reference data
    query_data <- downsampleSCE(sce = query_data,
                                max_cells = max_cells)
    reference_data <- downsampleSCE(sce = reference_data,
                                    max_cells = max_cells)

    # Extract reference PCA components and rotation matrix
    ref_mat <- reducedDim(reference_data, "PCA")[, pc_subset, drop = FALSE]
    rotation_mat <- attributes(reducedDim(
        reference_data, "PCA"))[["rotation"]][, pc_subset, drop = FALSE]
    PCA_genes <- rownames(rotation_mat)

    # Check if genes used for PCA are available in query data
    if (!all(PCA_genes %in% rownames(assay(query_data)))) {
        stop("Genes in reference PCA are not found in query data.")
    }

    # Center and scale query data based on reference for projection
    centering_vec <- apply(t(as.matrix(
        assay(reference_data, assay_name))), 2, mean)[PCA_genes]
    query_mat <- scale(t(as.matrix(
        assay(query_data, assay_name)))[, PCA_genes, drop = FALSE],
                       center = centering_vec, scale = FALSE) %*% rotation_mat

    # Returning output as a dataframe
    return(data.frame(rbind(ref_mat, query_mat),
                      dataset = c(rep("Reference", nrow(ref_mat)),
                                  rep("Query", nrow(query_mat))),
                      cell_type = c(ifelse(rep(is.null(ref_cell_type_col),
                                               nrow(ref_mat)),
                                           rep(NA, nrow(ref_mat)),
                                           reference_data[[ref_cell_type_col]]),
                                    ifelse(rep(is.null(query_cell_type_col),
                                               nrow(query_mat)),
                                           rep(NA, nrow(query_mat)),
                                           query_data[[query_cell_type_col]]))))
}
