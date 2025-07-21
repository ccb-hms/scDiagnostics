#' @title Project Query Data Onto SIR Space of Reference Data
#'
#' @description
#' This function projects a query \code{SingleCellExperiment} object onto the SIR (supervised independent
#' component) space of a reference \code{SingleCellExperiment} object. The SVD of the reference data is
#' computed on conditional means per cell type, and the query data is projected based on these reference
#' components.
#'
#' @details
#' The genes used for the projection (SVD) must be present in both the reference and query datasets.
#' The function first computes conditional means for each cell type in the reference data, then performs
#' SVD on these conditional means to obtain the rotation matrix used for projecting both the reference
#' and query datasets. The query data is centered and scaled based on the reference data.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col A character string specifying the column in the \code{colData} of \code{query_data}
#' that identifies the cell types.
#' @param ref_cell_type_col A character string specifying the column in the \code{colData} of \code{reference_data}
#' that identifies the cell types.
#' @param cell_types A character vector of cell types for which to compute conditional means in the reference data.
#' @param multiple_cond_means A logical value indicating whether to compute multiple conditional means per cell type
#' (through PCA and clustering). Defaults to \code{TRUE}.
#' @param cumulative_variance_threshold A numeric value between 0 and 1 specifying the variance threshold for PCA
#' when computing multiple conditional means. Defaults to \code{0.7}.
#' @param n_neighbor An integer specifying the number of nearest neighbors for clustering when computing multiple
#' conditional means. Defaults to \code{1}.
#' @param assay_name A character string specifying the assay name on which to perform computations. Defaults to \code{"logcounts"}.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A list containing:
#' \item{cond_means}{A matrix of the conditional means computed for the reference data.}
#' \item{rotation_mat}{The rotation matrix obtained from the SVD of the conditional means.}
#' \item{sir_projections}{A \code{data.frame} containing the SIR projections for both the reference and query datasets.}
#' \item{percent_var}{The percentage of variance explained by each component of the SIR projection.}
#'
#' @export
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Project the query data onto SIR space of reference
#' sir_output <- projectSIR(query_data = query_data,
#'                          reference_data = reference_data,
#'                          query_cell_type_col = "SingleR_annotation",
#'                          ref_cell_type_col = "expert_annotation")
#'
# Function to project query data onto PCA space of reference data
projectSIR <- function(query_data,
                       reference_data,
                       query_cell_type_col,
                       ref_cell_type_col,
                       cell_types = NULL,
                       multiple_cond_means = TRUE,
                       cumulative_variance_threshold = 0.7,
                       n_neighbor = 1,
                       assay_name = "logcounts",
                       max_cells = 2500){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  assay_name = assay_name)

    # Downsample query and reference data
    query_data <- downsampleSCE(sce = query_data,
                                max_cells = max_cells)
    reference_data <- downsampleSCE(sce = reference_data,
                                    max_cells = max_cells)

    # Check if cumulative_variance_threshold is between 0 and 1
    if (!is.numeric(cumulative_variance_threshold) ||
        cumulative_variance_threshold < 0 || cumulative_variance_threshold > 1) {
        stop("cumulative_variance_threshold must be a numeric value between 0 and 1.")
    }

    # Check if n_neighbor is a positive integer
    if (!is.numeric(n_neighbor) || n_neighbor <= 0 || n_neighbor != as.integer(n_neighbor)) {
        stop("n_neighbor must be a positive integer.")
    }

    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }

    # Compute conditional means for each cell type of reference data
    cond_means <- conditionalMeans(reference_data = reference_data,
                                   ref_cell_type_col = ref_cell_type_col,
                                   cell_types = cell_types,
                                   multiple_cond_means = multiple_cond_means,
                                   assay_name = assay_name,
                                   cumulative_variance_threshold = cumulative_variance_threshold,
                                   n_neighbor = n_neighbor)

    # Extract reference SVD components and rotation matrix
    SVD_genes <- rownames(reference_data)
    centering_vec <- apply(cond_means, 2, mean)[SVD_genes]
    svd_ref <- svd(scale(cond_means, center = centering_vec, scale = FALSE))
    rotation_mat <- svd_ref[["v"]]
    rownames(rotation_mat) <- SVD_genes

    # Check if genes used for PCA are available in query data
    if (!all(SVD_genes %in% rownames(assay(query_data, assay_name)))) {
        stop("Genes in reference SVD are not found in query data.")
    }

    # Compute percentage of variance explained
    percent_var <- svd_ref[["d"]]^2 /
        sum(svd_ref[["d"]]^2) * 100

    # Center and scale reference query data based on reference for projection
    ref_mat <- scale(t(as.matrix(
        assay(reference_data, assay_name)))[, SVD_genes, drop = FALSE],
        center = centering_vec, scale = FALSE) %*% rotation_mat
    query_mat <- scale(t(as.matrix(
        assay(query_data, assay_name)))[, SVD_genes, drop = FALSE],
        center = centering_vec, scale = FALSE) %*% rotation_mat
    sir_projections <- data.frame(
        rbind(ref_mat, query_mat),
        dataset = c(rep("Reference", nrow(ref_mat)),
                    rep("Query", nrow(query_mat))),
        cell_type = c(ifelse(rep(is.null(ref_cell_type_col),
                                 nrow(ref_mat)),
                             rep(NA, nrow(ref_mat)),
                             reference_data[[ref_cell_type_col]]),
                      ifelse(rep(is.null(query_cell_type_col),
                                 nrow(query_mat)),
                             rep(NA, nrow(query_mat)),
                             query_data[[query_cell_type_col]])))
    colnames(sir_projections)[seq_len(ncol(ref_mat))] <-
        paste0("SIR", seq_len(ncol(ref_mat)), " (",
               sprintf("%.1f%%", percent_var[seq_len(ncol(ref_mat))]), ")")

    # Returning output as a data frame
    return(list(cond_means = cond_means,
                rotation_mat = rotation_mat,
                sir_projections = sir_projections,
                percent_var = percent_var))
}

#' @title Compute Conditional Means for Cell Types
#'
#' @description
#' This function computes conditional means for each cell type in the reference data. It can compute
#' either a single conditional mean per cell type or multiple conditional means, depending on the
#' specified settings. Principal component analysis (PCA) is used for dimensionality reduction before
#' clustering when computing multiple conditional means.
#'
#' @details
#' The function offers two modes of operation:
#' - **Single conditional mean per cell type**: For each cell type, it computes the mean expression
#'   across all observations.
#' - **Multiple conditional means per cell type**: For each cell type, the function performs PCA to
#'   reduce dimensionality, followed by clustering to compute multiple conditional means.
#'
#' @param reference_data A \code{SingleCellExperiment} object containing the reference data, where rows
#' represent genes and columns represent cells.
#' @param ref_cell_type_col A character string specifying the column in \code{colData(reference_data)} that contains the cell type labels.
#' @param cell_types A character vector of cell types for which to compute conditional means.
#' @param multiple_cond_means A logical value indicating whether to compute multiple conditional means per cell type.
#' Defaults to \code{FALSE}.
#' @param assay_name A character string specifying the name of the assay to use for the computation. Defaults to \code{"logcounts"}.
#' @param cumulative_variance_threshold A numeric value between 0 and 1 specifying the variance threshold
#' for PCA when computing multiple conditional means. Defaults to \code{0.7}.
#' @param n_neighbor An integer specifying the number of nearest neighbors for clustering when computing multiple conditional means.
#' Defaults to \code{1}.
#'
#' @return A numeric matrix with the conditional means for each cell type. If \code{multiple_cond_means = TRUE}, the matrix
#' will contain multiple rows for each cell type, representing the different conditional means computed via clustering.
#'
#' @keywords internal
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to compute Ledoit-Wolf covariance matrix
conditionalMeans <- function(reference_data,
                             ref_cell_type_col,
                             cell_types,
                             multiple_cond_means = FALSE,
                             assay_name = "logcounts",
                             cumulative_variance_threshold = 0.7,
                             n_neighbor = 1) {

    # Compute conditional means for each cell type of reference data
    if (multiple_cond_means) {

        # Matrix to store results
        cond_means <- matrix(nrow = 0, ncol = nrow(reference_data))
        colnames(cond_means) <- rownames(reference_data)

        for(cell_type in cell_types){

            # Compute multiple conditional means per cell type
            assay_mat <- scale(t(as.matrix(assay(
                reference_data[, which(reference_data[[ref_cell_type_col]] == cell_type)],
                assay_name))), center = TRUE, scale = FALSE)
            assay_svd <- svd(assay_mat)
            cumulative_variance <- cumsum(assay_svd$d^2) / sum(assay_svd$d^2)
            n_components <- min(which(cumulative_variance >= cumulative_variance_threshold))
            projections <- assay_mat %*%
                assay_svd$v[, seq_len(n_components)]
            clusters <- bluster::clusterRows(
                projections,
                BLUSPARAM = bluster::TwoStepParam(
                    second = bluster::NNGraphParam(
                        k = n_neighbor)))
            cluster_means <- do.call(
                rbind, lapply(unique(clusters),
                              function(cl) colMeans(assay_mat[clusters == cl,, drop = FALSE])))
            rownames(cluster_means) <- rep(cell_type, nrow(cluster_means))
            cond_means  <- rbind(cond_means, cluster_means)
        }

    } else {

        # Compute a single conditional mean per cell type
        cond_means <- lapply(cell_types, function(x)
            apply(as.matrix(assay(
                reference_data[, which(reference_data[[ref_cell_type_col]] == x)],
                assay_name)), 1, mean))
        cond_means <- do.call(rbind, cond_means)
        rownames(cond_means) <- cell_types

    }

    # Return conditional means
    return(cond_means)
}

