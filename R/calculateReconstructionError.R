#' @title Calculate PCA Reconstruction Errors for Out-of-Distribution Anomaly Detection
#'
#' @description
#' This function detects "out-of-distribution" anomalies by calculating the PCA reconstruction error
#' (Sum of Squared Errors) for each cell. It projects cells into a reference PCA space, attempts to
#' reconstruct their original gene expression profile based solely on reference PCA rules, and measures
#' the residual difference.
#'
#' @details
#' PCA creates a low-dimensional summary of biological variation. By computing a local PCA space
#' specifically for each reference cell type, the algorithm learns the strict biological rules governing
#' that specific cell state.
#'
#' If a query cell contains a novel biological state (e.g., a viral infection, unique drug response,
#' or it is actually an unrepresented cell subtype hiding in the cluster), it will express genes that
#' the local reference PCA ignores. When the query cell is projected into the reference PCA space and
#' mathematically reconstructed, those novel gene expressions are lost.
#'
#' By subtracting the reconstructed matrix from the original matrix, this function isolates the
#' "Residuals" (biology that the reference cannot explain). The Sum of Squared Errors (SSE) of
#' these residuals serves as a highly sensitive anomaly score for novel biological states.
#'
#' @param reference_data A \code{\link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_data An optional \code{\link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' If NULL, the reconstruction errors are computed for the reference data alone. Default is NULL.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param cell_types A character vector specifying the cell types to analyze. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to use in the reconstruction. Default is 1:5.
#' @param n_hvgs An integer specifying the number of highly variable genes to calculate for each cell type's local PCA space. Default is 100.
#' @param mad_multiplier A numeric value specifying the number of Median Absolute Deviations (MADs)
#' above the reference median to use as the anomaly cutoff. Default is 2.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#'
#' @return A list containing the following components for each cell type and the combined data:
#' \item{reference_reconstruction_errors}{Reconstruction error (SSE) for each cell in the reference data.}
#' \item{reference_anomaly}{Logical vector indicating whether each reference cell is classified as an anomaly.}
#' \item{query_reconstruction_errors}{Reconstruction error (SSE) for each cell in the query data (if provided).}
#' \item{query_anomaly}{Logical vector indicating whether each query cell is classified as an anomaly.}
#' \item{applied_threshold}{The numeric threshold applied to determine anomalies for that cell type.}
#' \item{var_explained}{Proportion of variance explained by the retained principal components for that cell type's local PCA.}
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateReconstructionErrorObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Calculate PCA reconstruction errors
#' recon_output <- calculateReconstructionError(
#'     reference_data = reference_data,
#'     query_data = query_data,
#'     ref_cell_type_col = "expert_annotation",
#'     query_cell_type_col = "SingleR_annotation",
#'     pc_subset = 1:5,
#'     n_hvgs = 100,
#'     mad_multiplier = 2
#' )
#'
#' # Plot the output for a specific cell type
#' plot(recon_output,
#'      cell_type = "CD4",
#'      plot_type = "violin")
#'
#' @importFrom methods is
#' @importFrom stats na.omit median mad prcomp
#'
# Function to calculate PCA Reconstruction Errors
calculateReconstructionError <- function(reference_data,
                                         query_data = NULL,
                                         ref_cell_type_col,
                                         query_cell_type_col = NULL,
                                         cell_types = NULL,
                                         pc_subset = 1:5,
                                         n_hvgs = 100,
                                         mad_multiplier = 2,
                                         assay_name = "logcounts",
                                         max_cells_query = 5000,
                                         max_cells_ref = 5000) {

    # Reconstruction Error strictly requires PCA
    if (is.null(pc_subset)) {
        stop("'pc_subset' cannot be NULL for Reconstruction Error calculations. PCA is required.")
    }

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  pc_subset_ref = NULL,
                  assay_name = assay_name,
                  max_cells_query = max_cells_query,
                  max_cells_ref = max_cells_ref)

    # Check if n_hvgs is a positive integer
    if (!is.numeric(n_hvgs) || length(n_hvgs) != 1 || n_hvgs <= 0 || n_hvgs != as.integer(n_hvgs)) {
        stop("\'n_hvgs\' must be a single positive integer.")
    }
    n_hvgs <- as.integer(n_hvgs)

    # Check for scran dependency (needed for local HVGs)
    if (!requireNamespace("scran", quietly = TRUE)) {
        stop("Package 'scran' is required to calculate local highly variable genes. Please install it.")
    }

    # Convert cell type columns to character if needed
    reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                convert_cols = ref_cell_type_col)
    if(!is.null(query_data)){
        query_data <- convertColumnsToCharacter(sce_object = query_data,
                                                convert_cols = query_cell_type_col)
    }

    # Input check for mad_multiplier
    if (!is.numeric(mad_multiplier) || mad_multiplier <= 0) {
        stop("\'mad_multiplier\' must be a positive numeric value.")
    }

    # Select and validate cell types
    cell_types <- selectCellTypes(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = query_cell_type_col,
                                  ref_cell_type_col = ref_cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = NULL)

    # Downsample datasets
    reference_data <- downsampleSCE(sce_object = reference_data,
                                    max_cells = max_cells_ref,
                                    cell_types = cell_types,
                                    cell_type_col = ref_cell_type_col)

    if(!is.null(query_data)){
        query_data <- downsampleSCE(sce_object = query_data,
                                    max_cells = max_cells_query,
                                    cell_types = cell_types,
                                    cell_type_col = query_cell_type_col)
    }

    # Extract metadata vectors
    ref_cell_types <- reference_data[[ref_cell_type_col]]
    if(!is.null(query_data)){
        query_cell_types <- query_data[[query_cell_type_col]]
    }

    # ________________________________________________________________
    # CORE MATH: Local PCA Computation & Reconstruction per Cell Type
    # ________________________________________________________________

    output <- list()
    cell_types_list <- as.list(cell_types)
    cell_types_list[["Combined"]] <- cell_types

    for (cell_type in cell_types_list) {

        list_name <- ifelse(length(cell_type) == 1, cell_type, "Combined")

        # 1. Subset Reference Cells
        ref_subset_idx <- which(ref_cell_types %in% cell_type)

        # Check if enough cells exist to perform PCA
        if (length(ref_subset_idx) < max(pc_subset) + 1) {
            warning(paste("Not enough reference cells in", list_name, "to compute PCA. Skipping."))
            next
        }

        ref_sce_sub <- reference_data[, ref_subset_idx]

        # 2. Find Local HVGs for this cell type
        var_stats <- scran::modelGeneVar(ref_sce_sub, assay.type = assay_name)
        n_hvgs_actual <- min(n_hvgs, nrow(ref_sce_sub))
        local_hvgs <- scran::getTopHVGs(var_stats, n = n_hvgs_actual)

        # Ensure genes exist in query
        if (!is.null(query_data) && !all(local_hvgs %in% rownames(SummarizedExperiment::assay(query_data, assay_name)))) {
            warning(paste("Some local HVGs for", list_name, "are missing from the query data. Skipping."))
            next
        }

        # 3. Extract and Center Local Reference Matrix
        ref_assay <- as.matrix(SummarizedExperiment::assay(ref_sce_sub, assay_name)[local_hvgs, , drop = FALSE])
        centering_vec <- Matrix::rowMeans(ref_assay)

        ref_transposed <- t(ref_assay)
        ref_centered <- sweep(ref_transposed, 2, centering_vec, "-")

        # 4. Compute Local PCA (Using highly optimized base prcomp)
        pca_res <- stats::prcomp(ref_centered, center = FALSE, scale. = FALSE)

        # Adjust pc_subset if the cell type has fewer available PCs than requested
        max_pc_avail <- ncol(pca_res$rotation)
        current_pc_subset <- pc_subset[pc_subset <= max_pc_avail]

        rotation_mat <- pca_res$rotation[, current_pc_subset, drop = FALSE]

        # Calculate local variance explained
        percent_var <- (pca_res$sdev[current_pc_subset]^2 / sum(pca_res$sdev^2)) * 100

        # 5. Project and Reconstruct Reference
        ref_pca_scores <- ref_centered %*% rotation_mat
        ref_reconstructed <- ref_pca_scores %*% t(rotation_mat)

        ref_residuals <- ref_centered - ref_reconstructed
        ref_errors_subset <- rowSums(ref_residuals^2)
        names(ref_errors_subset) <- colnames(ref_sce_sub)

        # 6. Calculate dynamic threshold based on Reference subset (MAD only)
        ref_median <- stats::median(ref_errors_subset, na.rm = TRUE)
        ref_mad <- stats::mad(ref_errors_subset, na.rm = TRUE)
        cutoff <- ref_median + (mad_multiplier * ref_mad)

        # 7. Project and Reconstruct Query (if available)
        if (!is.null(query_data)) {
            query_subset_idx <- which(query_cell_types %in% cell_type)

            if (length(query_subset_idx) > 0) {
                query_sce_sub <- query_data[, query_subset_idx]
                query_assay <- as.matrix(SummarizedExperiment::assay(query_sce_sub, assay_name)[local_hvgs, , drop = FALSE])

                query_transposed <- t(query_assay)
                query_centered <- sweep(query_transposed, 2, centering_vec, "-")

                query_pca_scores <- query_centered %*% rotation_mat
                query_reconstructed <- query_pca_scores %*% t(rotation_mat)

                query_residuals <- query_centered - query_reconstructed
                query_errors_subset <- rowSums(query_residuals^2)
                names(query_errors_subset) <- colnames(query_sce_sub)
            } else {
                query_errors_subset <- numeric(0)
            }
        }

        # 8. Store results
        output[[list_name]] <- list()
        output[[list_name]][["reference_reconstruction_errors"]] <- ref_errors_subset
        output[[list_name]][["reference_anomaly"]] <- ref_errors_subset > cutoff

        output[[list_name]][["reference_mat_subset"]] <- t(ref_assay)

        if (!is.null(query_data)) {
            output[[list_name]][["query_reconstruction_errors"]] <- query_errors_subset
            output[[list_name]][["query_anomaly"]] <- query_errors_subset > cutoff

            output[[list_name]][["query_mat_subset"]] <- t(query_assay)
        }

        output[[list_name]][["applied_threshold"]] <- cutoff
        output[[list_name]][["var_explained"]] <- percent_var
    }

    if (length(output) == 0) {
        warning("No cell types had sufficient cells to compute Reconstruction Errors.")
    }

    # Set the S3 class for plotting downstream
    class(output) <- c(class(output), "calculateReconstructionErrorObject")

    return(output)
}
