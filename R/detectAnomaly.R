#' @title PCA Anomaly Scores via Isolation Forests with Visualization
#'
#' @description
#' This function detects anomalies in single-cell data by projecting the data onto a PCA space and using an isolation forest
#' algorithm to identify anomalies.
#'
#' @details
#' This function projects the query data onto the PCA space of the reference data. An isolation forest is then built on the
#' reference data to identify anomalies in the query data based on their PCA projections. If no query dataset is provided by the user,
#' the anomaly scores are computed on the reference data itself. Anomaly scores for the data with all combined cell types are also
#' provided as part of the output.
#'
#' @param reference_data A \code{\link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_data An optional \code{\link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' If NULL, then the isolation forest anomaly scores are computed for the reference data. Default is NULL.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to use in the analysis. Default is 1:5
#' If set to \code{NULL} then no dimensionality reduction is performed and the assay data is used directly for computations.
#' @param n_hvgs An integer specifying the number of highly variable genes to retain when `pc_subset` is NULL.
#' If a query dataset is provided, the top `n_hvgs` are computed for both reference and query, and their union is used. Default is 1000.
#' @param n_tree An integer specifying the number of trees for the isolation forest. Default is 500
#' @param threshold_method A character string specifying the method to determine anomaly cutoffs.
#' Options are \code{"MAD"} (Median Absolute Deviation) or \code{"absolute"}. Default is \code{"MAD"}.
#' @param mad_multiplier A numeric value specifying the number of MADs above the reference median to use as the cutoff
#' when \code{threshold_method = "MAD"}. Default is 2.
#' @param anomaly_threshold A numeric value specifying the absolute threshold for identifying anomalies
#' when \code{threshold_method = "absolute"}. Default is 0.5.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#' @param ... Additional arguments passed to the `isolation.forest` function.
#'
#' @return A list containing the following components for each cell type and the combined data:
#' \item{anomaly_scores}{Anomaly scores for each cell in the query data.}
#' \item{anomaly}{Logical vector indicating whether each cell is classified as an anomaly.}
#' \item{reference_mat_subset}{PCA projections of the reference data.}
#' \item{query_mat_subset}{PCA projections of the query data (if provided).}
#' \item{var_explained}{Proportion of variance explained by the retained principal components.}
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.detectAnomalyObject}}
#'
#' @references
#' \itemize{
#'   \item Liu, F. T., Ting, K. M., & Zhou, Z. H. (2008). Isolation forest. In 2008 Eighth IEEE International Conference on Data Mining (pp. 413-422). IEEE.
#'   \item \href{https://cran.r-project.org/web/packages/isotree/isotree.pdf}{isotree: Isolation-Based Outlier Detection}
#' }
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Store PCA anomaly data
#' anomaly_output <- detectAnomaly(reference_data = reference_data,
#'                                 query_data = query_data,
#'                                 ref_cell_type_col = "expert_annotation",
#'                                 query_cell_type_col = "SingleR_annotation",
#'                                 pc_subset = 1:5,
#'                                 n_tree = 500,
#'                                 threshold_method = "MAD",
#'                                 mad_multiplier = 2)
#'
#' # Plot the output for a cell type
#' plot(anomaly_output,
#'      cell_type = "CD4",
#'      pc_subset = 1:3,
#'      data_type = "query")
#'
#' @importFrom methods is
#' @importFrom stats na.omit predict qnorm median mad
#' @importFrom utils tail
#'
# Function to perform diagnostics using isolation forest with PCA and visualization
detectAnomaly <- function(reference_data,
                          query_data = NULL,
                          ref_cell_type_col,
                          query_cell_type_col = NULL,
                          cell_types = NULL,
                          pc_subset = 1:5,
                          n_hvgs = 1000,
                          n_tree = 500,
                          threshold_method = c("MAD", "absolute"),
                          mad_multiplier = 2,
                          anomaly_threshold = 0.5,
                          assay_name = "logcounts",
                          max_cells_query = 5000,
                          max_cells_ref = 5000,
                          ...) {

    # Match arguments
    threshold_method <- match.arg(threshold_method)

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name,
                  max_cells_query = max_cells_query,
                  max_cells_ref = max_cells_ref)

    # Convert cell type columns to character if needed
    reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                convert_cols = ref_cell_type_col)
    if(!is.null(query_data)){
        query_data <- convertColumnsToCharacter(sce_object = query_data,
                                                convert_cols = query_cell_type_col)
    }

    # Check if n_hvgs is a positive integer
    if (is.null(pc_subset)) {
        if (!is.numeric(n_hvgs) || length(n_hvgs) != 1 || n_hvgs <= 0 || n_hvgs != as.integer(n_hvgs)) {
            stop("\'n_hvgs\' must be a single positive integer.")
        }
        n_hvgs <- as.integer(n_hvgs)
    }

    # Check if n_tree is a positive integer
    if (!is.numeric(n_tree) || n_tree <= 0 || n_tree != as.integer(n_tree)) {
        stop("\'n_tree\' must be a positive integer.")
    }

    # Input check for mad_multiplier
    if (!is.numeric(mad_multiplier) || mad_multiplier <= 0) {
        stop("\'mad_multiplier\' must be a positive numeric value.")
    }

    # Input check for anomaly_threshold
    if (!is.numeric(anomaly_threshold) || anomaly_threshold <= 0 ||
        anomaly_threshold >= 1) {
        stop("\'anomaly_threshold\' must be a positive number greater than 0 and less than 1.")
    }

    # Select cell types
    cell_types <- selectCellTypes(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = query_cell_type_col,
                                  ref_cell_type_col = ref_cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = NULL)

    # Get data from reference and query datasets
    if(!is.null(pc_subset)){
        if(!is.null(query_data)){
            pca_output <- projectPCA(query_data = query_data,
                                     reference_data = reference_data,
                                     query_cell_type_col = query_cell_type_col,
                                     ref_cell_type_col = ref_cell_type_col,
                                     cell_types = cell_types,
                                     pc_subset = pc_subset,
                                     assay_name = assay_name,
                                     max_cells_ref = max_cells_ref,
                                     max_cells_query = max_cells_query)
            query_mat <- pca_output[pca_output[["dataset"]] == "Query",]
            reference_mat <- pca_output[pca_output[["dataset"]] == "Reference",]

            # Extract cell type information from PCA output
            query_cell_types <- query_mat[["cell_type"]]
            reference_cell_types <- reference_mat[["cell_type"]]

            # Remove the metadata columns to keep only PCA coordinates
            pc_cols <- grep("^PC", colnames(pca_output), value = TRUE)
            query_mat <- query_mat[, pc_cols, drop = FALSE]
            reference_mat <- reference_mat[, pc_cols, drop = FALSE]
        } else {
            reference_data <- downsampleSCE(sce_object = reference_data,
                                            max_cells = max_cells_ref,
                                            cell_types =  cell_types,
                                            cell_type_col = ref_cell_type_col)
            reference_mat <- reducedDim(reference_data, "PCA")[, pc_subset]
            reference_cell_types <- reference_data[[ref_cell_type_col]]
        }
    } else {
        # PCA is NULL: Use HVGs to avoid the Curse of Dimensionality in Isolation Forests

        reference_data <- downsampleSCE(sce_object = reference_data,
                                        max_cells = max_cells_ref,
                                        cell_types =  cell_types,
                                        cell_type_col = ref_cell_type_col)

        # Check for scran dependency
        if (!requireNamespace("scran", quietly = TRUE)) {
            stop("Package 'scran' is required when pc_subset = NULL. Please install it with: BiocManager::install('scran')")
        }

        # Get Reference HVGs
        var_ref <- scran::modelGeneVar(reference_data, assay.type = assay_name)
        hvg_ref <- scran::getTopHVGs(var_ref, n = n_hvgs)

        hvg_combined <- hvg_ref

        if(!is.null(query_data)){
            # Downsample query
            query_data <- downsampleSCE(sce_object = query_data,
                                        max_cells = max_cells_query,
                                        cell_types =  cell_types,
                                        cell_type_col = query_cell_type_col)

            # Get Query HVGs
            var_query <- scran::modelGeneVar(query_data, assay.type = assay_name)
            hvg_query <- scran::getTopHVGs(var_query, n = n_hvgs)

            # Combine and take unique (Union of Ref and Query HVGs)
            hvg_combined <- unique(c(hvg_ref, hvg_query))
        }

        # Subset matrices to the combined HVGs
        reference_mat <- t(as.matrix(assay(reference_data, assay_name)[hvg_combined, ]))
        reference_cell_types <- reference_data[[ref_cell_type_col]]

        if(!is.null(query_data)){
            query_mat <- t(as.matrix(assay(query_data, assay_name)[hvg_combined, ]))
            query_cell_types <- query_data[[query_cell_type_col]]
        }
    }

    # List to store output
    output <- list()

    # Cell types list for anomaly detection
    cell_types_list <- as.list(cell_types)
    cell_types_list[["Combined"]] <- cell_types

    for (cell_type in cell_types_list) {

        # Filter reference and query PCA data for the current cell type
        reference_mat_subset <- na.omit(reference_mat[which(
            reference_cell_types %in% cell_type),])

        # Build isolation forest on reference PCA data for this cell type
        isolation_forest <- isotree::isolation.forest(reference_mat_subset,
                                                      ntree = n_tree,)

        # Calculate anomaly scores for query data (scaled by reference path length)
        reference_anomaly_scores <- predict(isolation_forest,
                                            newdata = reference_mat_subset,
                                            type = "score")

        # --- NEW: Calculate Dynamic Cutoff based on Reference Scores ---
        if (threshold_method == "MAD") {
            ref_median <- median(reference_anomaly_scores)
            ref_mad <- mad(reference_anomaly_scores)
            cutoff <- ref_median + (mad_multiplier * ref_mad)
        } else {
            cutoff <- anomaly_threshold
        }

        if(!is.null(query_data)){
            query_mat_subset <- na.omit(query_mat[which(
                query_cell_types %in% cell_type),])
            query_anomaly_scores <- predict(isolation_forest,
                                            newdata = query_mat_subset,
                                            type = "score")
        }

        # Store cell type anomaly scores and PCA data
        list_name <- ifelse(length(cell_type) == 1, cell_type, "Combined")
        output[[list_name]] <- list()
        output[[list_name]][["reference_anomaly_scores"]] <- reference_anomaly_scores
        output[[list_name]][["reference_anomaly"]] <- reference_anomaly_scores > cutoff
        output[[list_name]][["reference_mat_subset"]] <- reference_mat_subset

        if(!is.null(query_data)){
            output[[list_name]][["query_mat_subset"]] <- query_mat_subset
            output[[list_name]][["query_anomaly_scores"]] <- query_anomaly_scores
            output[[list_name]][["query_anomaly"]] <- query_anomaly_scores > cutoff
        }

        if(!is.null(pc_subset)){
            output[[list_name]][["var_explained"]] <- attributes(reducedDim(reference_data, "PCA"))[["percentVar"]][pc_subset]
        }

        # Optionally, store the threshold used so the user (or plotting function) knows what was applied
        output[[list_name]][["applied_threshold"]] <- cutoff
    }

    # Set the class of the output
    class(output) <- c(class(output), "detectAnomalyObject")

    # Return anomaly, PCA data and optional PCA anomaly plots for each cell type
    return(output)
}
