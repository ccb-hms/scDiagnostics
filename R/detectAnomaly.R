#' @title PCA Anomaly Scores via Isolation Forests with Visualization
#'
#' @description
#' This function detects anomalies in single-cell data by projecting the data onto a PCA space and using an isolation forest
#' algorithm to identify anomalies.
#'
#' @details This function projects the query data onto the PCA space of the reference data. An isolation forest is then built on the
#' reference data to identify anomalies in the query data based on their PCA projections. If no query dataset is provided by the user,
#' the anomaly scores are computed on the reference data itself. Anomaly scores for the data with all combined cell types are also
#' provided as part of the output.
#'
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_data An optional \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' If NULL, then the isolation forest anomaly scores are computed for the reference data. Default is NULL.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to use in the analysis. Default is 1:5
#' If set to \code{NULL} then no dimensionality reduction is performed and the assay data is used directly for computations.
#' @param n_tree An integer specifying the number of trees for the isolation forest. Default is 500
#' @param anomaly_threshold A numeric value specifying the threshold for identifying anomalies, Default is 0.6.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
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
#'                                 anomaly_threshold = 0.6)
#'
#' # Plot the output for a cell type
#' plot(anomaly_output,
#'      cell_type = "CD4",
#'      pc_subset = 1:3,
#'      data_type = "query")
#'
#' @importFrom methods is
#' @importFrom stats na.omit predict qnorm
#' @importFrom utils tail
#'
# Function to perform diagnostics using isolation forest with PCA and visualization
detectAnomaly <- function(reference_data,
                          query_data = NULL,
                          ref_cell_type_col,
                          query_cell_type_col = NULL,
                          cell_types = NULL,
                          pc_subset = 1:5,
                          n_tree = 500,
                          anomaly_threshold = 0.6,
                          assay_name = "logcounts",
                          ...) {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Get common cell types if they are not specified by user
    reference_labels <- reference_data[[ref_cell_type_col]]
    if(is.null(cell_types)){
        if(is.null(query_data)){
            cell_types <- c(as.list(na.omit(unique(reference_labels))),
                            list(na.omit(unique(reference_labels))))
        } else{
            query_labels <- query_data[[query_cell_type_col]]
            cell_types <- c(as.list(na.omit(intersect(unique(reference_labels),
                                                      unique(query_labels)))),
                            list(na.omit(intersect(unique(reference_labels),
                                                   unique(query_labels)))))
        }
    }

    # Check if n_tree is a positive integer
    if (!is.numeric(n_tree) || n_tree <= 0 || n_tree != as.integer(n_tree)) {
        stop("\'n_tree\' must be a positive integer.")
    }

    # Input check for anomaly_threshold
    if (!is.numeric(anomaly_threshold) || anomaly_threshold <= 0 ||
        anomaly_threshold >= 1) {
        stop("\'anomaly_threshold\' must be a positive number greater than 0 and less than 1.")
    }

    # Get data from reference and query datasets
    if(!is.null(pc_subset)){
        reference_mat <- reducedDim(reference_data, "PCA")[, pc_subset]
        if(!is.null(query_data)){
            pca_output <- projectPCA(query_data = query_data,
                                     reference_data = reference_data,
                                     query_cell_type_col = query_cell_type_col,
                                     ref_cell_type_col = ref_cell_type_col,
                                     pc_subset = pc_subset,
                                     assay_name = assay_name)
            query_mat <- pca_output[pca_output[["dataset"]] == "Query",
                                    paste0("PC", pc_subset)]
        }
    } else{
        reference_mat <- t(as.matrix(assay(reference_data, assay_name)))
        if(!is.null(query_data)){
            query_mat <- t(as.matrix(assay(query_data, assay_name)))
        }
    }

    # Extract reference and query labels
    reference_labels <- reference_data[[ref_cell_type_col]]
    query_labels <- query_data[[query_cell_type_col]]

    # List to store output
    output <- list()

    for (cell_type in cell_types) {

        # Filter reference and query PCA data for the current cell type
        reference_mat_subset <- na.omit(reference_mat[which(
            reference_labels %in% cell_type),])

        # Build isolation forest on reference PCA data for this cell type
        isolation_forest <- isotree::isolation.forest(reference_mat_subset,
                                                      ntree = n_tree)

        # Calculate anomaly scores for query data (scaled by reference path length)
        reference_anomaly_scores <- predict(isolation_forest,
                                            newdata = reference_mat_subset,
                                            type = "score")
        if(!is.null(query_data)){
            query_mat_subset <- na.omit(query_mat[which(
                query_labels %in% cell_type),])
            query_anomaly_scores <- predict(isolation_forest,
                                            newdata = query_mat_subset,
                                            type = "score")
        }

        # Store cell type anomaly scores and PCA data
        list_name <- ifelse(length(cell_type) == 1, cell_type, "Combined")
        output[[list_name]] <- list()
        output[[list_name]][["reference_anomaly_scores"]] <-
            reference_anomaly_scores
        output[[list_name]][["reference_anomaly"]] <-
            reference_anomaly_scores > anomaly_threshold
        output[[list_name]][["reference_mat_subset"]] <- reference_mat_subset
        if(!is.null(query_data)){
            output[[list_name]][["query_mat_subset"]] <- query_mat_subset
            output[[list_name]][["query_anomaly_scores"]] <-
                query_anomaly_scores
            output[[list_name]][["query_anomaly"]] <-
                query_anomaly_scores > anomaly_threshold
        }
        if(!is.null(pc_subset))
            output[[list_name]][["var_explained"]] <-
            attributes(reducedDim(
                reference_data, "PCA"))[["percentVar"]][pc_subset]
    }

    # Set the class of the output
    class(output) <- c(class(output), "detectAnomalyObject")

    # Return anomaly, PCA data and optional PCA anomaly plots for each cell type
    return(output)
}
