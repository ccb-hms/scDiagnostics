#' @title Compute Cell Distances Between Reference and Query Data
#'
#' @description
#' This function computes the distances within the reference dataset and the distances from each query cell to all
#' reference cells for each cell type. It uses PCA for dimensionality reduction and Euclidean distance for distance calculation.
#'
#' @details
#' The function first performs PCA on the reference dataset and projects the query dataset onto the same PCA space.
#' It then computes pairwise Euclidean distances within the reference dataset for each cell type, as well as distances from each
#' query cell to all reference cells of a particular cell type. The results are stored in a list, with one entry per cell type.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data}
#' that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data}
#' that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default 1:5.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#'
#' @return A list containing distance data for each cell type. Each entry in the list contains:
#' \describe{
#'   \item{ref_distances}{A vector of all pairwise distances within the reference subset for the cell type.}
#'   \item{query_to_ref_distances}{A matrix of distances from each query cell to all reference cells for the cell type.}
#' }
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateCellDistancesObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Plot the PC data
#' distance_data <- calculateCellDistances(query_data = query_data,
#'                                         reference_data = reference_data,
#'                                         query_cell_type_col = "SingleR_annotation",
#'                                         ref_cell_type_col = "expert_annotation",
#'                                         pc_subset = 1:10)
#'
#' # Identify outliers for CD4
#' cd4_anomalies <- detectAnomaly(reference_data = reference_data,
#'                                query_data = query_data,
#'                                query_cell_type_col = "SingleR_annotation",
#'                                ref_cell_type_col = "expert_annotation",
#'                                pc_subset = 1:10,
#'                                n_tree = 500,
#'                                anomaly_threshold = 0.5)
#' cd4_top6_anomalies <- names(sort(cd4_anomalies$CD4$query_anomaly_scores, decreasing = TRUE)[1:6])
#'
#' # Plot the densities of the distances
#' plot(distance_data, ref_cell_type = "CD4", cell_names = cd4_top6_anomalies)
#' plot(distance_data, ref_cell_type = "CD8", cell_names = cd4_top6_anomalies)
#'
# Function to compute distances within reference data and between query data and reference cells
calculateCellDistances <- function(query_data,
                                   reference_data,
                                   query_cell_type_col,
                                   ref_cell_type_col,
                                   cell_types = NULL,
                                   pc_subset = 1:5,
                                   assay_name = "logcounts") {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             pc_subset = pc_subset,
                             assay_name = assay_name)

    # Create a list to store distance data for each cell type
    distance_data <- vector("list", length = length(cell_types))
    names(distance_data) <- cell_types

    # Function to compute Euclidean distance between a vector and each row of a matrix
    .compute_distances <- function(matrix, vector) {

        # Apply the distance function to each row of the matrix
        distances <- apply(matrix, 1, function(row) {
            sqrt(sum((row - vector) ^ 2))
        })

        return(distances)
    }

    for (cell_type in cell_types) {

        # Subset principal component scores for current cell type
        ref_subset_scores <- pca_output[which(
            pca_output[["dataset"]] == "Reference" &
                pca_output[["cell_type"]] == cell_type), pc_subset]
        query_subset_scores <- pca_output[pca_output[["dataset"]] == "Query",
                                          pc_subset]

        # Compute all pairwise distances within the reference subset
        ref_distances <- as.vector(dist(ref_subset_scores))

        # Compute distances from each query cell to all reference cells
        query_to_ref_distances <- apply(
            query_subset_scores, 1, function(query_cell,
                                             ref_subset_scores) {
                .compute_distances(ref_subset_scores,
                                   query_cell)
                }, ref_subset_scores = ref_subset_scores)

        # Store the distances
        distance_data[[cell_type]] <- list(
            ref_distances = ref_distances,
            query_to_ref_distances = t(query_to_ref_distances)
        )
    }

    # Add class of object
    class(distance_data) <- c(class(distance_data),
                              "calculateCellDistancesObject")

    # Return the distance data
    return(distance_data)
}
