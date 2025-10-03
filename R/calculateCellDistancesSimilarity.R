#' @title Function to Calculate Bhattacharyya Coefficients and Hellinger Distances
#'
#' @description
#' This function computes Bhattacharyya coefficients and Hellinger distances to quantify the similarity of density
#' distributions between query cells and reference data for each cell type.
#'
#' @details
#' This function first computes distance data using the \code{calculateCellDistances} function, which calculates
#' pairwise distances between cells within the reference data and between query cells and reference cells in the PCA space.
#' Bhattacharyya coefficients and Hellinger distances are calculated to quantify the similarity of density distributions between query
#' cells and reference data for each cell type. Bhattacharyya coefficient measures the similarity of two probability distributions,
#' while Hellinger distance measures the distance between two probability distributions.
#'
#' Bhattacharyya coefficients range between 0 and 1. A value closer to 1 indicates higher similarity between distributions, while a value
#' closer to 0 indicates lower similarity
#'
#' Hellinger distances range between 0 and 1. A value closer to 0 indicates higher similarity between distributions, while a value
#' closer to 1 indicates lower similarity.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data}
#' that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data}
#' that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are
#'                   included.
#' @param cell_names_query A character vector specifying the names of the query cells for which to compute distance measures.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is 1:5.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
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
#' # Get overlap measures
#' overlap_measures <- calculateCellDistancesSimilarity(query_data = query_data,
#'                                                      reference_data = reference_data,
#'                                                      cell_names_query = cd4_top6_anomalies,
#'                                                      query_cell_type_col = "SingleR_annotation",
#'                                                      ref_cell_type_col = "expert_annotation",
#'                                                      pc_subset = 1:10)
#' overlap_measures
#'
# Function to compute Bhattacharyya coefficients and Hellinger distances
calculateCellDistancesSimilarity <- function(query_data,
                                             reference_data,
                                             query_cell_type_col,
                                             ref_cell_type_col,
                                             cell_types = NULL,
                                             cell_names_query,
                                             pc_subset = 1:5,
                                             assay_name = "logcounts",
                                             max_cells_ref = 5000) {

    # Format the query cell names - remove "Query_" prefix if present
    cell_names_query <- gsub("^Query_", "", cell_names_query)

    # Check standard input arguments (now with proper cell_types)
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_names_query = cell_names_query,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name,
                  max_cells_ref = max_cells_ref)

    # Convert cell type columns to character if needed
    query_data <- convertColumnsToCharacter(sce_object = query_data,
                                            convert_cols = query_cell_type_col)
    reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                convert_cols = ref_cell_type_col)

    # Select cell types
    cell_types <- selectCellTypes(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = query_cell_type_col,
                                  ref_cell_type_col = ref_cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = NULL)

    # Subset query data first
    query_data_subset <- query_data[, cell_names_query, drop = FALSE]

    # Compute distance data
    distance_data <- calculateCellDistances(
        query_data = query_data_subset,
        reference_data = reference_data,
        query_cell_type_col = query_cell_type_col,
        ref_cell_type_col = ref_cell_type_col,
        cell_types = cell_types,
        pc_subset = pc_subset,
        assay_name = assay_name,
        max_cells_ref = max_cells_ref)

    # Update the names of the query cells
    cell_names_query <- paste0("Query_", cell_names_query)

    # Initialize empty lists to store results
    bhattacharyya_list <- hellinger_list <-
        vector("list", length = length(distance_data))
    names(bhattacharyya_list) <- names(hellinger_list) <- names(distance_data)

    # Iterate over each cell type
    for (cell_type in names(distance_data)) {

        # Extract distances within the reference dataset for the current cell type
        ref_distances <- distance_data[[cell_type]][["ref_distances"]]

        # Compute density of reference distances
        ref_density <- density(ref_distances)

        # Initialize an empty vector to store overlap measures for the current cell type
        bhattacharyya_coef <- numeric(length(cell_names_query))
        hellinger_dist <- numeric(length(cell_names_query))

        # Iterate over each cell
        for (i in seq_len(length(cell_names_query))) {

            # Extract distances from the current cell to reference cells
            cell_distances <-
                distance_data[[cell_type]][["query_to_ref_distances"]][cell_names_query[i], , drop = FALSE]

            # Compute density of cell distances
            cell_density <- density(cell_distances)

            # Create a common grid for evaluating densities
            common_grid <- seq(min(min(ref_density[["x"]]),
                                   min(cell_density[["x"]]), 0),
                               max(max(ref_density[["x"]]),
                                   max(cell_density[["x"]])),
                               length.out = 1000)

            # Interpolate densities onto the common grid
            ref_density_interp <- approxfun(ref_density[["x"]],
                                            ref_density[["y"]])(common_grid)
            ref_density_interp[is.na(ref_density_interp)] <- 0
            cell_density_interp <- approxfun(cell_density[["x"]],
                                             cell_density[["y"]])(common_grid)
            cell_density_interp[is.na(cell_density_interp)] <- 0

            # Compute and store Bhattacharyya coefficient/Hellinger distance
            bhattacharyya_coef[i] <- sum(
                sqrt(ref_density_interp * cell_density_interp) *
                    mean(diff(common_grid)))
            hellinger_dist[i] <- sqrt(
                1 - sum(sqrt(ref_density_interp * cell_density_interp)) *
                    mean(diff(common_grid)))
        }

        # Store overlap measures for the current cell type
        bhattacharyya_list[[cell_type]] <- bhattacharyya_coef
        hellinger_list[[cell_type]] <- hellinger_dist
    }

    # Return list with overlap measures
    bhattacharyya_coef <- data.frame(Cell = cell_names_query, bhattacharyya_list)
    hellinger_dist <- data.frame(Cell = cell_names_query, hellinger_list)
    return(list(bhattacharyya_coef = bhattacharyya_coef,
                hellinger_dist = hellinger_dist))
}


