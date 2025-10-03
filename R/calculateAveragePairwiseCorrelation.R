#' @title Compute Average Pairwise Correlation between Cell Types
#'
#' @description
#' Computes the average pairwise correlations between specified cell types
#' in single-cell gene expression data.
#'
#' @details
#' This function operates on \code{\linkS4class{SingleCellExperiment}} objects,
#' ideal for single-cell analysis workflows. It calculates pairwise correlations between query and
#' reference cells using a specified correlation method, then averages these correlations for each
#' cell type pair. This function aids in assessing the similarity between cells in reference and query datasets,
#' providing insights into the reliability of cell type annotations in single-cell gene expression data.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to use in the analysis. Default is 1:10.
#' If set to \code{NULL} then no dimensionality reduction is performed and the assay data is used directly for computations.
#' @param correlation_method The correlation method to use for calculating pairwise correlations.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000
#'
#' @return A matrix containing the average pairwise correlation values.
#'         Rows and columns are labeled with the cell types. Each element
#'         in the matrix represents the average correlation between a pair
#'         of cell types.
#'
#' @seealso \code{\link{plot.calculateAveragePairwiseCorrelationObject}}
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Compute pairwise correlations
#' cor_matrix_avg <- calculateAveragePairwiseCorrelation(query_data = query_data,
#'                                                       reference_data = reference_data,
#'                                                       query_cell_type_col = "SingleR_annotation",
#'                                                       ref_cell_type_col = "expert_annotation",
#'                                                       cell_types = c("CD4", "CD8", "B_and_plasma"),
#'                                                       pc_subset = 1:10,
#'                                                       correlation_method = "spearman")
#'
#' # Visualize correlation output
#' plot(cor_matrix_avg)
#'
#' @import SingleCellExperiment
#' @importFrom stats cor
#'
#' @export
#'
# Function to calculate average pairwise correlation between cell types
calculateAveragePairwiseCorrelation <- function(
        query_data,
        reference_data,
        query_cell_type_col,
        ref_cell_type_col,
        cell_types = NULL,
        pc_subset = 1:10,
        correlation_method = c("spearman", "pearson"),
        assay_name = "logcounts",
        max_cells_query = 5000,
        max_cells_ref = 5000) {

    # Match correlation method argument
    correlation_method <- match.arg(correlation_method)

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
                                  dual_only = TRUE,
                                  n_cell_types = NULL)

    # Function to compute correlation between two cell types
    .computeCorrelation <- function(type1, type2) {

        if(!is.null(pc_subset)){
            # Project query data onto PCA space of reference data
            pca_output <- projectPCA(
                query_data = query_data,
                reference_data = reference_data,
                query_cell_type_col = query_cell_type_col,
                ref_cell_type_col = ref_cell_type_col,
                cell_types = cell_types,
                pc_subset = pc_subset,
                assay_name = assay_name,
                max_cells_ref = max_cells_ref,
                max_cells_query = max_cells_query)
            ref_mat <- pca_output[which(
                pca_output[["dataset"]] == "Reference" &
                    pca_output[["cell_type"]] == type2),
                paste0("PC", pc_subset)]
            query_mat <- pca_output[which(
                pca_output[["dataset"]] == "Query" &
                    pca_output[["cell_type"]] == type1),
                paste0("PC", pc_subset)]
        } else{

            # Subset query and reference data to the specified cell type
            query_subset <- query_data[, which(
                query_data[[query_cell_type_col]] == type1), drop = FALSE]
            ref_subset <- reference_data[, which(
                reference_data[[ref_cell_type_col]] == type2), drop = FALSE]

            query_mat <- t(as.matrix(assay(query_subset, assay_name)))
            ref_mat <- t(as.matrix(assay(ref_subset, assay_name)))
        }

        cor_matrix <- cor(t(query_mat), t(ref_mat),
                          method = correlation_method)
        return(mean(cor_matrix))
    }

    # Use outer to compute pairwise correlations
    cor_matrix_avg <- outer(cell_types, cell_types,
                            Vectorize(.computeCorrelation))

    # Assign cell type names to rows and columns
    rownames(cor_matrix_avg) <- paste0("Query-", cell_types)
    colnames(cor_matrix_avg) <- paste0("Ref-", cell_types)

    # Update class of output
    class(cor_matrix_avg) <- c(class(cor_matrix_avg),
                               "calculateAveragePairwiseCorrelationObject")

    return(cor_matrix_avg)
}
