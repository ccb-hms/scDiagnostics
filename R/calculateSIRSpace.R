#' @title Calculate Sliced Inverse Regression (SIR) Space for Different Cell Types
#'
#' @description
#' This function calculates the SIR space projections for different cell types in the query and reference datasets.
#'
#' @details
#' The function projects the query dataset onto the SIR space of the reference dataset based on shared cell types.
#' It computes conditional means for the reference dataset, extracts the SVD components, and performs the projection
#' of both the query and reference data. It uses the `projectSIR` function to perform the actual projection and
#' allows the user to specify particular cell types for analysis.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing the numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing the numeric expression matrix for the reference cells.
#' @param query_cell_type_col A character string specifying the column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col A character string specifying the column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the analysis. If NULL, all common cell types between the query and reference data will be used.
#' @param multiple_cond_means Logical. Whether to compute conditional means for multiple conditions in the reference dataset. Default is TRUE.
#' @param cumulative_variance_threshold A numeric value specifying the cumulative variance threshold for selecting principal components. Default is 0.7.
#' @param n_neighbor A numeric value specifying the number of neighbors for computing the SIR space. Default is 1.
#' @param assay_name A character string specifying the name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#'
#' @return A list containing the SIR projections, rotation matrix, and percentage of variance explained for the given cell types.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateSIRSpaceObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Compute important variables for all pairwise cell comparisons
#' sir_output <- calculateSIRSpace(reference_data = reference_data,
#'                                 query_data = query_data,
#'                                 query_cell_type_col = "expert_annotation",
#'                                 ref_cell_type_col = "expert_annotation",
#'                                 multiple_cond_means = TRUE,
#'                                 cumulative_variance_threshold = 0.9,
#'                                 n_neighbor = 1)
#'
#' # Generate plots SIR projections
#' plot(sir_output,
#'      sir_subset = 1:5,
#'      cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'      lower_facet = "scatter",
#'      diagonal_facet = "boxplot",
#'      upper_facet = "blank")
#'
#' # Plot top loadings
#' plot(sir_output,
#'      sir_subset = 1:5,
#'      plot_type = "loadings",
#'      n_top = 10)
#'
# Function to plot cell types in SIR space
calculateSIRSpace <- function(query_data,
                              reference_data,
                              query_cell_type_col,
                              ref_cell_type_col,
                              cell_types = NULL,
                              multiple_cond_means = TRUE,
                              cumulative_variance_threshold = 0.7,
                              n_neighbor = 1,
                              assay_name = "logcounts",
                              max_cells_query = 5000,
                              max_cells_ref = 5000){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  assay_name = assay_name,
                  max_cells_query = max_cells_query,
                  max_cells_ref = max_cells_ref)

    # Check if cumulative_variance_threshold is between 0 and 1
    if (!is.numeric(cumulative_variance_threshold) ||
        cumulative_variance_threshold < 0 || cumulative_variance_threshold > 1) {
        stop("cumulative_variance_threshold must be a numeric value between 0 and 1.")
    }

    # Check if n_neighbor is a positive integer
    if (!is.numeric(n_neighbor) || n_neighbor <= 0 ||
        n_neighbor != as.integer(n_neighbor)) {
        stop("n_neighbor must be a positive integer.")
    }

    # Select cell types
    cell_types <- selectCellTypes(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = query_cell_type_col,
                                  ref_cell_type_col = ref_cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = NULL)

    # Convert cell type columns to character if needed
    query_data <- convertColumnsToCharacter(sce_object = query_data,
                                            convert_cols = query_cell_type_col)
    reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                convert_cols = ref_cell_type_col)

    # Downsample query and reference data (with cell type filtering)
    query_data <- downsampleSCE(sce_object = query_data,
                                cell_types = cell_types,
                                cell_type_col = query_cell_type_col,
                                max_cells = max_cells_query)
    reference_data <- downsampleSCE(sce_object = reference_data,
                                    cell_types = cell_types,
                                    cell_type_col = ref_cell_type_col,
                                    max_cells = max_cells_ref)

    # Get the projected PCA data
    sir_output <- projectSIR(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             cell_types = cell_types,
                             multiple_cond_means = multiple_cond_means,
                             assay_name = assay_name,
                             cumulative_variance_threshold = cumulative_variance_threshold,
                             n_neighbor = n_neighbor)

    # Return SIR projections output
    class(sir_output) <- c(class(sir_output), "calculateSIRSpaceObject")
    return(sir_output)
}


