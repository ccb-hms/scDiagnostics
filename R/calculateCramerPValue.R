#' @title Calculate Cramer Test P-Values for Two-Sample Comparison of Multivariate ECDFs
#'
#' @description
#' This function performs the Cramer test for comparing multivariate empirical cumulative distribution functions (ECDFs)
#' between two samples.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Projects the data into the PCA space.
#'   \item Subsets the data to the specified cell types and principal components.
#'   \item Performs the Cramer test for each cell type using the \code{cramer.test} function in the \code{cramer} package.
#' }
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#'
#' @return A named vector of p-values from the Cramer test for each cell type.
#'
#' @references Baringhaus, L., & Franz, C. (2004). "On a new multivariate two-sample test".
#' Journal of Multivariate Analysis, 88(1), 190-206.
#'
#' @export
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Plot the PC data (with query data)
#' cramer_test <- calculateCramerPValue(reference_data = reference_data,
#'                                      query_data = query_data,
#'                                      ref_cell_type_col = "expert_annotation",
#'                                      query_cell_type_col = "SingleR_annotation",
#'                                      cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                                      pc_subset = 1:5)
#' cramer_test
#'
# Function to perform Cramer test for two-sample comparison of multivariate ECDFs
calculateCramerPValue <- function(query_data,
                                  reference_data,
                                  query_cell_type_col,
                                  ref_cell_type_col,
                                  cell_types = NULL,
                                  pc_subset = 1:5,
                                  assay_name = "logcounts",
                                  max_cells_query = 5000,
                                  max_cells_ref = 5000) {

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

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             cell_types = cell_types,
                             pc_subset = pc_subset,
                             assay_name = assay_name,
                             max_cells_ref = max_cells_ref,
                             max_cells_query = max_cells_query)

    # Set data for Cramer test
    cell_list <- split(pca_output, pca_output[["cell_type"]])

    # Set the PC variables
    pc_vars <- paste0("PC", pc_subset)

    # Adding regression summaries for each cell type
    cramer_test <- vector("list", length = length(cell_list))
    names(cramer_test) <- cell_types
    for(cell_type in cell_types){

        dataset_ind <- cell_list[[cell_type]][, "dataset"] == "Reference"
        cramer_test[[cell_type]] <- cramer::cramer.test(
            as.matrix(cell_list[[cell_type]][dataset_ind, pc_vars]),
            as.matrix(cell_list[[cell_type]][!dataset_ind, pc_vars]),
            kernel = "phiBahr")
    }
    p_values <- unlist(lapply(cramer_test, function(t) t[["p.value"]]))

    # Return Cramer test output
    return(p_values)
}
