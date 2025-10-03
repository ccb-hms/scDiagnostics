#' @title Downsample SingleCellExperiment Objects
#'
#' @description
#' This internal function downsamples \code{\linkS4class{SingleCellExperiment}} objects while preserving
#' reducedDims coordinate information (PCA, UMAP, t-SNE, etc.). Optionally, it can also subset by cell types
#' before downsampling.
#'
#' @details
#' The function can perform cell type filtering followed by random downsampling without replacement
#' when the number of cells exceeds the specified threshold. All reducedDims coordinates are preserved
#' through standard sce_object subsetting operations. This function does not preserve PCA rotation matrices or
#' other model-specific attributes.
#'
#' @param sce_object A \code{\linkS4class{SingleCellExperiment}} object to potentially downsample.
#' May contain PCA, UMAP, TSNE, or other reducedDims.
#' @param cell_type_col The column name in colData that contains cell type information.
#' Required if cell_types is not NULL. Default is NULL.
#' @param cell_types A character vector specifying which cell types to retain. If NULL,
#' no cell type filtering is performed. Default is NULL.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells,
#' it is returned unchanged. If NULL, no downsampling is performed (all cells are kept).
#' Default is 2500.
#' @param seed Random seed for reproducible downsampling. If NULL, no seed is set.
#' Default is NULL.
#'
#' @keywords internal
#'
#' @return A \code{\linkS4class{SingleCellExperiment}} object with at most max_cells cells
#' and optionally filtered by cell types. ReducedDims coordinates are preserved through
#' standard subsetting. If max_cells is NULL, the original object is returned unchanged.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to downsample sce_object objects
downsampleSCE <- function(sce_object,
                          cell_type_col = NULL,
                          cell_types = NULL,
                          max_cells = 2500,
                          seed = NULL) {

    # Check standard input arguments
    argumentCheck(query_data = sce_object,
                  query_cell_type_col = cell_type_col,
                  max_cells_query = max_cells)

    # Convert cell type columns to character if needed
    sce_object <- convertColumnsToCharacter(sce_object = sce_object,
                                            convert_cols = cell_type_col)

    # Select cell types
    cell_types <- selectCellTypes(query_data = sce_object,
                                  query_cell_type_col = cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = NULL)

    # Check seed argument
    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed)) {
            stop("'seed' must be a single integer or NULL.")
        }
        seed <- as.integer(seed)
    }

    # Set seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Subset by cell types if requested
    if (!is.null(cell_types)) {
        cell_type_mask <- sce_object[[cell_type_col]] %in% cell_types
        if (sum(cell_type_mask) == 0) {
            warning("No cells found for the specified cell types. Returning empty sce_object object.")
            return(sce_object[, FALSE])  # Return empty sce_object with same genes but no cells
        }
        sce_object <- sce_object[, cell_type_mask]
    }

    # Check if downsampling is needed
    n_cells <- ncol(sce_object)

    # If max_cells is NULL, return object without downsampling (but potentially cell-type filtered)
    if (is.null(max_cells)) {
        return(sce_object)
    }

    # If we already have fewer cells than max_cells, return as is
    if (n_cells <= max_cells) {
        return(sce_object)
    }

    # Step 3: Downsample cells
    cell_indices <- sample(n_cells, max_cells)
    sce_object_subset <- sce_object[, cell_indices]

    # Return the subsetted object (reducedDims coordinates are automatically preserved)
    return(sce_object_subset)
}
