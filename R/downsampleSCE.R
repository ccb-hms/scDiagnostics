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
#' through standard SCE subsetting operations. This function does not preserve PCA rotation matrices or
#' other model-specific attributes.
#'
#' @param sce A \code{\linkS4class{SingleCellExperiment}} object to potentially downsample.
#' May contain PCA, UMAP, TSNE, or other reducedDims.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells,
#' it is returned unchanged. If NULL, no downsampling is performed (all cells are kept).
#' Default is 2500.
#' @param cell_types A character vector specifying which cell types to retain. If NULL,
#' no cell type filtering is performed. Default is NULL.
#' @param cell_type_col The column name in colData that contains cell type information.
#' Required if cell_types is not NULL. Default is NULL.
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
# Function to downsample SCE objects
downsampleSCE <- function(sce,
                          max_cells = 2500,
                          cell_types = NULL,
                          cell_type_col = NULL,
                          seed = NULL) {

    # Check if sce is a SingleCellExperiment object
    if (!is(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment object.")
    }

    # Check max_cells argument
    if (!is.null(max_cells)) {
        if (!is.numeric(max_cells) || length(max_cells) != 1 || max_cells <= 0 || max_cells != as.integer(max_cells)) {
            stop("'max_cells' must be a single positive integer or NULL.")
        }
        max_cells <- as.integer(max_cells)
    }

    # Check cell_types and cell_type_col arguments
    if (!is.null(cell_types)) {
        if (!is.character(cell_types)) {
            stop("'cell_types' must be a character vector or NULL.")
        }
        if (is.null(cell_type_col) || !is.character(cell_type_col) || length(cell_type_col) != 1) {
            stop("'cell_type_col' must be a single character string when 'cell_types' is provided.")
        }
        if (!cell_type_col %in% names(SummarizedExperiment::colData(sce))) {
            stop(paste("Column '", cell_type_col, "' not found in sce colData"))
        }
    }

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
        cell_type_mask <- sce[[cell_type_col]] %in% cell_types
        if (sum(cell_type_mask) == 0) {
            warning("No cells found for the specified cell types. Returning empty SCE object.")
            return(sce[, FALSE])  # Return empty SCE with same genes but no cells
        }
        sce <- sce[, cell_type_mask]
    }

    # Check if downsampling is needed
    n_cells <- ncol(sce)

    # If max_cells is NULL, return object without downsampling (but potentially cell-type filtered)
    if (is.null(max_cells)) {
        return(sce)
    }

    # If we already have fewer cells than max_cells, return as is
    if (n_cells <= max_cells) {
        return(sce)
    }

    # Step 3: Downsample cells
    cell_indices <- sample(n_cells, max_cells)
    sce_subset <- sce[, cell_indices]

    # Return the subsetted object (reducedDims coordinates are automatically preserved)
    return(sce_subset)
}
