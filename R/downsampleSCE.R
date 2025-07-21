#' @title Downsample SingleCellExperiment Objects
#'
#' @description
#' This internal function downsamples \code{\linkS4class{SingleCellExperiment}} objects while preserving
#' all reducedDims information, including PCA, UMAP, t-SNE, and any other dimensionality
#' reduction results with their associated attributes.
#'
#' @details
#' The function performs random downsampling without replacement when the number of
#' cells exceeds the specified threshold. All reducedDims are preserved by subsetting
#' the cell coordinates and manually restoring associated attributes that would
#' otherwise be lost during subsetting operations.
#'
#' Supported reducedDims and their preserved attributes:
#' \itemize{
#'   \item PCA: rotation, percentVar, varExplained
#'   \item UMAP: any custom attributes from uwot or other UMAP implementations
#'   \item TSNE: any custom attributes from Rtsne or other t-SNE implementations
#'   \item Any other reducedDim: all associated attributes
#' }
#'
#' @param sce A \code{\linkS4class{SingleCellExperiment}} object to potentially downsample.
#' Must contain PCA in reducedDims. May also contain UMAP, TSNE, or other reducedDims.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells,
#' it is returned unchanged. Default is 2500.
#' @param seed Random seed for reproducible downsampling. If NULL, no seed is set.
#' Default is NULL.
#'
#' @keywords internal
#'
#' @return A \code{\linkS4class{SingleCellExperiment}} object with at most max_cells cells.
#' If downsampling occurred, all reducedDims (PCA, UMAP, TSNE, etc.) and their
#' attributes are preserved.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to downsample SCE objects
downsampleSCE <- function(sce, max_cells = 2500, seed = NULL) {

    # Check if sce is a SingleCellExperiment object
    if (!is(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment object.")
    }

    # Check max_cells argument
    if (!is.numeric(max_cells) || length(max_cells) != 1 || max_cells <= 0 || max_cells != as.integer(max_cells)) {
        stop("'max_cells' must be a single positive integer.")
    }
    max_cells <- as.integer(max_cells)

    # Check seed argument
    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed)) {
            stop("'seed' must be a single integer or NULL.")
        }
        seed <- as.integer(seed)
    }

    # Check if downsampling is needed
    n_cells <- ncol(sce)
    if (n_cells <= max_cells) {
        return(sce)
    }

    # Set seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Sample cells
    cell_indices <- sample(n_cells, max_cells)

    # Store all reducedDims attributes before subsetting
    reduced_dims_attrs <- list()
    for (dim_name in reducedDimNames(sce)) {
        reduced_dims_attrs[[dim_name]] <-
            attributes(reducedDim(sce, dim_name))
    }

    # Perform subsetting
    sce_subset <- sce[, cell_indices]

    # Restore reducedDims attributes for all dimensionality reductions
    for (dim_name in names(reduced_dims_attrs)) {
        if (dim_name %in% reducedDimNames(sce_subset)) {

            # Get the subsetted coordinates
            dim_coords <- reducedDim(sce_subset, dim_name)
            original_attrs <- reduced_dims_attrs[[dim_name]]

            # Restore all attributes except automatic ones (dim, dimnames)
            for (attr_name in names(original_attrs)) {
                if (!attr_name %in% c("dim", "dimnames")) {
                    attr(dim_coords, attr_name) <-
                        original_attrs[[attr_name]]
                }
            }

            # Put back into SCE
            reducedDim(sce_subset, dim_name) <- dim_coords
        }
    }

    return(sce_subset)
}
