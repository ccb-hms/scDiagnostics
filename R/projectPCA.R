#' @title Project Query Data Onto PCA Space of Reference Data
#'
#' @description
#' This function projects a query singleCellExperiment object onto the PCA space of a reference
#' singleCellExperiment object. The PCA analysis on the reference data is assumed to be pre-computed
#' and stored within the object. Optionally filters by cell types and downsamples the results.
#'
#' @details
#' This function assumes that the "PCA" element exists within the \code{reducedDims} of the reference data
#' (obtained using \code{reducedDim(reference_data)}) and that the genes used for PCA are present in both
#' the reference and query data. It performs centering and scaling of the query data based on the reference
#' data before projection using the FULL datasets to maintain proper mean centering. Cell type filtering
#' and downsampling are performed AFTER projection to preserve the statistical properties of the PCA space.
#' Cell names from the original SCE objects are preserved as rownames in the output.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix
#' for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix
#' for the reference cells.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data}
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data}
#' that identifies the cell types.
#' @param cell_types A character vector specifying which cell types to retain in the output. If NULL,
#' no cell type filtering is performed. Default is NULL.
#' @param pc_subset A numeric vector specifying the subset of principal components (PCs) to compare. Default is 1:10.
#' @param assay_name Name of the assay on which to perform computations. Defaults to \code{"logcounts"}.
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is NULL.
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is NULL.
#'
#' @return A \code{data.frame} containing the projected data in rows (reference and query data combined),
#' optionally filtered by cell types and downsampled. Rownames preserve the original cell names from
#' the SCE objects.
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
#' # Project the query data onto PCA space of reference
#' pca_output <- projectPCA(query_data = query_data,
#'                          reference_data = reference_data,
#'                          query_cell_type_col = "SingleR_annotation",
#'                          ref_cell_type_col = "expert_annotation",
#'                          pc_subset = 1:10)
#'
#' # Project with cell type filtering and balanced downsampling
#' pca_output_filtered <- projectPCA(query_data = query_data,
#'                                   reference_data = reference_data,
#'                                   query_cell_type_col = "SingleR_annotation",
#'                                   ref_cell_type_col = "expert_annotation",
#'                                   pc_subset = 1:5,
#'                                   cell_types = c("CD4", "CD8"),
#'                                   max_cells_ref = 1000,
#'                                   max_cells_query = 1000)
#'
# Function to project query data onto PCA space of reference data
projectPCA <- function(query_data,
                       reference_data,
                       query_cell_type_col,
                       ref_cell_type_col,
                       cell_types = NULL,
                       pc_subset = 1:10,
                       assay_name = "logcounts",
                       max_cells_ref = NULL,
                       max_cells_query = NULL){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Validate cell_types argument
    if (!is.null(cell_types)) {
        if (!is.character(cell_types)) {
            stop("\'cell_types\' must be a character vector or NULL.")
        }
    }

    # Validate max_cells_ref argument
    if (!is.null(max_cells_ref)) {
        if (!is.numeric(max_cells_ref) || length(max_cells_ref) != 1 ||
            max_cells_ref <= 0 || max_cells_ref != as.integer(max_cells_ref)) {
            stop("\'max_cells_ref\' must be a single positive integer or NULL.")
        }
        max_cells_ref <- as.integer(max_cells_ref)
    }

    # Validate max_cells_query argument
    if (!is.null(max_cells_query)) {
        if (!is.numeric(max_cells_query) || length(max_cells_query) != 1 ||
            max_cells_query <= 0 || max_cells_query != as.integer(max_cells_query)) {
            stop("\'max_cells_query\' must be a single positive integer or NULL.")
        }
        max_cells_query <- as.integer(max_cells_query)
    }

    # Extract reference PCA components and rotation matrix
    ref_mat <- reducedDim(reference_data, "PCA")[, pc_subset, drop = FALSE]
    rotation_mat <- attributes(reducedDim(
        reference_data, "PCA"))[["rotation"]][, pc_subset, drop = FALSE]
    PCA_genes <- rownames(rotation_mat)

    # Check if genes used for PCA are available in query data
    if (!all(PCA_genes %in% rownames(assay(query_data)))) {
        stop("Genes in reference PCA are not found in query data.")
    }

    # Center query data using full reference dataset (maintains original PCA centering)
    ref_assay <- assay(reference_data, assay_name)
    centering_vec <- Matrix::rowMeans(ref_assay)[PCA_genes]
    query_assay_subset <- assay(query_data, assay_name)[PCA_genes, , drop = FALSE]

    # Ensure gene order matches
    centering_vec <- centering_vec[rownames(rotation_mat)]
    query_assay_subset <- query_assay_subset[rownames(rotation_mat), , drop = FALSE]

    # Use sparse-compatible operations
    query_transposed <- Matrix::t(query_assay_subset)
    query_centered <- sweep(query_transposed, 2, centering_vec, "-")
    query_mat <- as.matrix(query_centered %*% rotation_mat)


    # Get cell names from the original SCE objects
    ref_cell_names <- colnames(reference_data)
    query_cell_names <- colnames(query_data)

    # Create full output dataframe with all cells
    full_output <- data.frame(
        rbind(ref_mat, query_mat),
        dataset = c(rep("Reference", nrow(ref_mat)),
                    rep("Query", nrow(query_mat))),
        cell_type = c(ifelse(rep(is.null(ref_cell_type_col), nrow(ref_mat)),
                             rep(NA, nrow(ref_mat)),
                             reference_data[[ref_cell_type_col]]),
                      ifelse(rep(is.null(query_cell_type_col), nrow(query_mat)),
                             rep(NA, nrow(query_mat)),
                             query_data[[query_cell_type_col]])),
        stringsAsFactors = FALSE
    )

    # Set rownames to preserve cell names
    rownames(full_output) <- c(ref_cell_names, query_cell_names)

    # Filter by cell types if specified
    if (!is.null(cell_types)) {
        # Check if specified cell types exist in the data
        available_types <- unique(full_output[["cell_type"]])
        available_types <- available_types[!is.na(available_types)]

        missing_types <- setdiff(cell_types, available_types)
        if (length(missing_types) > 0) {
            warning("Cell types not found in data: ", paste(missing_types, collapse = ", "))
        }

        # Filter to specified cell types
        cell_type_mask <- full_output[["cell_type"]] %in% cell_types
        if (sum(cell_type_mask) == 0) {
            warning("No cells found for the specified cell types. Returning empty data frame.")
            return(full_output[FALSE, ])
        }

        full_output <- full_output[which(cell_type_mask), ]
    }

    # Downsample each dataset separately if specified
    ref_indices <- which(full_output[["dataset"]] == "Reference")
    query_indices <- which(full_output[["dataset"]] == "Query")

    selected_indices <- c()

    # Downsample reference cells
    if (!is.null(max_cells_ref) && length(ref_indices) > max_cells_ref) {
        selected_ref_indices <- sample(ref_indices, max_cells_ref)
        selected_indices <- c(selected_indices, selected_ref_indices)
    } else {
        selected_indices <- c(selected_indices, ref_indices)
    }

    # Downsample query cells
    if (!is.null(max_cells_query) && length(query_indices) > max_cells_query) {
        selected_query_indices <- sample(query_indices, max_cells_query)
        selected_indices <- c(selected_indices, selected_query_indices)
    } else {
        selected_indices <- c(selected_indices, query_indices)
    }

    # Apply the selection (rownames are automatically preserved)
    full_output <- full_output[selected_indices, ]

    # Do NOT reset rownames - keep the original cell names
    return(full_output)
}
