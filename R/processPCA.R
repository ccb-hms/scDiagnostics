#' @title Process PCA for SingleCellExperiment Objects
#'
#' @description
#' This function ensures that a \code{\linkS4class{SingleCellExperiment}} object has valid PCA computed using
#' highly variable genes when needed. It only performs downsampling when PCA computation
#' is required, preserving existing valid PCA computations without modification.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Checks if PCA exists and is valid in the provided \code{\linkS4class{SingleCellExperiment}} object
#'   \item Validates PCA integrity including rotation matrix, percentVar, gene consistency, and dimensions
#'   \item If PCA is valid, returns the object unchanged (no downsampling)
#'   \item If PCA is missing or invalid and dataset is large, downsamples before computing PCA
#'   \item Computes PCA using highly variable genes when PCA is missing or invalid
#'   \item Utilizes scran for HVG selection and scater for PCA computation (soft dependencies)
#' }
#'
#' The downsampling strategy uses random sampling without replacement and only occurs
#' when PCA computation is necessary. This preserves expensive pre-computed PCA results
#' while ensuring computational efficiency for new PCA computations.
#'
#' PCA validation includes checking for:
#' \itemize{
#'   \item Presence of PCA in reducedDims
#'   \item Existence of rotation matrix and percentVar attributes
#'   \item Gene consistency between rotation matrix and current assay
#'   \item Dimension consistency between PCA coordinates and cell count
#' }
#'
#' @param sce_object A \code{\linkS4class{SingleCellExperiment}} object to process.
#' @param assay_name Name of the assay to use for HVG selection and PCA computation.
#' Should contain log-normalized expression values. Default is "logcounts".
#' @param n_hvgs Number of highly variable genes to select for PCA computation. Default is 2000.
#' @param max_cells Maximum number of cells to retain if downsampling is needed for PCA computation.
#' If NULL, no downsampling is performed. Default is NULL.
#'
#' @return A \code{\linkS4class{SingleCellExperiment}} object with valid PCA in the reducedDims slot,
#' including rotation matrix and percentVar attributes. Will have original cell count if PCA was valid,
#' or at most max_cells if PCA was computed.
#'
#' @note
#' This function requires the scran and scater packages for HVG selection and PCA computation.
#' These packages should be installed via BiocManager::install(c("scran", "scater")).
#'
#' Objects with existing valid PCA are returned unchanged to preserve expensive pre-computations.
#' Only datasets requiring PCA computation are subject to downsampling.
#'
#' @examples
#' # Load and prepare dataset
#' library(TENxPBMCData)
#' library(scuttle)
#'
#' pbmc_data <- TENxPBMCData("pbmc3k")
#' pbmc_subset <- pbmc_data[, 1:500]
#' pbmc_subset <- logNormCounts(pbmc_subset)
#'
#' # Remove any existing PCA
#' reducedDims(pbmc_subset) <- list()
#'
#' # Process dataset - will compute PCA using HVGs
#' processed_data <- processPCA(sce_object = pbmc_subset, n_hvgs = 1000)
#'
#' # Check results
#' "PCA" %in% reducedDimNames(processed_data)  # Should be TRUE
#' ncol(processed_data)  # Should be 500 (unchanged)
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to process SCE objects with PCA computation
processPCA <- function(sce_object,
                       assay_name = "logcounts",
                       n_hvgs = 2000,
                       max_cells = NULL) {

    # Validate input
    if (!is(sce_object, "SingleCellExperiment")) {
        stop("'sce_object' must be a SingleCellExperiment object.")
    }

    if (!is.null(assay_name) &&
        !(assay_name %in% SummarizedExperiment::assayNames(sce_object))) {
        stop("'sce_object' does not contain the specified assay.")
    }

    # Check max_cells
    if (!is.null(max_cells)) {
        if (!is.numeric(max_cells) || max_cells <= 0 ||
            max_cells != as.integer(max_cells)) {
            stop("'max_cells' must be a positive integer.")
        }
    }

    # Check n_hvgs
    if (!is.numeric(n_hvgs) || length(n_hvgs) != 1 ||
        n_hvgs <= 0 || n_hvgs != as.integer(n_hvgs)) {
        stop("'n_hvgs' must be a single positive integer.")
    }
    n_hvgs <- as.integer(n_hvgs)

    # Helper function to validate PCA integrity
    .validatePCA <- function(sce, assay_name) {
        # Check if PCA exists
        if (!"PCA" %in% reducedDimNames(sce)) {
            return(FALSE)
        }

        # Check if rotation matrix exists
        pca_attrs <- attributes(reducedDim(sce, "PCA"))
        rotation_mat <- pca_attrs[["rotation"]]
        if (is.null(rotation_mat)) {
            return(FALSE)
        }

        # Check if percentVar exists
        percent_var <- pca_attrs[["percentVar"]]
        if (is.null(percent_var)) {
            return(FALSE)
        }

        # Check if genes in rotation matrix exist in the assay
        pca_genes <- rownames(rotation_mat)
        assay_genes <- rownames(assay(sce, assay_name))
        if (!all(pca_genes %in% assay_genes)) {
            return(FALSE)
        }

        # Check dimension consistency
        pca_coords <- reducedDim(sce, "PCA")
        if (nrow(pca_coords) != ncol(sce)) {
            return(FALSE)
        }

        return(TRUE)
    }

    # Helper function to compute PCA with HVGs
    .computePCAWithHvgs <- function(sce, n_hvgs, assay_name) {

        # Check if required packages are available
        if (!requireNamespace("scran", quietly = TRUE)) {
            stop("Package 'scran' is required but not installed. Please install it with: BiocManager::install('scran')")
        }
        if (!requireNamespace("scater", quietly = TRUE)) {
            stop("Package 'scater' is required but not installed. Please install it with: BiocManager::install('scater')")
        }

        # Get HVGs
        var_stats <- scran::modelGeneVar(sce, assay.type = assay_name)
        hvg_genes <- scran::getTopHVGs(var_stats, n = n_hvgs)

        message("Using ", length(hvg_genes), " highly variable genes for PCA computation")

        # Compute PCA on HVGs
        sce <- scater::runPCA(sce,
                              assay.type = assay_name,
                              subset_row = hvg_genes)

        return(sce)
    }

    # Check if PCA is valid
    n_cells <- ncol(sce_object)
    has_valid_pca <- .validatePCA(sce_object, assay_name)

    if (has_valid_pca) {
        # PCA exists AND is valid - return unchanged (no downsampling)
        message("Data already has valid PCA - returning unchanged (", n_cells, " cells)")
        return(sce_object)
    } else {
        # PCA is missing or invalid - needs to be computed
        if ("PCA" %in% reducedDimNames(sce_object)) {
            message("Data has invalid PCA - recomputing...")
        } else {
            message("Data missing PCA - computing...")
        }

        # Downsample if needed
        processed_sce <- sce_object
        if (!is.null(max_cells) && n_cells > max_cells) {
            message("Downsampling data from ", n_cells, " to ",
                    max_cells, " cells before PCA computation")
            cell_indices <- sample(n_cells, max_cells)
            processed_sce <- sce_object[, cell_indices]
        }

        # Compute PCA
        message("Computing PCA...")
        processed_sce <- .computePCAWithHvgs(processed_sce,
                                             n_hvgs,
                                             assay_name)

        return(processed_sce)
    }
}
