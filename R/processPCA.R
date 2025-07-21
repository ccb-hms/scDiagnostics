#' @title Process PCA for SingleCellExperiment Objects
#'
#' @description
#' This function ensures that SingleCellExperiment objects have PCA computed using
#' highly variable genes when needed. It only performs downsampling when PCA computation
#' is required, preserving existing PCA computations without modification.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Checks if PCA exists in the provided SingleCellExperiment objects
#'   \item If PCA already exists, returns the object unchanged (no downsampling)
#'   \item If PCA is missing and dataset is large, downsamples before computing PCA
#'   \item Computes PCA using highly variable genes when PCA is missing
#'   \item When both datasets are provided and one has existing PCA, uses the genes from the existing PCA for the other dataset
#'   \item When both datasets lack PCA, uses the intersection of highly variable genes for consistent PCA spaces
#'   \item Utilizes scran for HVG selection and scater for PCA computation (soft dependencies)
#' }
#'
#' The downsampling strategy uses random sampling without replacement and only occurs
#' when PCA computation is necessary. This preserves expensive pre-computed PCA results
#' while ensuring computational efficiency for new PCA computations.
#'
#' For datasets where one has existing PCA and the other doesn't, the function uses the genes
#' from the existing PCA rotation matrix to ensure compatible PCA spaces. When both datasets
#' lack PCA, it uses the intersection of highly variable genes from both datasets.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object for query data.
#' If NULL, no processing is performed on query data. Default is NULL.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object for reference data.
#' If NULL, no processing is performed on reference data. Default is NULL.
#' @param max_cells Maximum number of cells to retain per dataset when PCA computation is required.
#' Objects with more cells will be randomly downsampled to this number before PCA computation.
#' Objects with existing PCA are never downsampled. Default is 5000.
#' @param assay_name Name of the assay to use for HVG selection and PCA computation.
#' Should contain log-normalized expression values. Default is "logcounts".
#' @param n_hvgs Number of highly variable genes to select for PCA computation.
#' When both datasets lack PCA, this number is used for each dataset before
#' taking the intersection. Default is 2000.
#'
#' @return When only one dataset is provided, returns the processed SingleCellExperiment object directly.
#' When both datasets are provided, returns a list containing:
#' \itemize{
#'   \item query_data: Processed query data
#'   \item reference_data: Processed reference data
#' }
#' Each returned object will have:
#' \itemize{
#'   \item PCA in reducedDims slot
#'   \item Original cell count if PCA existed, or at most max_cells if PCA was computed
#' }
#'
#' @note
#' This function requires the scran and scater packages for HVG selection and PCA computation.
#' These packages should be installed via BiocManager::install(c("scran", "scater")).
#'
#' Objects with existing PCA are returned unchanged to preserve expensive pre-computations.
#' Only datasets requiring PCA computation are subject to downsampling.
#'
#' @examples
#' # Example 1: Process single dataset
#' library(TENxPBMCData)
#' library(scuttle)
#'
#' # Load and prepare dataset
#' pbmc_data <- TENxPBMCData("pbmc3k")
#' pbmc_subset <- pbmc_data[, 1:500]  # Use subset for example
#' pbmc_subset <- logNormCounts(pbmc_subset)
#'
#' # Remove any existing PCA
#' reducedDims(pbmc_subset) <- list()
#'
#' # Process dataset - will compute PCA using HVGs
#' processed_data <- processPCA(query_data = pbmc_subset, n_hvgs = 1000)
#'
#' # Check results
#' "PCA" %in% reducedDimNames(processed_data)  # Should be TRUE
#' ncol(processed_data)  # Should be 500 (unchanged)
#'
#' # Example 2: Process both query and reference datasets
#' # Create reference dataset
#' pbmc_ref <- pbmc_data[, 501:1000]
#' pbmc_ref <- logNormCounts(pbmc_ref)
#' reducedDims(pbmc_ref) <- list()
#'
#' # Process both datasets - will use intersection of HVGs
#' result <- processPCA(
#'   query_data = pbmc_subset,
#'   reference_data = pbmc_ref,
#'   max_cells = 10000,
#'   n_hvgs = 800
#' )
#'
#' # Check results
#' "PCA" %in% reducedDimNames(result$query_data)      # Should be TRUE
#' "PCA" %in% reducedDimNames(result$reference_data)  # Should be TRUE
#' ncol(result$query_data)      # Should be 500
#' ncol(result$reference_data)  # Should be 500
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to process SCE objects with PCA computation
processPCA <- function(query_data = NULL,
                       reference_data = NULL,
                       max_cells = 5000,
                       n_hvgs = 2000,
                       assay_name = "logcounts") {

    # Validate input using argumentCheck for what it can handle
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  assay_name = assay_name)

    # Additional argument checks not covered by argumentCheck
    if (is.null(query_data) && is.null(reference_data)) {
        stop("At least one of query_data or reference_data must be provided.")
    }

    # Check max_cells
    if (!is.numeric(max_cells) || length(max_cells) != 1 || max_cells <= 0 || max_cells != as.integer(max_cells)) {
        stop("'max_cells' must be a single positive integer.")
    }
    max_cells <- as.integer(max_cells)

    # Check n_hvgs
    if (!is.numeric(n_hvgs) || length(n_hvgs) != 1 || n_hvgs <= 0 || n_hvgs != as.integer(n_hvgs)) {
        stop("'n_hvgs' must be a single positive integer.")
    }
    n_hvgs <- as.integer(n_hvgs)

    # Helper function to compute PCA with HVGs
    .computePCAWithHvgs <- function(sce, hvg_genes = NULL) {

        # Check if required packages are available
        if (!requireNamespace("scran", quietly = TRUE)) {
            stop("Package 'scran' is required but not installed. Please install it with: BiocManager::install('scran')")
        }
        if (!requireNamespace("scater", quietly = TRUE)) {
            stop("Package 'scater' is required but not installed. Please install it with: BiocManager::install('scater')")
        }

        # Select HVGs if not provided
        if (!is.null(hvg_genes)) {
            # Check which genes actually exist in the dataset
            available_genes <- intersect(hvg_genes, rownames(sce))

            if (length(available_genes) < length(hvg_genes)) {
                missing_count <- length(hvg_genes) - length(available_genes)
                message("Note: ", missing_count, " genes from existing PCA not found in dataset. Using ", length(available_genes), " available genes.")
            }

            if (length(available_genes) < 50) {
                stop("Too few genes (", length(available_genes), ") available for PCA computation. Consider using datasets with more gene overlap.")
            }

            hvg_genes <- available_genes
        }

        # Compute PCA on HVGs
        sce <- scater::runPCA(sce,
                              assay.type = assay_name,
                              subset_row = hvg_genes)

        return(sce)
    }

    # Helper function to extract genes from existing PCA
    .getGenesFromPCA <- function(sce) {
        if ("PCA" %in% reducedDimNames(sce)) {
            rotation_matrix <- attr(reducedDim(sce, "PCA"), "rotation")
            if (!is.null(rotation_matrix)) {
                return(rownames(rotation_matrix))
            }
        }
        return(NULL)
    }

    # Process each dataset
    processed_query <- NULL
    processed_reference <- NULL
    genes_to_use <- NULL

    # Check PCA status for both datasets
    query_has_pca <- !is.null(query_data) && "PCA" %in% reducedDimNames(query_data)
    ref_has_pca <- !is.null(reference_data) && "PCA" %in% reducedDimNames(reference_data)

    # Determine which genes to use for PCA computation
    if (!is.null(query_data) && !is.null(reference_data)) {

        if (query_has_pca && !ref_has_pca) {
            # Use genes from query PCA for reference
            genes_to_use <- .getGenesFromPCA(query_data)
            message("Using ", length(genes_to_use),
                    " genes from existing query PCA for reference data")

        } else if (!query_has_pca && ref_has_pca) {
            # Use genes from reference PCA for query
            genes_to_use <- .getGenesFromPCA(reference_data)
            message("Using ", length(genes_to_use),
                    " genes from existing reference PCA for query data")

        } else if (!query_has_pca && !ref_has_pca) {
            # Both need PCA - use intersection of HVGs
            if (!requireNamespace("scran", quietly = TRUE)) {
                stop("Package 'scran' is required but not installed. Please install it with: BiocManager::install('scran')")
            }

            # Get HVGs for each dataset
            query_var <- scran::modelGeneVar(query_data, assay.type = assay_name)
            ref_var <- scran::modelGeneVar(reference_data, assay.type = assay_name)

            query_hvgs <- scran::getTopHVGs(query_var, n = n_hvgs)
            ref_hvgs <- scran::getTopHVGs(ref_var, n = n_hvgs)

            # Take intersection
            genes_to_use <- intersect(query_hvgs, ref_hvgs)

            if (length(genes_to_use) < 100) {
                warning("Only ", length(genes_to_use),
                        " common HVGs found between datasets. Consider increasing n_hvgs parameter.")
            }

            message("Using ", length(genes_to_use), " common HVGs for both datasets")
        }
        # If both have PCA, genes_to_use remains NULL (no new computation needed)
    }

    # Process query data
    if (!is.null(query_data)) {

        n_cells_query <- ncol(query_data)
        has_pca_query <- "PCA" %in% reducedDimNames(query_data)

        if (has_pca_query) {
            # PCA exists - return unchanged (no downsampling)
            message("Query data already has PCA - returning unchanged (", n_cells_query, " cells)")
            processed_query <- query_data
        } else {
            # PCA needs to be computed
            if (n_cells_query > max_cells) {
                # Downsample before computing PCA
                message("Downsampling query data from ", n_cells_query, " to ", max_cells, " cells before PCA computation")
                cell_indices <- sample(n_cells_query, max_cells)
                processed_query <- query_data[, cell_indices]
            } else {
                processed_query <- query_data
            }

            # Compute PCA
            message("Computing PCA for query data (", ncol(processed_query), " cells)")
            processed_query <- .computePCAWithHvgs(processed_query, genes_to_use)
        }
    }

    # Process reference data
    if (!is.null(reference_data)) {

        n_cells_ref <- ncol(reference_data)
        has_pca_ref <- "PCA" %in% reducedDimNames(reference_data)

        if (has_pca_ref) {
            # PCA exists - return unchanged (no downsampling)
            message("Reference data already has PCA - returning unchanged (", n_cells_ref, " cells)")
            processed_reference <- reference_data
        } else {
            # PCA needs to be computed
            if (n_cells_ref > max_cells) {
                # Downsample before computing PCA
                message("Downsampling reference data from ", n_cells_ref, " to ", max_cells, " cells before PCA computation")
                cell_indices <- sample(n_cells_ref, max_cells)
                processed_reference <- reference_data[, cell_indices]
            } else {
                processed_reference <- reference_data
            }

            # Compute PCA
            message("Computing PCA for reference data (", ncol(processed_reference), " cells)")
            processed_reference <- .computePCAWithHvgs(processed_reference, genes_to_use)
        }
    }

    # Return logic based on what was provided
    if (!is.null(query_data) && !is.null(reference_data)) {
        # Both datasets provided - return list
        return(list(
            query_data = processed_query,
            reference_data = processed_reference
        ))
    } else if (!is.null(query_data)) {
        # Only query provided - return SCE directly
        return(processed_query)
    } else if (!is.null(reference_data)) {
        # Only reference provided - return SCE directly
        return(processed_reference)
    }
}
