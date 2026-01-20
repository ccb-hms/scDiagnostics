#' @title Calculate Top Loading Gene Expression Shifts
#'
#' @description
#' This function identifies genes with the highest loadings for specified principal components
#' and performs statistical tests to detect distributional differences between query and reference data.
#' It also calculates the proportion of variance explained by each principal component within
#' specific cell types. Optionally, it can detect anomalous cells using isolation forests.
#'
#' @details
#' This function extracts the top loading genes for each specified principal component from the reference
#' PCA space and performs distributional comparisons between query and reference data. For each gene,
#' it performs statistical tests to identify genes that may be causing PC-specific alignment issues
#' between datasets. A key feature is the calculation of cell-type-specific variance explained by
#' global PCs, providing a more nuanced view of how major biological axes affect individual populations.
#' When anomaly detection is enabled, isolation forests are used to identify anomalous cells based on
#' their PCA projections.
#'
#' When \code{anomaly_comparison = TRUE}, the statistical analysis focuses specifically on
#' comparing non-anomalous reference cells against anomalous query cells. This can help
#' identify genes that are differentially expressed between "normal" reference cells and
#' potentially problematic query cells, providing insights into what makes certain query
#' cells anomalous.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to analyze. If NULL, all common cell types are used.
#' @param pc_subset A numeric vector specifying which principal components to analyze. Default is 1:5.
#' @param n_top_loadings Number of top loading genes to analyze per PC. Default is 50.
#' @param genes_to_analyze A character vector specifying genes to analyze. If NULL (default),
#'                         genes are selected based on top loadings from specified principal components (see \code{n_top_loadings}). Default is NULL.
#' @param p_value_threshold P-value threshold for statistical significance. Default is 0.05.
#' @param adjust_method Method for multiple testing correction. Default is "fdr".
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param detect_anomalies Logical indicating whether to perform anomaly detection using isolation forests.
#'                         Default is FALSE.
#' @param anomaly_comparison Logical indicating whether to perform statistical comparisons
#'                           between non-anomalous reference cells and anomalous query cells instead of all-vs-all
#'                           comparisons. When TRUE, only non-anomalous reference cells are compared against only
#'                           anomalous query cells for each cell type. Requires detect_anomalies = TRUE. Default is FALSE.
#' @param anomaly_threshold A numeric value specifying the threshold for identifying anomalies when
#'                          \code{detect_anomalies} is TRUE. Default is 0.6.
#' @param n_tree An integer specifying the number of trees for the isolation forest when
#'               \code{detect_anomalies} is TRUE. Default is 500.
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#'
#' @return A list containing:
#' \itemize{
#'   \item PC results: Named elements for each PC (e.g., "PC1", "PC2") containing data frames with gene-level analysis results.
#'   \item expression_data: Matrix of expression values for all analyzed genes (genes × cells).
#'   \item cell_metadata: Data frame with columns: cell_id, dataset, cell_type, original_index, and optionally anomaly_status.
#'   \item gene_metadata: Data frame with columns: gene, pc, loading for all analyzed genes.
#'   \item percent_var: Named numeric vector of global percent variance explained for each analyzed PC.
#'   \item cell_type_variance: A data frame detailing the percent of variance a global PC explains within specific cell types for both query and reference datasets.
#'   \item anomaly_results: If \code{detect_anomalies} is TRUE, contains the full output from \code{detectAnomaly}.
#' }
#'
#' The `cell_type_variance` data frame contains columns: pc, cell_type, dataset, percent_variance.
#' When anomaly detection is enabled, `cell_metadata` includes an additional `anomaly_status` column.
#'
#' @export
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateGeneShiftsObject}}, \code{\link{detectAnomaly}}
#'
#' @importFrom stats wilcox.test var p.adjust na.omit setNames
#'
# Function to calculate expression shifts for genes with top loadings
calculateGeneShifts <- function(query_data,
                                reference_data,
                                query_cell_type_col,
                                ref_cell_type_col,
                                cell_types = NULL,
                                pc_subset = 1:5,
                                n_top_loadings = 50,
                                genes_to_analyze = NULL,
                                p_value_threshold = 0.05,
                                adjust_method = "fdr",
                                assay_name = "logcounts",
                                detect_anomalies = FALSE,
                                anomaly_comparison = FALSE,
                                anomaly_threshold = 0.6,
                                n_tree = 500,
                                max_cells_query = 5000,
                                max_cells_ref = 5000) {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  ref_cell_type_col = ref_cell_type_col,
                  query_cell_type_col = query_cell_type_col,
                  pc_subset_query = pc_subset,
                  assay_name = assay_name,
                  max_cells_query = max_cells_query,
                  max_cells_ref = max_cells_ref)

    # Convert cell type columns to character if needed
    query_data <- convertColumnsToCharacter(sce_object = query_data,
                                            convert_cols = query_cell_type_col)
    reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                convert_cols = ref_cell_type_col)

    # Input validation
    if (!is.numeric(n_top_loadings) || length(n_top_loadings) != 1 || n_top_loadings <= 0) {
        stop("n_top_loadings must be a positive integer")
    }

    # Validate genes_to_analyze parameter
    if (!is.null(genes_to_analyze)) {
        if (!is.character(genes_to_analyze) || length(genes_to_analyze) == 0) {
            stop("genes_to_analyze must be a non-empty character vector or NULL")
        }
    }

    if (!is.numeric(p_value_threshold) || length(p_value_threshold) != 1 ||
        p_value_threshold < 0 || p_value_threshold > 1) {
        stop("p_value_threshold must be between 0 and 1")
    }

    # Validate anomaly detection parameters
    if (!is.logical(detect_anomalies) || length(detect_anomalies) != 1) {
        stop("detect_anomalies must be a logical value")
    }

    if (detect_anomalies) {
        if (!is.numeric(anomaly_threshold) || length(anomaly_threshold) != 1 ||
            anomaly_threshold <= 0 || anomaly_threshold >= 1) {
            stop("anomaly_threshold must be a numeric value between 0 and 1")
        }

        if (!is.numeric(n_tree) || length(n_tree) != 1 || n_tree <= 0 || n_tree != as.integer(n_tree)) {
            stop("n_tree must be a positive integer")
        }
    }

    # Validate anomaly comparison parameter
    if (!is.logical(anomaly_comparison) || length(anomaly_comparison) != 1) {
        stop("anomaly_comparison must be a logical value")
    }

    if (anomaly_comparison && !detect_anomalies) {
        stop("anomaly_comparison = TRUE requires detect_anomalies = TRUE")
    }

    # Helper function to filter cells by anomaly status
    .filterCellsByAnomalyStatus <- function(cell_indices, cell_names, anomaly_results,
                                            cell_type, dataset_type, keep_anomalous = TRUE) {
        if (is.null(anomaly_results) || !cell_type %in% names(anomaly_results)) {
            return(cell_indices)  # Return all cells if no anomaly data
        }

        # Get anomaly status for this cell type
        if (dataset_type == "reference") {
            anomaly_names <- names(anomaly_results[[cell_type]][["reference_anomaly"]])
            anomaly_status <- anomaly_results[[cell_type]][["reference_anomaly"]]
            # Strip "Reference_" prefix to match cell names
            anomaly_names_clean <- gsub("^Reference_", "", anomaly_names)
        } else {
            anomaly_names <- names(anomaly_results[[cell_type]][["query_anomaly"]])
            anomaly_status <- anomaly_results[[cell_type]][["query_anomaly"]]
            # Strip "Query_" prefix to match cell names
            anomaly_names_clean <- gsub("^Query_", "", anomaly_names)
        }

        if (is.null(anomaly_names) || length(anomaly_names) == 0) {
            return(cell_indices)  # Return all cells if no anomaly data
        }

        # Map cell indices to cell names and then to anomaly status
        current_cell_names <- cell_names[cell_indices]

        # Find which cells are anomalous
        cells_are_anomalous <- current_cell_names %in% anomaly_names_clean[anomaly_status]

        # Filter based on keep_anomalous parameter
        if (keep_anomalous) {
            filtered_indices <- cell_indices[cells_are_anomalous]
        } else {
            filtered_indices <- cell_indices[!cells_are_anomalous]
        }

        return(filtered_indices)
    }

    # Select cell types
    cell_types <- selectCellTypes(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = query_cell_type_col,
                                  ref_cell_type_col = ref_cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = TRUE,
                                  n_cell_types = NULL)

    # Ensure cell names exist for anomaly detection mapping
    # Store original cell names or create them if they don't exist
    original_ref_names <- colnames(reference_data)
    original_query_names <- colnames(query_data)

    if (is.null(original_ref_names) || any(is.na(original_ref_names)) ||
        length(unique(original_ref_names)) != length(original_ref_names)) {
        original_ref_names <- paste0("REF_CELL_", seq_len(ncol(reference_data)))
        colnames(reference_data) <- original_ref_names
    }

    if (is.null(original_query_names) || any(is.na(original_query_names)) ||
        length(unique(original_query_names)) != length(original_query_names)) {
        original_query_names <- paste0("QUERY_CELL_", seq_len(ncol(query_data)))
        colnames(query_data) <- original_query_names
    }

    # Get PCA results
    pca_attrs <- attributes(SingleCellExperiment::reducedDim(reference_data, "PCA"))
    pca_rotation <- pca_attrs[["rotation"]]
    percent_var <- pca_attrs[["percentVar"]]
    if (is.null(pca_rotation) || is.null(percent_var)) {
        stop("PCA rotation matrix and/or percentVar not found in reference_data's reducedDim 'PCA' attributes.")
    }
    if (max(pc_subset) > ncol(pca_rotation)) {
        stop("pc_subset contains indices greater than the number of PCs available in the reference.")
    }
    pc_percent_var <- percent_var[pc_subset]
    names(pc_percent_var) <- paste0("PC", pc_subset)

    # Perform anomaly detection FIRST (before any downsampling)
    # This ensures the rotation matrix is intact
    # Anomaly detection always uses all genes in the datasets (genome-wide)
    anomaly_results <- NULL
    if (detect_anomalies) {
        tryCatch({
            anomaly_results <- detectAnomaly(
                reference_data = reference_data,
                query_data = query_data,
                ref_cell_type_col = ref_cell_type_col,
                query_cell_type_col = query_cell_type_col,
                cell_types = cell_types,
                pc_subset = pc_subset,
                n_tree = n_tree,
                anomaly_threshold = anomaly_threshold,
                assay_name = assay_name,
                max_cells_ref = NULL,
                max_cells_query = NULL
            )
        }, error = function(e) {
            warning("Anomaly detection failed: ", e[["message"]], ". Continuing without anomaly detection.")
            anomaly_results <<- NULL
        })
    }

    # Now downsample the data (with cell type filtering)
    query_data <- downsampleSCE(sce_object = query_data,
                                max_cells = max_cells_query,
                                cell_types = cell_types,
                                cell_type_col = query_cell_type_col)
    reference_data <- downsampleSCE(sce_object = reference_data,
                                    max_cells = max_cells_ref,
                                    cell_types = cell_types,
                                    cell_type_col = ref_cell_type_col)

    # Check if cell type columns exist after downsampling
    if (!query_cell_type_col %in% names(SummarizedExperiment::colData(query_data))) {
        stop(paste("Column '", query_cell_type_col, "' not found in query_data colData"))
    }

    if (!ref_cell_type_col %in% names(SummarizedExperiment::colData(reference_data))) {
        stop(paste("Column '", ref_cell_type_col, "' not found in reference_data colData"))
    }

    # Get common genes
    common_genes <- intersect(rownames(query_data), rownames(reference_data))
    common_genes <- intersect(common_genes, rownames(pca_rotation))
    if (length(common_genes) == 0) {
        stop("No common genes found between query, reference, and PCA rotation matrix.")
    }
    pca_rotation <- pca_rotation[common_genes, ]

    # Prepare cell indices (data is already filtered by cell type and downsampled)
    query_cell_indices <- split(seq_len(ncol(query_data)), query_data[[query_cell_type_col]])
    ref_cell_indices <- split(seq_len(ncol(reference_data)), reference_data[[ref_cell_type_col]])
    query_cell_indices <- query_cell_indices[names(query_cell_indices) %in% cell_types]
    ref_cell_indices <- ref_cell_indices[names(ref_cell_indices) %in% cell_types]
    available_cell_types <- intersect(names(query_cell_indices), names(ref_cell_indices))
    if (length(available_cell_types) == 0) {
        warning("No common cell types with sufficient cells found between datasets.")
        return(list())
    }

    # Determine which genes to analyze
    if (!is.null(genes_to_analyze)) {
        # Use user-specified genes
        genes_to_use <- genes_to_analyze
        # Check if all specified genes are available in the data
        available_genes <- genes_to_use[genes_to_use %in% rownames(query_data) &
                                            genes_to_use %in% rownames(reference_data)]
        if (length(available_genes) == 0) {
            stop("None of the specified genes in genes_to_analyze are found in both query and reference data.")
        }
        if (length(available_genes) < length(genes_to_use)) {
            warning("Some genes in genes_to_analyze were not found in both datasets. Using ",
                    length(available_genes), " out of ", length(genes_to_use), " genes.")
        }
        genes_to_use <- available_genes

        # Create gene metadata for user-specified genes
        gene_metadata_list <- list()
        pc_name_placeholder <- paste0("PC", pc_subset[1])  # Use first PC as placeholder
        gene_metadata_list[[pc_name_placeholder]] <- data.frame(
            gene = genes_to_use,
            pc = NA_integer_,
            loading = NA_real_,
            stringsAsFactors = FALSE
        )

        all_top_genes <- genes_to_use
        pc_results <- stats::setNames(vector("list", length(pc_subset)), paste0("PC", pc_subset))

    } else {
        # Use original top loading gene selection logic
        all_top_genes <- character(0)
        gene_metadata_list <- list()
        pc_results <- stats::setNames(vector("list", length(pc_subset)), paste0("PC", pc_subset))

        for (pc in pc_subset) {
            pc_name <- paste0("PC", pc)
            pc_loadings <- pca_rotation[, pc]
            top_loading_indices <- order(abs(pc_loadings),
                                         decreasing = TRUE)[1:min(n_top_loadings, length(pc_loadings))]
            top_genes <- names(pc_loadings)[top_loading_indices]
            top_loadings_vals <- pc_loadings[top_loading_indices]

            all_top_genes <- unique(c(all_top_genes, top_genes))
            gene_metadata_list[[pc_name]] <- data.frame(gene = top_genes,
                                                        pc = pc,
                                                        loading = top_loadings_vals,
                                                        stringsAsFactors = FALSE)
        }
    }

    # Gene-level statistical analysis
    if (!is.null(genes_to_analyze)) {
        # Analysis for user-specified genes
        genes_to_use <- all_top_genes

        pc_result_list <- list()
        for (ct in available_cell_types) {
            query_cells_ct <- query_cell_indices[[ct]]
            ref_cells_ct <- ref_cell_indices[[ct]]

            # Apply anomaly filtering if requested
            if (anomaly_comparison && !is.null(anomaly_results)) {
                # Filter reference cells: keep only NON-anomalous
                ref_cells_ct <- .filterCellsByAnomalyStatus(
                    ref_cells_ct,
                    colnames(reference_data),
                    anomaly_results,
                    ct,
                    "reference",
                    keep_anomalous = FALSE
                )

                # Filter query cells: keep only ANOMALOUS
                query_cells_ct <- .filterCellsByAnomalyStatus(
                    query_cells_ct,
                    colnames(query_data),
                    anomaly_results,
                    ct,
                    "query",
                    keep_anomalous = TRUE
                )

                # Check if we have enough cells after filtering
                if (length(query_cells_ct) < 3) {
                    warning("Cell type '", ct, "' has fewer than 3 anomalous query cells. Skipping statistical analysis.")
                    next
                }
                if (length(ref_cells_ct) < 3) {
                    warning("Cell type '", ct, "' has fewer than 3 non-anomalous reference cells. Skipping statistical analysis.")
                    next
                }
            } else {
                # Original check for minimum cell counts
                if (length(query_cells_ct) < 3 || length(ref_cells_ct) < 3) next
            }

            query_expr_ct <- SummarizedExperiment::assay(query_data,
                                                         assay_name)[genes_to_use,
                                                                     query_cells_ct, drop = FALSE]
            ref_expr_ct <- SummarizedExperiment::assay(reference_data,
                                                       assay_name)[genes_to_use,
                                                                   ref_cells_ct, drop = FALSE]

            # For user-specified genes, use NA for loadings
            loadings_placeholder <- rep(NA_real_, length(genes_to_use))
            gene_results <- processGenesSimple(genes_to_use,
                                               loadings_placeholder,
                                               query_expr_ct,
                                               ref_expr_ct,
                                               ct)
            if (length(gene_results) > 0)
                pc_result_list <- c(pc_result_list, gene_results)
        }

        if (length(pc_result_list) > 0) {
            df <- do.call(rbind, pc_result_list)
            df[["p_adjusted"]] <- stats::p.adjust(df[["p_value"]], method = adjust_method)
            df[["significant"]] <- df[["p_adjusted"]] <= p_value_threshold
            df <- df[order(df[["p_adjusted"]]), ]
            rownames(df) <- NULL
            # Place results in first PC slot since genes_to_analyze aren't tied to specific PCs
            pc_results[[1]] <- df
        } else {
            pc_results[[1]] <- data.frame()
        }

    } else {
        # Original PC-specific analysis logic
        for (pc in pc_subset) {
            pc_name <- paste0("PC", pc)
            pc_loadings <- pca_rotation[, pc]
            top_loading_indices <- order(abs(pc_loadings),
                                         decreasing = TRUE)[1:min(n_top_loadings, length(pc_loadings))]
            top_genes <- names(pc_loadings)[top_loading_indices]
            top_loadings_vals <- pc_loadings[top_loading_indices]

            pc_result_list <- list()
            for (ct in available_cell_types) {
                query_cells_ct <- query_cell_indices[[ct]]
                ref_cells_ct <- ref_cell_indices[[ct]]

                # Apply anomaly filtering if requested
                if (anomaly_comparison && !is.null(anomaly_results)) {
                    # Filter reference cells: keep only NON-anomalous
                    ref_cells_ct <- .filterCellsByAnomalyStatus(
                        ref_cells_ct,
                        colnames(reference_data),
                        anomaly_results,
                        ct,
                        "reference",
                        keep_anomalous = FALSE
                    )

                    # Filter query cells: keep only ANOMALOUS
                    query_cells_ct <- .filterCellsByAnomalyStatus(
                        query_cells_ct,
                        colnames(query_data),
                        anomaly_results,
                        ct,
                        "query",
                        keep_anomalous = TRUE
                    )

                    # Check if we have enough cells after filtering
                    if (length(query_cells_ct) < 3) {
                        warning("Cell type '", ct, "' has fewer than 3 anomalous query cells. Skipping statistical analysis.")
                        next
                    }
                    if (length(ref_cells_ct) < 3) {
                        warning("Cell type '", ct, "' has fewer than 3 non-anomalous reference cells. Skipping statistical analysis.")
                        next
                    }
                } else {
                    # Original check for minimum cell counts
                    if (length(query_cells_ct) < 3 || length(ref_cells_ct) < 3) next
                }

                query_expr_ct <- SummarizedExperiment::assay(query_data,
                                                             assay_name)[top_genes,
                                                                         query_cells_ct, drop = FALSE]
                ref_expr_ct <- SummarizedExperiment::assay(reference_data,
                                                           assay_name)[top_genes,
                                                                       ref_cells_ct, drop = FALSE]

                gene_results <- processGenesSimple(top_genes,
                                                   top_loadings_vals,
                                                   query_expr_ct,
                                                   ref_expr_ct,
                                                   ct)
                if (length(gene_results) > 0)
                    pc_result_list <- c(pc_result_list, gene_results)
            }

            if (length(pc_result_list) > 0) {
                df <- do.call(rbind, pc_result_list)
                df[["p_adjusted"]] <- stats::p.adjust(df[["p_value"]], method = adjust_method)
                df[["significant"]] <- df[["p_adjusted"]] <= p_value_threshold
                pc_results[[pc_name]] <- df[order(df[["p_adjusted"]]), ]
                rownames(pc_results[[pc_name]]) <- NULL
            } else {
                pc_results[[pc_name]] <- data.frame()
            }
        }
    }

    gene_metadata <- do.call(rbind, gene_metadata_list)
    rownames(gene_metadata) <- NULL

    # Cell-type-specific variance explained by global PCs
    var_explained_list <- list()
    full_ref_assay <- SummarizedExperiment::assay(reference_data,
                                                  assay_name)[common_genes, ]
    full_query_assay <- SummarizedExperiment::assay(query_data,
                                                    assay_name)[common_genes, ]

    for (ct in available_cell_types) {
        ref_cells_ct <- ref_cell_indices[[ct]]
        query_cells_ct <- query_cell_indices[[ct]]

        # Apply the same anomaly filtering for variance calculation if anomaly_comparison is TRUE
        if (anomaly_comparison && !is.null(anomaly_results)) {
            ref_cells_ct <- .filterCellsByAnomalyStatus(
                ref_cells_ct,
                colnames(reference_data),
                anomaly_results,
                ct,
                "reference",
                keep_anomalous = FALSE
            )

            query_cells_ct <- .filterCellsByAnomalyStatus(
                query_cells_ct,
                colnames(query_data),
                anomaly_results,
                ct,
                "query",
                keep_anomalous = TRUE
            )
        }

        ref_expr_ct <- full_ref_assay[, ref_cells_ct, drop = FALSE]
        query_expr_ct <- full_query_assay[, query_cells_ct, drop = FALSE]

        total_var_ref <- if(ncol(ref_expr_ct) > 1) sum(apply(as.matrix(ref_expr_ct), 1, stats::var)) else 0
        total_var_query <- if(ncol(query_expr_ct) > 1) sum(apply(as.matrix(query_expr_ct), 1, stats::var)) else 0

        for (pc in pc_subset) {
            pc_name <- paste0("PC", pc)
            loadings_pc <- pca_rotation[, pc, drop = FALSE]

            if (total_var_ref > .Machine[["double.eps"]]) {
                ref_scores <- crossprod(ref_expr_ct, loadings_pc)
                var_scores_ref <- stats::var(as.vector(ref_scores))
                pct_var_ref <- (var_scores_ref / total_var_ref) * 100
                var_explained_list[[length(var_explained_list) + 1]] <- data.frame(
                    pc = pc_name, cell_type = ct, dataset = "Reference",
                    percent_variance = pct_var_ref, stringsAsFactors = FALSE)
            }
            if (total_var_query > .Machine[["double.eps"]]) {
                query_scores <- crossprod(query_expr_ct, loadings_pc)
                var_scores_query <- stats::var(as.vector(query_scores))
                pct_var_query <- (var_scores_query / total_var_query) * 100
                var_explained_list[[length(var_explained_list) + 1]] <- data.frame(
                    pc = pc_name, cell_type = ct, dataset = "Query",
                    percent_variance = pct_var_query, stringsAsFactors = FALSE)
            }
        }
    }
    cell_type_variance_df <- do.call(rbind, var_explained_list)
    rownames(cell_type_variance_df) <- NULL

    # Data for plotting (use all cells, not filtered by anomaly status)
    query_relevant_cells <- unlist(query_cell_indices[available_cell_types], use.names = FALSE)
    ref_relevant_cells <- unlist(ref_cell_indices[available_cell_types], use.names = FALSE)

    query_expr_plot <- SummarizedExperiment::assay(query_data,
                                                   assay_name)[all_top_genes,
                                                               query_relevant_cells, drop = FALSE]
    ref_expr_plot <- SummarizedExperiment::assay(reference_data,
                                                 assay_name)[all_top_genes,
                                                             ref_relevant_cells, drop = FALSE]
    expression_data <- cbind(ref_expr_plot, query_expr_plot)

    # Create cell metadata
    cell_metadata <- rbind(
        data.frame(cell_id = colnames(ref_expr_plot), dataset = "Reference",
                   cell_type = reference_data[[ref_cell_type_col]][ref_relevant_cells],
                   original_index = ref_relevant_cells, stringsAsFactors = FALSE),
        data.frame(cell_id = colnames(query_expr_plot), dataset = "Query",
                   cell_type = query_data[[query_cell_type_col]][query_relevant_cells],
                   original_index = query_relevant_cells, stringsAsFactors = FALSE)
    )

    # Map anomaly status using cell names
    if (!is.null(anomaly_results)) {
        # Initialize anomaly status
        cell_metadata[["anomaly_status"]] <- "Normal"

        # Process anomaly results for each cell type
        for (ct in available_cell_types) {
            if (ct %in% names(anomaly_results)) {

                # Map reference anomaly status using cell names
                if ("reference_anomaly" %in% names(anomaly_results[[ct]])) {
                    ref_anomaly_names <- names(anomaly_results[[ct]][["reference_anomaly"]])
                    ref_anomaly_status <- anomaly_results[[ct]][["reference_anomaly"]]

                    if (!is.null(ref_anomaly_names) && length(ref_anomaly_names) > 0) {
                        # Find which cells are anomalous and strip "Reference_" prefix
                        anomalous_ref_names <- ref_anomaly_names[ref_anomaly_status]
                        anomalous_ref_names_clean <- gsub("^Reference_", "", anomalous_ref_names)

                        ref_mask <- cell_metadata[["dataset"]] == "Reference" &
                            cell_metadata[["cell_id"]] %in% anomalous_ref_names_clean
                        cell_metadata[["anomaly_status"]][ref_mask] <- "Anomaly"
                    }
                }

                # Map query anomaly status using cell names
                if ("query_anomaly" %in% names(anomaly_results[[ct]])) {
                    query_anomaly_names <- names(anomaly_results[[ct]][["query_anomaly"]])
                    query_anomaly_status <- anomaly_results[[ct]][["query_anomaly"]]

                    if (!is.null(query_anomaly_names) && length(query_anomaly_names) > 0) {
                        # Find which cells are anomalous and strip "Query_" prefix
                        anomalous_query_names <- query_anomaly_names[query_anomaly_status]
                        anomalous_query_names_clean <- gsub("^Query_", "", anomalous_query_names)

                        query_mask <- cell_metadata[["dataset"]] == "Query" &
                            cell_metadata[["cell_id"]] %in% anomalous_query_names_clean
                        cell_metadata[["anomaly_status"]][query_mask] <- "Anomaly"
                    }
                }
            }
        }
    }

    # Assemble final results
    final_results <- c(
        pc_results,
        list(
            expression_data = expression_data,
            cell_metadata = cell_metadata,
            gene_metadata = gene_metadata,
            pca_rotation = pca_rotation,
            percent_var = pc_percent_var,
            cell_type_variance = cell_type_variance_df
        )
    )

    # Add analysis type information
    final_results[["analysis_type"]] <- if (anomaly_comparison && !is.null(anomaly_results)) {
        "non_anomalous_reference_vs_anomalous_query"
    } else {
        "all_reference_vs_all_query"
    }

    # Add anomaly results if available
    if (!is.null(anomaly_results)) {
        final_results[["anomaly_results"]] <- anomaly_results
    }

    class(final_results) <- c(class(final_results), "calculateGeneShiftsObject")
    return(final_results)
}

#' @title Process Genes for Statistical Analysis
#'
#' @description
#' Internal helper function that performs statistical testing on top loading genes
#' to compare their expression distributions between query and reference datasets.
#'
#' @param top_genes Character vector of gene names to analyze.
#' @param top_loadings Numeric vector of PC loading values corresponding to the genes.
#' @param query_expr_matrix Numeric matrix of expression values for query data (genes × cells).
#' @param ref_expr_matrix Numeric matrix of expression values for reference data (genes × cells).
#' @param cell_type Character string specifying the cell type being analyzed.
#'
#' @keywords internal
#'
#' @return A list of data frames, where each data frame contains results for one gene.
#'
# Helper function to process genes for statistical analysis
processGenesSimple <- function(top_genes,
                               top_loadings,
                               query_expr_matrix,
                               ref_expr_matrix,
                               cell_type) {

    result_list <- list()

    for (i in seq_along(top_genes)) {
        gene <- top_genes[i]
        loading <- top_loadings[i]

        query_expr <- query_expr_matrix[gene, ]
        ref_expr <- ref_expr_matrix[gene, ]

        query_expr <- query_expr[!is.na(query_expr)]
        ref_expr <- ref_expr[!is.na(ref_expr)]

        if (length(query_expr) < 3 || length(ref_expr) < 3) {
            next
        }

        mean_query <- mean(query_expr)
        mean_reference <- mean(ref_expr)

        test_result <- suppressWarnings(
            stats::wilcox.test(query_expr, ref_expr, alternative = "two.sided")
        )

        result_row <- data.frame(
            gene = gene,
            loading = loading,
            cell_type = cell_type,
            p_value = test_result[["p.value"]],
            mean_query = mean_query,
            mean_reference = mean_reference,
            stringsAsFactors = FALSE
        )

        result_list[[length(result_list) + 1]] <- result_row
    }

    return(result_list)
}
