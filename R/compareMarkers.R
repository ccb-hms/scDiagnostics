#' @title Compare Marker Gene Expression between Query and Reference Data
#'
#' @description
#' This function identifies marker genes for each cell type in both query and reference datasets
#' using the standard Bioconductor approach (Wilcoxon rank-sum test), and compares their expression
#' patterns to assess annotation quality. It can optionally filter query cells based on anomaly detection
#' results and restrict analysis to specific cell types.
#'
#' @details
#' The function performs the following steps:
#' 1. Optionally performs anomaly detection and filters query cells based on results.
#' 2. Identifies marker genes for each cell type in both datasets using \code{findMarkers} approach.
#' 3. Reference markers are always computed using all reference cells for each cell type.
#' 4. Query markers are computed using filtered cells (anomalous/non-anomalous) if specified.
#' 5. Compares the overlap of top marker genes between corresponding cell types.
#' 6. Evaluates the expression consistency of reference markers in query data.
#' 7. Provides quality scores based on marker gene concordance.
#'
#' Marker genes are identified using Wilcoxon rank-sum tests comparing each cell type against all others.
#' High overlap and consistent expression of markers indicate good annotation quality.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are
#'                   included.
#' @param n_markers Number of top marker genes to consider for each cell type. Default is 50.
#' @param min_cells Minimum number of cells required per cell type for marker identification. Default is 10.
#' @param anomaly_filter Character string specifying how to filter query cells based on anomaly detection.
#'                       Options: "none" (default), "anomalous_only", "non_anomalous_only".
#' @param assay_name Name of the assay to use for computations. Default is "logcounts".
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#' @param ... Additional arguments passed to the \code{detectAnomaly} function.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{marker_overlap}{Matrix showing overlap of top markers between query and reference for each cell type.}
#'   \item{expression_consistency}{Matrix showing expression consistency of reference markers in query data.}
#'   \item{quality_scores}{Named vector of quality assessments for each cell type.}
#'   \item{markers_query}{List of marker gene results for each cell type in query data.}
#'   \item{markers_ref}{List of marker gene results for each cell type in reference data.}
#'   \item{common_cell_types}{Vector of cell types present in both datasets.}
#'   \item{n_cells_query}{Named vector of cell counts per type in query data.}
#'   \item{n_cells_ref}{Named vector of cell counts per type in reference data.}
#'   \item{anomaly_filter_used}{Character string indicating the anomaly filter applied.}
#'   \item{selected_cell_types}{Character vector of cell types analyzed.}
#'   \item{anomaly_output}{Output from anomaly detection if performed.}
#' }
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.compareMarkersObject}}, \code{\link{detectAnomaly}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Compare marker genes
#' marker_comparison <- compareMarkers(query_data = query_data,
#'                                     reference_data = reference_data,
#'                                     query_cell_type_col = "expert_annotation",
#'                                     ref_cell_type_col = "expert_annotation")
#'
#' # With anomaly filtering
#' marker_comparison_filtered <- compareMarkers(query_data = query_data,
#'                                             reference_data = reference_data,
#'                                             query_cell_type_col = "expert_annotation",
#'                                             ref_cell_type_col = "expert_annotation",
#'                                             anomaly_filter = "non_anomalous_only")
#'
#' # Visualize results
#' plot(marker_comparison)
#'
#' @importFrom stats var wilcox.test
#' @importFrom utils head
#'
# Function to compare marker genes between query and reference datasets
compareMarkers <- function(query_data,
                           reference_data,
                           query_cell_type_col,
                           ref_cell_type_col,
                           cell_types = NULL,
                           n_markers = 50,
                           min_cells = 10,
                           anomaly_filter = c("none", "anomalous_only", "non_anomalous_only"),
                           assay_name = "logcounts",
                           max_cells_query = 5000,
                           max_cells_ref = 5000,
                           ...){

    # Match arguments
    anomaly_filter <- match.arg(anomaly_filter)

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
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

    # Ensure cell names exist for anomaly detection mapping
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

    # Perform anomaly detection FIRST (before any downsampling) if needed
    anomaly_output <- NULL
    if (anomaly_filter != "none") {
        tryCatch({
            anomaly_output <- detectAnomaly(reference_data = reference_data,
                                            query_data = query_data,
                                            ref_cell_type_col = ref_cell_type_col,
                                            query_cell_type_col = query_cell_type_col,
                                            cell_types = cell_types,
                                            max_cells_ref = NULL,      # No downsampling in anomaly detection
                                            max_cells_query = NULL,    # No downsampling in anomaly detection
                                            ...)
        }, error = function(e) {
            warning("Anomaly detection failed: ", e[["message"]], ". Proceeding without anomaly filtering.")
            anomaly_output <<- NULL
        })
    }

    # Downsample query and reference data (with cell type filtering)
    query_data <- downsampleSCE(sce_object = query_data,
                                max_cells = max_cells_query,
                                cell_types = cell_types,
                                cell_type_col = query_cell_type_col)
    reference_data <- downsampleSCE(sce_object = reference_data,
                                    max_cells = max_cells_ref,
                                    cell_types = cell_types,
                                    cell_type_col = ref_cell_type_col)

    # Get cell type information after downsampling
    query_cell_types_orig <- colData(query_data)[[query_cell_type_col]]
    ref_cell_types_orig <- colData(reference_data)[[ref_cell_type_col]]

    # Reference data is NEVER filtered - always use all reference cells
    reference_data_filtered <- reference_data
    ref_cell_types <- ref_cell_types_orig

    # Query data filtering based on anomaly detection results
    query_data_filtered <- query_data
    query_cell_types <- query_cell_types_orig

    if (anomaly_filter != "none" && !is.null(anomaly_output)) {

        # Function to filter ONLY query cells based on anomaly results using cell names
        .filterQueryCellsByAnomalyNames <- function(data, anomaly_data, filter_type) {
            current_cell_names <- colnames(data)
            cells_to_keep <- rep(TRUE, ncol(data))
            current_cell_types <- colData(data)[[query_cell_type_col]]

            for (cell_type in unique(current_cell_types)) {
                if (cell_type %in% names(anomaly_data)) {
                    # Get indices of cells of this type in current (downsampled) data
                    cell_indices <- which(current_cell_types == cell_type)
                    cell_names_this_type <- current_cell_names[cell_indices]

                    if ("query_anomaly" %in% names(anomaly_data[[cell_type]])) {
                        # Get anomaly status from original detection
                        anomaly_names <- names(anomaly_data[[cell_type]][["query_anomaly"]])
                        anomaly_status <- anomaly_data[[cell_type]][["query_anomaly"]]

                        if (!is.null(anomaly_names) && length(anomaly_names) > 0) {
                            # Map current cells to their anomaly status
                            for (i in seq_along(cell_indices)) {
                                cell_name <- cell_names_this_type[i]
                                cell_idx <- cell_indices[i]

                                # Find this cell in the anomaly results
                                anomaly_idx <- match(cell_name, anomaly_names)
                                if (!is.na(anomaly_idx)) {
                                    is_anomalous <- anomaly_status[anomaly_idx]

                                    if (filter_type == "anomalous_only") {
                                        cells_to_keep[cell_idx] <- is_anomalous
                                    } else if (filter_type == "non_anomalous_only") {
                                        cells_to_keep[cell_idx] <- !is_anomalous
                                    }
                                } else {
                                    # Cell not found in anomaly results (shouldn't happen)
                                    # Keep it by default
                                    if (filter_type == "anomalous_only") {
                                        cells_to_keep[cell_idx] <- FALSE
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return(cells_to_keep)
        }

        # Filter ONLY query data using cell names
        query_keep <- .filterQueryCellsByAnomalyNames(query_data_filtered, anomaly_output, anomaly_filter)

        if (sum(query_keep) == 0) {
            warning("Anomaly filtering removed all query cells. Proceeding without filtering.")
        } else {
            query_data_filtered <- query_data_filtered[, query_keep]
            query_cell_types <- query_cell_types[query_keep]
        }
    }

    # Get unique cell types with sufficient cells
    query_types <- names(table(query_cell_types))[table(query_cell_types) >= min_cells]
    ref_types <- names(table(ref_cell_types))[table(ref_cell_types) >= min_cells]
    common_cell_types <- intersect(query_types, ref_types)

    # Filter by user-specified cell types if provided
    if (!is.null(cell_types)) {
        common_cell_types <- intersect(common_cell_types, cell_types)
    }

    if (length(common_cell_types) == 0) {
        stop("No common cell types with sufficient cells found between query and reference data")
    }

    # Function to find markers using Wilcoxon test (standard Bioconductor approach)
    .findMarkers <- function(expr_matrix, cell_types, target_type) {
        target_cells <- cell_types == target_type
        other_cells <- cell_types != target_type

        if (sum(target_cells) < 3 || sum(other_cells) < 3) {
            return(data.frame(gene = character(0), pval = numeric(0),
                              logFC = numeric(0), stringsAsFactors = FALSE))
        }

        # Perform Wilcoxon test for each gene
        pvals <- rep(1, nrow(expr_matrix))
        logFCs <- rep(0, nrow(expr_matrix))

        for (i in seq_len(nrow(expr_matrix))) {
            target_expr <- expr_matrix[i, target_cells]
            other_expr <- expr_matrix[i, other_cells]

            # Skip if no variation
            if (var(target_expr) == 0 && var(other_expr) == 0) {
                next
            }

            # Wilcoxon test with warning suppression
            test_result <- tryCatch({
                wilcox.test(target_expr, other_expr, alternative = "greater", exact = FALSE)
            }, error = function(e) {
                list(p.value = 1)
            })

            pvals[i] <- test_result[["p.value"]]
            logFCs[i] <- mean(target_expr) - mean(other_expr)
        }

        # Adjust p-values
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Create results data frame
        results <- data.frame(
            gene = rownames(expr_matrix),
            pval = pvals,
            adj_pval = adj_pvals,
            logFC = logFCs,
            stringsAsFactors = FALSE
        )

        # Filter and sort
        results <- results[results[["adj_pval"]] < 0.05 & results[["logFC"]] > 0, ]
        results <- results[order(results[["adj_pval"]], decreasing = FALSE), ]

        return(results)
    }

    # Extract expression matrices
    query_matrix <- as.matrix(assay(query_data_filtered, assay_name))
    ref_matrix <- as.matrix(assay(reference_data_filtered, assay_name))  # Uses ALL reference cells

    # Find markers for each cell type
    markers_query <- list()
    markers_ref <- list()

    for (cell_type in common_cell_types) {
        # Query markers: use filtered cells (anomalous/non-anomalous if specified)
        markers_query[[cell_type]] <- .findMarkers(query_matrix, query_cell_types, cell_type)

        # Reference markers: ALWAYS use all reference cells
        markers_ref[[cell_type]] <- .findMarkers(ref_matrix, ref_cell_types, cell_type)
    }

    # Calculate marker overlap
    marker_overlap <- rep(0, length(common_cell_types))
    names(marker_overlap) <- common_cell_types

    # Calculate expression consistency
    expression_consistency <- rep(0, length(common_cell_types))
    names(expression_consistency) <- common_cell_types

    # Cell counts - convert table to named numeric vector
    query_table <- table(query_cell_types)  # Filtered query counts
    ref_table <- table(ref_cell_types)      # All reference counts

    n_cells_query <- as.numeric(query_table[common_cell_types])
    names(n_cells_query) <- common_cell_types

    n_cells_ref <- as.numeric(ref_table[common_cell_types])
    names(n_cells_ref) <- common_cell_types

    for (cell_type in common_cell_types) {
        # Get top markers
        query_top <- head(markers_query[[cell_type]][["gene"]], n_markers)
        ref_top <- head(markers_ref[[cell_type]][["gene"]], n_markers)

        # Calculate overlap (Jaccard index)
        if (length(query_top) > 0 && length(ref_top) > 0) {
            overlap <- length(intersect(query_top, ref_top))
            union_size <- length(union(query_top, ref_top))
            marker_overlap[cell_type] <- overlap / union_size
        }

        # Calculate expression consistency of reference markers in query
        if (length(ref_top) > 0) {
            ref_markers_available <- intersect(ref_top, rownames(query_matrix))
            if (length(ref_markers_available) > 0) {
                target_cells <- query_cell_types == cell_type
                other_cells <- query_cell_types != cell_type

                # Check if reference markers are still upregulated in query
                consistent_markers <- 0
                for (marker in ref_markers_available) {
                    target_mean <- mean(query_matrix[marker, target_cells])
                    other_mean <- mean(query_matrix[marker, other_cells])
                    if (target_mean > other_mean) {
                        consistent_markers <- consistent_markers + 1
                    }
                }
                expression_consistency[cell_type] <- consistent_markers / length(ref_markers_available)
            }
        }
    }

    # Quality assessment based on MINIMUM of both metrics (most conservative)
    overlap_quality <- rep("Poor", length(common_cell_types))
    overlap_quality[marker_overlap >= 0.7] <- "Good"
    overlap_quality[marker_overlap >= 0.4 & marker_overlap < 0.7] <- "Moderate"

    consistency_quality <- rep("Poor", length(common_cell_types))
    consistency_quality[expression_consistency >= 0.7] <- "Good"
    consistency_quality[expression_consistency >= 0.4 & expression_consistency < 0.7] <- "Moderate"

    # Overall quality is the WORST of the two
    quality_scores <- rep("Poor", length(common_cell_types))
    names(quality_scores) <- common_cell_types

    for (i in seq_along(common_cell_types)) {
        overlap_qual <- overlap_quality[i]
        consistency_qual <- consistency_quality[i]

        # Take the worst quality (Poor > Moderate > Good in terms of "badness")
        if (overlap_qual == "Poor" || consistency_qual == "Poor") {
            quality_scores[i] <- "Poor"
        } else if (overlap_qual == "Moderate" || consistency_qual == "Moderate") {
            quality_scores[i] <- "Moderate"
        } else {
            quality_scores[i] <- "Good"
        }
    }

    # Prepare output
    output <- list(
        marker_overlap = marker_overlap,
        expression_consistency = expression_consistency,
        quality_scores = quality_scores,
        markers_query = markers_query,
        markers_ref = markers_ref,
        common_cell_types = common_cell_types,
        n_cells_query = n_cells_query,
        n_cells_ref = n_cells_ref,
        anomaly_filter_used = anomaly_filter,
        selected_cell_types = common_cell_types,
        anomaly_output = anomaly_output
    )

    class(output) <- c(class(output), "compareMarkersObject")
    return(output)
}
