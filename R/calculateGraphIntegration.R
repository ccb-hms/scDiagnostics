#' @title Calculate Graph Community Integration Diagnostics
#'
#' @description
#' This function performs graph-based community detection to identify annotation inconsistencies
#' by detecting query-only communities, true cross-cell-type mixing patterns, and local
#' annotation inconsistencies based on immediate neighborhood analysis.
#'
#' @details
#' The function performs three types of analysis: (1) Communities containing only query cells,
#' (2) Communities where query cells are mixed with reference cells of different cell types
#' WITHOUT any reference cells of the same type, and (3) Local analysis of each query cell's
#' immediate neighbors to detect annotation inconsistencies even within mixed communities.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param cell_types A character vector specifying the cell types to include in the analysis. If NULL, all cell types are included.
#' @param pc_subset A vector specifying the subset of principal components to use in the analysis. Default is 1:10.
#' @param k_neighbors An integer specifying the number of nearest neighbors for graph construction. Default is 30.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param resolution Resolution parameter for Leiden clustering. Default is 0.15 for fewer, larger communities.
#' @param min_cells_per_community Minimum number of cells required for a community to be analyzed. Default is 10.
#' @param min_cells_per_celltype Minimum number of cells required per cell type for inclusion. Default is 20.
#' @param high_query_prop_threshold Minimum proportion of query cells to consider a community "query-only". Default is 0.9.
#' @param cross_type_threshold Minimum proportion needed to flag cross-cell-type mixing. Default is 0.1.
#' @param local_consistency_threshold Minimum proportion of reference neighbors that should support a query cell's annotation. Default is 0.6.
#' @param local_confidence_threshold Minimum confidence difference needed to suggest re-annotation. Default is 0.2.
#'
#' @return A list containing:
#' \item{high_query_prop_analysis}{Analysis of communities with only query cells}
#' \item{cross_type_mixing}{Analysis of communities with true query-reference cross-cell-type mixing}
#' \item{local_annotation_inconsistencies}{Local neighborhood-based annotation inconsistencies}
#' \item{local_inconsistency_summary}{Summary of local inconsistencies by cell type}
#' \item{community_composition}{Detailed composition of each community}
#' \item{annotation_consistency}{Summary of annotation consistency issues}
#' \item{overall_metrics}{Overall diagnostic metrics}
#' \item{graph_info}{Graph structure information for plotting}
#' \item{parameters}{Analysis parameters used}
#' The list is assigned the class \code{"calculateGraphIntegration"}.
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
#' # Remove a cell type (Myeloid)
#' library(scater)
#' library(SingleR)
#' reference_data <- reference_data[, reference_data$expert_annotation != "Myeloid"]
#' reference_data <- runPCA(reference_data, ncomponents = 50)
#' SingleR_annotation <- SingleR(query_data, reference_data,
#'                               labels = reference_data$expert_annotation)
#' query_data$SingleR_annotation <- SingleR_annotation$labels
#'
#' # Check annotation data
#' table(Expert = query_data$expert_annotation, SingleR = query_data$SingleR_annotation)
#'
#' # Run comprehensive annotation consistency diagnostics
#' graph_diagnostics <- calculateGraphIntegration(
#'     query_data = query_data,
#'     reference_data = reference_data,
#'     query_cell_type_col = "SingleR_annotation",
#'     ref_cell_type_col = "expert_annotation",
#'    pc_subset = 1:10,
#'    k_neighbors = 30,
#'     resolution = 0.1,
#'     high_query_prop_threshold = 0.9,
#'     cross_type_threshold = 0.15,
#'     local_consistency_threshold = 0.6,
#'     local_confidence_threshold = 0.2
#' )
#'
#' # Look at main output
#' graph_diagnostics$overall_metrics
#'
#' # Network graph showing all issue types (color by cell type)
#' plot(graph_diagnostics, plot_type = "community_network", color_by = "cell_type")
#'
#' # Network graph showing all issue types
#' plot(graph_diagnostics, plot_type = "cell_network",
#'      max_nodes = 2000, color_by = "community_type")
#'
#' # Network graph showing all issue types
#' plot(graph_diagnostics, plot_type = "community_data")
#'
#' # Summary bar chart of all issues by cell type
#' plot(graph_diagnostics, plot_type = "summary")
#'
#' # Focus on local annotation inconsistencies
#' plot(graph_diagnostics, plot_type = "local_issues")
#'
#' # Overall annotation issues overview
#' plot(graph_diagnostics, plot_type = "annotation_issues")
#'
#' @importFrom stats dist
#'
# Function to calculate graph integration diagnostics
calculateGraphIntegration <- function(query_data,
                                      reference_data,
                                      query_cell_type_col,
                                      ref_cell_type_col,
                                      cell_types = NULL,
                                      pc_subset = 1:10,
                                      k_neighbors = 30,
                                      assay_name = "logcounts",
                                      resolution = 0.1,
                                      min_cells_per_community = 10,
                                      min_cells_per_celltype = 20,
                                      high_query_prop_threshold = 0.9,
                                      cross_type_threshold = 0.15,
                                      local_consistency_threshold = 0.6,
                                      local_confidence_threshold = 0.2) {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Check additional parameters
    if (!is.numeric(k_neighbors) || k_neighbors <= 0 ||
        k_neighbors != as.integer(k_neighbors)) {
        stop("\'k_neighbors\' must be a positive integer.")
    }

    if (!is.numeric(resolution) || resolution <= 0) {
        stop("\'resolution\' must be a positive number.")
    }

    if (!is.numeric(high_query_prop_threshold) || high_query_prop_threshold <= 0.5 ||
        high_query_prop_threshold > 1) {
        stop("\'high_query_prop_threshold\' must be between 0.5 and 1.")
    }

    if (!is.numeric(cross_type_threshold) || cross_type_threshold <= 0 ||
        cross_type_threshold >= 0.5) {
        stop("\'cross_type_threshold\' must be between 0 and 0.5.")
    }

    if (!is.numeric(local_consistency_threshold) || local_consistency_threshold <= 0 ||
        local_consistency_threshold > 1) {
        stop("\'local_consistency_threshold\' must be between 0 and 1.")
    }

    if (!is.numeric(local_confidence_threshold) || local_confidence_threshold <= 0 ||
        local_confidence_threshold > 1) {
        stop("\'local_confidence_threshold\' must be between 0 and 1.")
    }

    # Get common cell types if not specified
    if (is.null(cell_types)) {
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }

    # Filter cell types by minimum cell count
    ref_counts <- table(reference_data[[ref_cell_type_col]])
    query_counts <- table(query_data[[query_cell_type_col]])

    valid_cell_types <- intersect(
        names(ref_counts)[ref_counts >= min_cells_per_celltype],
        names(query_counts)[query_counts >= min_cells_per_celltype]
    )

    cell_types <- intersect(cell_types, valid_cell_types)

    if (length(cell_types) == 0) {
        stop("No cell types meet the minimum cell count requirement.")
    }

    # Get PCA projection
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             pc_subset = pc_subset,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             assay_name = assay_name)

    # Filter to selected cell types
    pca_filtered <- pca_output[pca_output[["cell_type"]] %in% cell_types, ]

    if (nrow(pca_filtered) == 0) {
        stop("No cells remain after filtering by cell types.")
    }

    # Extract PCA coordinates
    pc_coords <- as.matrix(pca_filtered[, paste0("PC", pc_subset)])

    # Construct k-nearest neighbor graph
    .constructKNNGraph <- function(coords, k) {
        knn_result <- FNN::get.knn(coords, k = k)
        edges <- do.call(rbind, lapply(seq_len(nrow(coords)), function(i) {
            neighbors <- knn_result[["nn.index"]][i, ]
            cbind(i, neighbors)
        }))
        graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
        graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
        return(list(graph = graph, edges = edges))
    }

    # Perform Leiden clustering
    .performLeidenClustering <- function(graph, resolution) {
        clustering <- igraph::cluster_leiden(graph, resolution = resolution)
        membership_vec <- igraph::membership(clustering)
        mod_value <- igraph::modularity(graph, membership_vec)
        return(list(membership = membership_vec, modularity = mod_value))
    }

    # Main analysis function for annotation consistency
    .analyzeAnnotationConsistency <- function(membership, cell_info,
                                              query_threshold, cross_threshold) {
        communities <- unique(membership)

        # Initialize results
        high_query_prop_communities <- list()
        cross_type_communities <- list()
        all_communities <- list()

        for (comm in communities) {
            comm_cells <- which(membership == comm)
            comm_info <- cell_info[comm_cells, ]

            # Skip small communities
            if (nrow(comm_info) < min_cells_per_community) next

            # Basic counts
            dataset_counts <- table(comm_info[["dataset"]])
            celltype_counts <- table(comm_info[["cell_type"]])

            total_cells <- nrow(comm_info)
            n_reference <- ifelse("Reference" %in% names(dataset_counts),
                                  as.numeric(dataset_counts[["Reference"]]), 0)
            n_query <- ifelse("Query" %in% names(dataset_counts),
                              as.numeric(dataset_counts[["Query"]]), 0)

            query_prop <- n_query / total_cells

            # Community summary
            comm_summary <- list(
                community = comm,
                total_cells = total_cells,
                n_reference = n_reference,
                n_query = n_query,
                query_proportion = query_prop,
                cell_types = names(celltype_counts),
                cell_type_counts = as.list(celltype_counts),
                dominant_celltype = names(celltype_counts)[which.max(celltype_counts)]
            )

            all_communities[[length(all_communities) + 1]] <- comm_summary

            # Check for query-only communities
            if (query_prop >= query_threshold) {
                query_analysis <- list(
                    community = comm,
                    total_query_cells = n_query,
                    query_cell_types = names(
                        table(comm_info[["cell_type"]][comm_info[["dataset"]] == "Query"])),
                    query_celltype_counts = table(
                        comm_info[["cell_type"]][comm_info[["dataset"]] == "Query"]),
                    is_homogeneous = length(unique(comm_info[["cell_type"]])) == 1,
                    mixed_types = if(length(unique(comm_info[["cell_type"]])) > 1) {
                        paste(names(table(comm_info[["cell_type"]])), collapse = ", ")
                    } else { NA }
                )
                high_query_prop_communities[[length(high_query_prop_communities) + 1]] <-
                    query_analysis
            }

            # Check for TRUE cross-cell-type mixing only
            if (n_reference > 0 && n_query > 0) {
                # Get cell types for each dataset
                ref_celltypes <- table(
                    comm_info[["cell_type"]][comm_info[["dataset"]] == "Reference"])
                query_celltypes <- table(
                    comm_info[["cell_type"]][comm_info[["dataset"]] == "Query"])

                ref_types <- names(ref_celltypes)
                query_types <- names(query_celltypes)

                # Only flag TRUE mismatches: query cells with NO matching reference type
                cross_mixing_detected <- FALSE
                cross_details <- list()

                for (qt in query_types) {
                    # Check if this query type has NO corresponding reference cells
                    if (!(qt %in% ref_types)) {
                        # This query type is ONLY mixing with different reference types
                        query_count <- as.numeric(query_celltypes[[qt]])
                        query_prop_in_community <- query_count / n_query

                        # Only flag if significant proportion
                        if (query_prop_in_community >= cross_threshold) {
                            cross_mixing_detected <- TRUE

                            # Find which reference types this query type is mixing with
                            ref_mixing_with <- ref_types
                            ref_counts_mixing <- ref_celltypes[ref_mixing_with]

                            cross_details[[paste("Query", qt, "with Ref",
                                                 paste(ref_mixing_with,
                                                       collapse = "+"))]] <- list(
                                                           query_type = qt,
                                                           reference_types_mixing_with = ref_mixing_with,
                                                           reference_counts_mixing_with = ref_counts_mixing,
                                                           query_cells = query_count,
                                                           query_proportion = query_prop_in_community,
                                                           is_true_mismatch = TRUE
                                                       )
                        }
                    }
                }

                if (cross_mixing_detected) {
                    cross_analysis <- list(
                        community = comm,
                        total_cells = total_cells,
                        n_query = n_query,
                        n_reference = n_reference,
                        cross_mixing_details = cross_details,
                        query_types_with_no_ref_match = setdiff(query_types,
                                                                ref_types),
                        reference_types_in_community = ref_types,
                        query_types_in_community = query_types
                    )
                    cross_type_communities[[length(cross_type_communities) + 1]] <-
                        cross_analysis
                }
            }
        }

        return(list(
            high_query_prop = high_query_prop_communities,
            cross_type = cross_type_communities,
            all_communities = all_communities
        ))
    }

    # NEW: Local annotation consistency analysis
    .analyzeLocalAnnotationConsistency <- function(
        membership, cell_info, knn_edges,
        local_threshold, confidence_threshold) {
        # For each query cell, check if its immediate neighbors support its annotation
        query_cells <- which(cell_info[["dataset"]] == "Query")
        local_inconsistencies <- list()

        for (query_idx in query_cells) {
            query_type <- cell_info[["cell_type"]][query_idx]
            query_community <- cell_info[["community"]][query_idx]

            # Find this query cell's direct neighbors from KNN graph
            neighbors_from <- knn_edges[knn_edges[, 1] == query_idx, 2]
            neighbors_to <- knn_edges[knn_edges[, 2] == query_idx, 1]
            all_neighbors <- unique(c(neighbors_from, neighbors_to))

            if (length(all_neighbors) == 0) next

            # Get neighbor information
            neighbor_info <- cell_info[all_neighbors, ]

            # Focus on reference neighbors only
            ref_neighbors <- neighbor_info[neighbor_info[["dataset"]] == "Reference", ]

            if (nrow(ref_neighbors) == 0) next

            # Calculate support for query cell's annotation among reference neighbors
            ref_neighbor_types <- table(ref_neighbors[["cell_type"]])
            total_ref_neighbors <- nrow(ref_neighbors)

            # Support for query's annotation
            support_for_query_type <- ifelse(query_type %in% names(ref_neighbor_types),
                                             as.numeric(ref_neighbor_types[[query_type]]),
                                             0)
            support_proportion <- support_for_query_type / total_ref_neighbors

            # Flag as inconsistent if support is below threshold
            if (support_proportion < local_threshold) {
                # Find the most supported type among reference neighbors
                most_supported_type <- names(ref_neighbor_types)[which.max(ref_neighbor_types)]
                most_supported_count <- max(ref_neighbor_types)
                most_supported_proportion <- most_supported_count / total_ref_neighbors

                # Only flag if there's a clear alternative with stronger support
                if (most_supported_proportion > support_proportion + confidence_threshold) {
                    inconsistency <- list(
                        query_cell_idx = query_idx,
                        query_annotation = query_type,
                        suggested_annotation = most_supported_type,
                        community = query_community,
                        ref_neighbors_total = total_ref_neighbors,
                        support_for_query = support_for_query_type,
                        support_proportion = support_proportion,
                        suggested_support = most_supported_count,
                        suggested_proportion = most_supported_proportion,
                        confidence_score = most_supported_proportion - support_proportion
                    )
                    local_inconsistencies[[length(local_inconsistencies) + 1]] <- inconsistency
                }
            }
        }

        return(local_inconsistencies)
    }

    # Build graph and perform clustering
    if (k_neighbors >= nrow(pc_coords)) {
        k_neighbors <- nrow(pc_coords) - 1
        warning(paste("k_neighbors reduced to", k_neighbors, "due to insufficient cells"))
    }

    graph_result <- .constructKNNGraph(pc_coords, k_neighbors)
    knn_graph <- graph_result[["graph"]]
    clustering_result <- .performLeidenClustering(knn_graph, resolution)

    # Prepare cell information
    cell_info <- data.frame(
        cell_id = seq_len(nrow(pca_filtered)),
        dataset = pca_filtered[["dataset"]],
        cell_type = pca_filtered[["cell_type"]],
        community = clustering_result[["membership"]],
        stringsAsFactors = FALSE
    )

    # Perform main community analysis
    analysis_results <- .analyzeAnnotationConsistency(
        clustering_result[["membership"]],
        cell_info,
        high_query_prop_threshold,
        cross_type_threshold
    )

    # Perform local annotation consistency analysis
    local_inconsistencies <- .analyzeLocalAnnotationConsistency(
        clustering_result[["membership"]],
        cell_info,
        graph_result[["edges"]],
        local_consistency_threshold,
        local_confidence_threshold
    )

    # Convert results to data frames for easier handling
    if (length(analysis_results[["high_query_prop"]]) > 0) {
        high_query_prop_df <- do.call(rbind, lapply(analysis_results[["high_query_prop"]], function(x) {
            data.frame(
                community = x[["community"]],
                total_query_cells = x[["total_query_cells"]],
                query_cell_types = paste(x[["query_cell_types"]], collapse = ", "),
                n_query_cell_types = length(x[["query_cell_types"]]),
                is_homogeneous = x[["is_homogeneous"]],
                mixed_types = ifelse(is.na(x[["mixed_types"]]), NA, x[["mixed_types"]]),
                dominant_query_type = names(x[["query_celltype_counts"]])[which.max(x[["query_celltype_counts"]])],
                stringsAsFactors = FALSE
            )
        }))
    } else {
        high_query_prop_df <- data.frame(
            community = integer(0), total_query_cells = integer(0),
            query_cell_types = character(0), n_query_cell_types = integer(0),
            is_homogeneous = logical(0), mixed_types = character(0),
            dominant_query_type = character(0)
        )
    }

    if (length(analysis_results[["cross_type"]]) > 0) {
        cross_type_df <- do.call(rbind, lapply(analysis_results[["cross_type"]], function(x) {
            # Summarize cross-mixing details - only TRUE mismatches
            cross_summary <- sapply(
                x[["cross_mixing_details"]], function(detail) {
                    ref_summary <- paste(
                        names(detail[["reference_counts_mixing_with"]]),
                        paste0("(", detail[["reference_counts_mixing_with"]], ")"),
                        collapse = "+")
                    paste0("Query ", detail[["query_type"]], " (",
                           detail[["query_cells"]], ") with Ref ",
                           ref_summary)
                })

            data.frame(
                community = x[["community"]],
                total_cells = x[["total_cells"]],
                n_query = x[["n_query"]],
                n_reference = x[["n_reference"]],
                query_types_no_ref_match = paste(x[["query_types_with_no_ref_match"]],
                                                 collapse = ", "),
                reference_types_available = paste(x[["reference_types_in_community"]],
                                                  collapse = ", "),
                cross_mixing_summary = paste(cross_summary, collapse = "; "),
                n_true_mismatches = length(x[["cross_mixing_details"]]),
                stringsAsFactors = FALSE
            )
        }))
    } else {
        cross_type_df <- data.frame(
            community = integer(0),
            total_cells = integer(0),
            n_query = integer(0),
            n_reference = integer(0),
            query_types_no_ref_match = character(0),
            reference_types_available = character(0),
            cross_mixing_summary = character(0),
            n_true_mismatches = integer(0)
        )
    }

    # Convert local inconsistencies to data frame
    if (length(local_inconsistencies) > 0) {
        local_inconsistency_df <- do.call(rbind, lapply(local_inconsistencies, function(x) {
            data.frame(
                query_cell_idx = x[["query_cell_idx"]],
                current_annotation = x[["query_annotation"]],
                suggested_annotation = x[["suggested_annotation"]],
                community = x[["community"]],
                ref_neighbors_total = x[["ref_neighbors_total"]],
                current_support = x[["support_for_query"]],
                current_support_prop = x[["support_proportion"]],
                suggested_support = x[["suggested_support"]],
                suggested_support_prop = x[["suggested_proportion"]],
                confidence_score = x[["confidence_score"]],
                stringsAsFactors = FALSE
            )
        }))

        # Add cell type information
        local_inconsistency_df[["cell_type"]] <-
            cell_info[["cell_type"]][local_inconsistency_df[["query_cell_idx"]]]
    } else {
        local_inconsistency_df <- data.frame(
            query_cell_idx = integer(0), current_annotation = character(0),
            suggested_annotation = character(0), community = integer(0),
            ref_neighbors_total = integer(0), current_support = integer(0),
            current_support_prop = numeric(0), suggested_support = integer(0),
            suggested_support_prop = numeric(0), confidence_score = numeric(0),
            cell_type = character(0)
        )
    }

    # Create community composition summary
    community_composition <-
        do.call(rbind, lapply(analysis_results[["all_communities"]], function(x) {
            data.frame(
                community = x[["community"]],
                total_cells = x[["total_cells"]],
                n_reference = x[["n_reference"]],
                n_query = x[["n_query"]],
                query_proportion = x[["query_proportion"]],
                cell_types = paste(x[["cell_types"]], collapse = ", "),
                n_cell_types = length(x[["cell_types"]]),
                dominant_celltype = x[["dominant_celltype"]],
                is_high_query_prop = x[["query_proportion"]] >= high_query_prop_threshold,
                has_true_cross_mixing = x[["community"]] %in% cross_type_df[["community"]],
                stringsAsFactors = FALSE
            )
        }))

    # Calculate annotation consistency summary
    annotation_summary <- data.frame(
        cell_type = cell_types,
        total_query_cells = sapply(
            cell_types, function(ct) sum(cell_info[["cell_type"]] == ct &
                                             cell_info[["dataset"]] == "Query")),
        high_query_prop_cells = sapply(cell_types, function(ct) {
            if (nrow(high_query_prop_df) == 0) return(0)
            sum(sapply(seq_len(nrow(high_query_prop_df)), function(i) {
                comm_id <- high_query_prop_df[["community"]][i]
                comm_cells <- cell_info[cell_info[["community"]] == comm_id, ]
                sum(comm_cells[["cell_type"]] == ct & comm_cells[["dataset"]] == "Query")
            }))
        }),
        true_cross_mixing_cells = sapply(cell_types, function(ct) {
            if (nrow(cross_type_df) == 0) return(0)
            sum(sapply(seq_len(nrow(cross_type_df)), function(i) {
                comm_id <- cross_type_df[["community"]][i]
                comm_cells <- cell_info[cell_info[["community"]] == comm_id, ]
                sum(comm_cells[["cell_type"]] == ct & comm_cells[["dataset"]] == "Query")
            }))
        }),
        stringsAsFactors = FALSE
    )

    annotation_summary[["query_isolation_rate"]] <-
        annotation_summary[["high_query_prop_cells"]] /
        annotation_summary[["total_query_cells"]]
    annotation_summary[["true_cross_mixing_rate"]] <-
        annotation_summary[["true_cross_mixing_cells"]] /
        annotation_summary[["total_query_cells"]]

    # Replace NaN with 0
    annotation_summary[["query_isolation_rate"]][
        is.nan(annotation_summary[["query_isolation_rate"]])] <- 0
    annotation_summary[["true_cross_mixing_rate"]][
        is.nan(annotation_summary[["true_cross_mixing_rate"]])] <- 0

    # Add local inconsistency summary by cell type
    local_inconsistency_summary <- data.frame(
        cell_type = cell_types,
        total_query_cells = sapply(cell_types, function(ct) sum(cell_info[["cell_type"]] == ct &
                                                                    cell_info[["dataset"]] == "Query")),
        locally_inconsistent_cells = sapply(cell_types, function(ct) {
            sum(local_inconsistency_df[["current_annotation"]] == ct)
        }),
        stringsAsFactors = FALSE
    )
    local_inconsistency_summary[["local_inconsistency_rate"]] <-
        local_inconsistency_summary[["locally_inconsistent_cells"]] /
        local_inconsistency_summary[["total_query_cells"]]
    local_inconsistency_summary[["local_inconsistency_rate"]][
        is.nan(local_inconsistency_summary[["local_inconsistency_rate"]])] <- 0

    # Calculate layout for graph plotting
    graph_layout <- igraph::layout_with_fr(knn_graph, dim = 2)
    graph_info <- list(layout = graph_layout, edges = graph_result[["edges"]],
                       n_nodes = nrow(pc_coords))

    # Overall metrics
    overall_metrics <- list(
        total_communities = nrow(community_composition),
        high_query_prop_communities = nrow(high_query_prop_df),
        true_cross_type_communities = nrow(cross_type_df),
        total_high_query_prop_cells = sum(high_query_prop_df[["total_query_cells"]]),
        total_true_cross_mixing_cells = sum(cross_type_df[["n_query"]]),
        total_locally_inconsistent_cells = nrow(local_inconsistency_df),
        modularity = clustering_result[["modularity"]],
        mean_query_isolation_rate = mean(annotation_summary[["query_isolation_rate"]]),
        mean_true_cross_mixing_rate = mean(annotation_summary[["true_cross_mixing_rate"]]),
        mean_local_inconsistency_rate = mean(local_inconsistency_summary[["local_inconsistency_rate"]])
    )

    # Parameters
    parameters <- list(
        pc_subset = pc_subset,
        k_neighbors = k_neighbors,
        resolution = resolution,
        min_cells_per_community = min_cells_per_community,
        min_cells_per_celltype = min_cells_per_celltype,
        high_query_prop_threshold = high_query_prop_threshold,
        cross_type_threshold = cross_type_threshold,
        local_consistency_threshold = local_consistency_threshold,
        local_confidence_threshold = local_confidence_threshold,
        cell_types_analyzed = cell_types
    )

    # Create output
    result <- list(
        high_query_prop_analysis = high_query_prop_df,
        cross_type_mixing = cross_type_df,
        local_annotation_inconsistencies = local_inconsistency_df,
        local_inconsistency_summary = local_inconsistency_summary,
        community_composition = community_composition,
        annotation_consistency = annotation_summary,
        overall_metrics = overall_metrics,
        graph_info = graph_info,
        parameters = parameters,
        cell_info = cell_info
    )

    class(result) <- c(class(result), "calculateGraphIntegrationObject")
    return(result)
}
