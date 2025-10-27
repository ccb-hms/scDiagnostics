#' @title Plot Graph-Based Integration Diagnostics
#'
#' @description
#' The S3 plot method generates visualizations of annotation consistency diagnostics,
#' including query-only communities, cross-cell-type mixing, and local annotation inconsistencies.
#'
#' @details
#' The S3 plot method creates optimized visualizations showing different types of annotation
#' issues including community-level and local neighborhood-level inconsistencies.
#'
#' @param x An object of class \code{calculateGraphIntegrationObject} containing the diagnostic results.
#' @param plot_type Character string specifying visualization type. Options: "community_network" (default),
#'                  "cell_network", "community_data", "summary", "local_issues", or "annotation_issues".
#' @param color_by Character string specifying the variable to use for coloring points/elements if `plot_type` is
#'                 "community_network" or "cell_network". Default is "cell_type".
#' @param max_nodes Maximum number of nodes to display for performance. Default is 2000.
#' @param point_size Point size for graph nodes. Default is 0.8.
#' @param exclude_reference_only Logical indicating whether to exclude reference-only communities/cells from visualization.
#'                               Default is FALSE.
#' @param ... Additional arguments passed to ggplot2 functions.
#'
#' @return A \code{ggplot} object showing integration diagnostics.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateGraphIntegration}}
#'
#' @rdname calculateGraphIntegration
#'
#' @importFrom stats reorder
#'
# Function to plot the graph integration diagnostics
plot.calculateGraphIntegrationObject <- function(
        x,
        plot_type = c("community_network", "cell_network", "community_data",
                      "summary", "local_issues", "annotation_issues"),
        color_by = c("cell_type", "community_type"),
        max_nodes = 2000,
        point_size = 0.8,
        exclude_reference_only = FALSE,
        ...) {

    plot_type <- match.arg(plot_type)
    color_by <- match.arg(color_by)

    colors <- list(
        high_query_prop = "#DC143C",
        cross_mixing = "#FF8C00",
        local_inconsistent = "#9932CC",
        well_integrated = "#2E8B57",
        reference_only = "#4169E1",
        small_community = "gray70",
        excellent = "#228B22",
        concerning = "#FF6347",
        problematic = "#DC143C"
    )

    pluralize <- function(count, singular, plural = paste0(singular, "s")) {
        if (count == 1) return(paste(count, singular))
        else return(paste(count, plural))
    }

    if (plot_type == "community_network") {
        community_composition <- x[["community_composition"]]

        if (nrow(community_composition) == 0) {
            p <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5,
                                  label = "No communities to display",
                                  size = 8, color = "gray50") +
                ggplot2::labs(title = "Community Network") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 14,
                                                                  face = "bold",
                                                                  hjust = 0.5))
        } else {
            # Filter out reference-only communities if requested
            if (exclude_reference_only) {
                reference_only_mask <- community_composition[["query_proportion"]] < 0.1
                community_composition <- community_composition[!reference_only_mask, ]
            }

            if (nrow(community_composition) == 0) {
                p <- ggplot2::ggplot() +
                    ggplot2::annotate("text", x = 0.5, y = 0.5,
                                      label = "No communities to display after filtering",
                                      size = 8, color = "gray50") +
                    ggplot2::labs(title = "Community Network") +
                    ggplot2::theme_void() +
                    ggplot2::theme(plot.title = ggplot2::element_text(size = 14,
                                                                      face = "bold",
                                                                      hjust = 0.5))
            } else {
                # Calculate community centroids in PCA space
                community_centroids <- data.frame(
                    community = integer(0),
                    x = numeric(0),
                    y = numeric(0),
                    stringsAsFactors = FALSE
                )

                for (comm in community_composition[["community"]]) {
                    comm_cells <- x[["cell_info"]][["community"]] == comm
                    if (sum(comm_cells) > 0) {
                        centroid_x <- mean(x[["graph_info"]][["layout"]][comm_cells, 1])
                        centroid_y <- mean(x[["graph_info"]][["layout"]][comm_cells, 2])
                        community_centroids <- rbind(community_centroids, data.frame(
                            community = comm,
                            x = centroid_x,
                            y = centroid_y,
                            stringsAsFactors = FALSE
                        ))
                    }
                }

                # Calculate adaptive threshold based on overall graph connectivity
                original_edges <- x[["graph_info"]][["edges"]]
                cell_communities <- x[["cell_info"]][["community"]]
                total_cells <- length(cell_communities)
                total_edges <- nrow(original_edges)

                # Base connectivity rate in the graph
                base_connectivity <- total_edges / choose(total_cells, 2)

                # Adaptive threshold: 2x the base rate, with bounds
                adaptive_threshold <- max(0.005, min(0.05, 2 * base_connectivity))

                # Create edge lookup for fast checking
                edge_lookup <- paste(original_edges[, 1],
                                     original_edges[, 2],
                                     sep = "_")
                edge_lookup <-
                    c(edge_lookup, paste(original_edges[, 2],
                                         original_edges[, 1],
                                         sep = "_"))
                edge_set <-
                    as.environment(as.list(setNames(rep(TRUE,
                                                        length(edge_lookup)),
                                                    edge_lookup)))

                edge_data <- data.frame(x = numeric(0), y = numeric(0),
                                        xend = numeric(0), yend = numeric(0),
                                        connection_strength = numeric(0))

                communities <- community_composition[["community"]]
                n_communities <- length(communities)

                # Check all pairs of communities
                if (n_communities > 1) {
                    for (i in seq_len(n_communities - 1)) {
                        for (j in seq(i + 1, n_communities)) {
                            comm1 <- communities[i]
                            comm2 <- communities[j]

                            cells_comm1 <- which(cell_communities == comm1)
                            cells_comm2 <- which(cell_communities == comm2)

                            if (length(cells_comm1) > 0 && length(cells_comm2) > 0) {
                                # Adaptive sampling based on community sizes
                                max_samples <- min(500, length(cells_comm1) * length(cells_comm2))
                                n_samples <- max(100, max_samples)

                                sampled_comm1 <- sample(cells_comm1, n_samples, replace = TRUE)
                                sampled_comm2 <- sample(cells_comm2, n_samples, replace = TRUE)

                                # Count connections
                                connections <- 0
                                for (k in seq_len(n_samples)) {
                                    pair_key <- paste(sampled_comm1[k], sampled_comm2[k], sep = "_")
                                    if (exists(pair_key, envir = edge_set)) {
                                        connections <- connections + 1
                                    }
                                }

                                connection_rate <- connections / n_samples

                                # Use adaptive threshold
                                if (connection_rate >= adaptive_threshold) {
                                    comm1_idx <- which(community_centroids$community == comm1)
                                    comm2_idx <- which(community_centroids$community == comm2)

                                    if (length(comm1_idx) > 0 && length(comm2_idx) > 0) {
                                        edge_data <- rbind(edge_data, data.frame(
                                            x = community_centroids$x[comm1_idx],
                                            y = community_centroids$y[comm1_idx],
                                            xend = community_centroids$x[comm2_idx],
                                            yend = community_centroids$y[comm2_idx],
                                            connection_strength = connection_rate
                                        ))
                                    }
                                }
                            }
                        }
                    }
                }

                # Prepare node data
                node_data <- merge(community_centroids,
                                   community_composition,
                                   by = "community")

                if (color_by == "cell_type") {
                    # Get all cell types and generate paired colors
                    all_cell_types <- sort(unique(x[["cell_info"]][["cell_type"]]))
                    cell_types_cases <- c()
                    for(cell_type_case in all_cell_types){
                        cell_types_cases <- c(cell_types_cases,
                                              paste(cell_type_case, "(Mixed)"),
                                              paste(cell_type_case, "(Pure)"))
                    }
                    paired_colors <- generateColors(cell_types_cases, paired = TRUE)

                    # Determine dominant cell type and whether community is mixed for each community
                    node_data[["dominant_cell_type"]] <- sapply(node_data[["community"]], function(comm) {
                        comm_cells <- x[["cell_info"]][x[["cell_info"]][["community"]] == comm, ]
                        cell_type_counts <- table(comm_cells[["cell_type"]])
                        names(cell_type_counts)[which.max(cell_type_counts)]
                    })
                    node_data[["is_mixed"]] <- node_data[["n_cell_types"]] > 1
                    node_data[["cell_type_case"]] <-
                        paste(node_data[["dominant_cell_type"]],
                              ifelse(node_data[["is_mixed"]], "(Mixed)", "(Pure)"))

                    # Assign colors based on dominant cell type and mixing status
                    node_data[["color"]] <- paired_colors[node_data[["cell_type_case"]]]

                    # Classify communities for shapes
                    node_data[["shape_type"]] <- "Well Integrated"
                    node_data[["shape_type"]][node_data[["community"]] %in%
                                                  x[["cross_type_mixing"]][["community"]]] <-
                        "Cross-Type Mixing"
                    node_data[["shape_type"]][node_data[["community"]] %in%
                                                  x[["high_query_prop_analysis"]][["community"]]] <-
                        "High Query Proportion"
                    node_data[["shape_type"]][node_data[["query_proportion"]] < 0.1 &
                                                  node_data[["shape_type"]] == "Well Integrated"] <-
                        "High Reference Proportion"

                    # Create color legend - only show colors that are actually used
                    used_types <- unique(node_data[["dominant_cell_type"]])
                    used_mixed <- any(node_data[["is_mixed"]])
                    used_pure <- any(!node_data[["is_mixed"]])

                    color_legend_labels <- names(paired_colors)
                    color_legend_values <- paired_colors

                    p <- ggplot2::ggplot() +
                        {if (nrow(edge_data) > 0) {
                            ggplot2::geom_segment(data = edge_data,
                                                  ggplot2::aes(
                                                      x = .data[["x"]], y = .data[["y"]],
                                                      xend = .data[["xend"]], yend = .data[["yend"]],
                                                      alpha = .data[["connection_strength"]]),
                                                  color = "gray60", linewidth = 0.8)
                        }} +
                        ggplot2::geom_point(data = node_data,
                                            ggplot2::aes(
                                                x = .data[["x"]], y = .data[["y"]],
                                                color = .data[["cell_type_case"]],
                                                shape = .data[["shape_type"]],
                                                size = .data[["total_cells"]]),
                                            alpha = 0.8) +
                        ggplot2::scale_color_manual(
                            values = color_legend_values,
                            name = "Dominant Cell Type"
                        ) +
                        ggplot2::scale_shape_manual(
                            values = c("Well Integrated" = 16,
                                       "Cross-Type Mixing" = 18,
                                       "High Query Proportion" = 17,
                                       "High Reference Proportion" = 15),
                            name = "Community Classification"
                        ) +
                        ggplot2::scale_size_continuous(
                            range = c(3, 10),
                            trans = "sqrt",
                            guide = "none"
                        ) +
                        ggplot2::scale_alpha_continuous(
                            range = c(0.3, 1.0),
                            guide = "none"
                        ) +
                        ggplot2::geom_text(
                            data = node_data,
                            ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                         label = .data[["community"]]),
                            color = "black", size = 3, fontface = "bold") +
                        ggplot2::labs(
                            title = "Community Network: Cell Type Composition",
                            subtitle = paste0("Showing ", nrow(node_data),
                                              " communities with ",
                                              nrow(edge_data),
                                              " connections"),
                            caption = "Node size proportional to # cells; Shape indicates community classification; Edge opacity proportional to connection strength"
                        ) +
                        ggplot2::theme_void() +
                        ggplot2::theme(
                            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                            plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
                            legend.position = "right"
                        )
                } else if (color_by == "community_type") {
                    # Original community type classification
                    node_data[["issue_type"]] <- "Well Integrated"
                    node_data[["issue_type"]][node_data[["community"]] %in%
                                                  x[["cross_type_mixing"]][["community"]]] <-
                        "Cross-Type Mixing"
                    node_data[["issue_type"]][node_data[["community"]] %in%
                                                  x[["high_query_prop_analysis"]][["community"]]] <-
                        "High Query Proportion"
                    node_data[["issue_type"]][node_data[["query_proportion"]] < 0.1 &
                                                  node_data[["issue_type"]] == "Well Integrated"] <-
                        "High Reference Proportion"

                    p <- ggplot2::ggplot() +
                        {if (nrow(edge_data) > 0) {
                            ggplot2::geom_segment(data = edge_data,
                                                  ggplot2::aes(
                                                      x = .data[["x"]], y = .data[["y"]],
                                                      xend = .data[["xend"]], yend = .data[["yend"]],
                                                      alpha = .data[["connection_strength"]]),
                                                  color = "gray60", linewidth = 0.8)
                        }} +
                        ggplot2::geom_point(data = node_data,
                                            ggplot2::aes(
                                                x = .data[["x"]], y = .data[["y"]],
                                                color = .data[["issue_type"]],
                                                size = .data[["total_cells"]]),
                                            alpha = 0.8) +
                        ggplot2::scale_color_manual(
                            values = c("High Query Proportion" = colors[["high_query_prop"]],
                                       "Cross-Type Mixing" = colors[["cross_mixing"]],
                                       "Well Integrated" = colors[["well_integrated"]],
                                       "High Reference Proportion" = colors[["reference_only"]]),
                            name = "Community Type"
                        ) +
                        ggplot2::scale_size_continuous(
                            range = c(3, 10),
                            trans = "sqrt",
                            guide = "none"
                        ) +
                        ggplot2::scale_alpha_continuous(
                            range = c(0.3, 1.0),
                            guide = "none"
                        ) +
                        ggplot2::geom_text(
                            data = node_data,
                            ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                         label = .data[["community"]]),
                            color = "black", size = 3, fontface = "bold") +
                        ggplot2::labs(
                            title = "Community Network: Inter-Community Connection Strength",
                            subtitle = paste0("Showing ", nrow(node_data), " communities with ",
                                              nrow(edge_data), " connections"),
                            caption = "Node size proportional to # cells in community; Edge opacity proportional to connection strength"
                        ) +
                        ggplot2::theme_void() +
                        ggplot2::theme(
                            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                            plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
                            legend.position = "right"
                        )

                }
            }
        }

    } else if (plot_type == "cell_network") {
        # Identify reference-only communities (communities with very low query proportion)
        reference_only_comms <- x[["community_composition"]][["community"]][
            x[["community_composition"]][["query_proportion"]] < 0.1
        ]

        # Filter out reference-only communities if requested
        if (exclude_reference_only) {
            keep_cells <- !(x[["cell_info"]][["community"]] %in% reference_only_comms)
            cell_info_filtered <- x[["cell_info"]][keep_cells, ]
            layout_filtered <- x[["graph_info"]][["layout"]][keep_cells, ]
            edges_filtered <- x[["graph_info"]][["edges"]][
                x[["graph_info"]][["edges"]][, 1] %in% which(keep_cells) &
                    x[["graph_info"]][["edges"]][, 2] %in% which(keep_cells),
            ]

            # Remap edge indices
            if (nrow(edges_filtered) > 0) {
                old_to_new <- rep(NA, nrow(x[["cell_info"]]))
                old_to_new[which(keep_cells)] <- seq_len(sum(keep_cells))
                edges_filtered[, 1] <- old_to_new[edges_filtered[, 1]]
                edges_filtered[, 2] <- old_to_new[edges_filtered[, 2]]
            }
        } else {
            cell_info_filtered <- x[["cell_info"]]
            layout_filtered <- x[["graph_info"]][["layout"]]
            edges_filtered <- x[["graph_info"]][["edges"]]
            keep_cells <- rep(TRUE, nrow(x[["cell_info"]]))
        }

        n_nodes <- nrow(layout_filtered)

        if (n_nodes > max_nodes) {
            local_inconsistent_cells <- x[["local_annotation_inconsistencies"]][["query_cell_idx"]]
            # Filter local inconsistent cells to only those remaining after reference_only filtering
            if (exclude_reference_only) {
                local_inconsistent_cells <- local_inconsistent_cells[
                    local_inconsistent_cells %in% which(keep_cells)
                ]
                # Remap indices
                old_to_new <- rep(NA, nrow(x[["cell_info"]]))
                old_to_new[which(keep_cells)] <- seq_len(sum(keep_cells))
                local_inconsistent_cells <- old_to_new[local_inconsistent_cells]
                local_inconsistent_cells <- local_inconsistent_cells[!is.na(local_inconsistent_cells)]
            }

            high_query_prop_comms <- x[["high_query_prop_analysis"]][["community"]]
            cross_mixing_comms <- x[["cross_type_mixing"]][["community"]]

            priority_1 <- which(cell_info_filtered[["community"]] %in%
                                    high_query_prop_comms &
                                    cell_info_filtered[["dataset"]] == "Query")
            priority_2 <- setdiff(local_inconsistent_cells, priority_1)
            priority_3 <- which(cell_info_filtered[["community"]] %in%
                                    cross_mixing_comms &
                                    cell_info_filtered[["dataset"]] == "Query" &
                                    !(seq_len(nrow(cell_info_filtered)) %in%
                                          c(priority_1, priority_2)))
            priority_4 <- which(cell_info_filtered[["community"]] %in%
                                    c(high_query_prop_comms, cross_mixing_comms) &
                                    cell_info_filtered[["dataset"]] == "Reference" &
                                    !(seq_len(nrow(cell_info_filtered)) %in% c(priority_1,
                                                                               priority_2,
                                                                               priority_3)))

            keep_nodes <- c()
            remaining_slots <- max_nodes
            min_slots_per_type <- 50

            if (length(priority_1) > 0) {
                take_p1 <- min(length(priority_1), remaining_slots)
                keep_nodes <- c(keep_nodes, priority_1[1:take_p1])
                remaining_slots <- remaining_slots - take_p1
            }

            if (length(priority_2) > 0 && remaining_slots > 0) {
                take_p2 <- min(length(priority_2), max(min_slots_per_type, remaining_slots))
                keep_nodes <- c(keep_nodes, priority_2[1:take_p2])
                remaining_slots <- remaining_slots - take_p2
            }

            if (length(priority_3) > 0 && remaining_slots > 0) {
                take_p3 <- min(length(priority_3), max(min_slots_per_type, remaining_slots))
                keep_nodes <- c(keep_nodes, priority_3[1:take_p3])
                remaining_slots <- remaining_slots - take_p3
            }

            if (length(priority_4) > 0 && remaining_slots > 0) {
                take_p4 <- min(length(priority_4), remaining_slots * 0.3)
                if (take_p4 > 0) {
                    keep_nodes <- c(keep_nodes, sample(priority_4, take_p4))
                    remaining_slots <- remaining_slots - take_p4
                }
            }

            if (remaining_slots > 0) {
                all_priority <- c(priority_1, priority_2,
                                  priority_3,
                                  priority_4)
                remaining_cells <- setdiff(seq_len(nrow(cell_info_filtered)),
                                           all_priority)
                if (length(remaining_cells) > 0) {
                    take_random <- min(length(remaining_cells), remaining_slots)
                    keep_nodes <- c(keep_nodes, sample(remaining_cells, take_random))
                }
            }

            cell_info_subset <- cell_info_filtered[keep_nodes, ]
            layout_subset <- layout_filtered[keep_nodes, ]

            edges_keep <- edges_filtered[, 1] %in% keep_nodes & edges_filtered[, 2] %in%
                keep_nodes
            edges_subset <- edges_filtered[edges_keep, ]

            if (nrow(edges_subset) > 0) {
                old_to_new <- rep(NA, n_nodes)
                old_to_new[keep_nodes] <- seq_along(keep_nodes)
                edges_subset[, 1] <- old_to_new[edges_subset[, 1]]
                edges_subset[, 2] <- old_to_new[edges_subset[, 2]]

                if (nrow(edges_subset) > 3000) {
                    edges_subset <- edges_subset[sample(nrow(edges_subset), 3000), ]
                }
            }
        } else {
            cell_info_subset <- cell_info_filtered
            layout_subset <- layout_filtered
            edges_subset <- edges_filtered
            keep_nodes <- seq_len(nrow(cell_info_filtered))
            if (nrow(edges_subset) > 3000) {
                edges_subset <- edges_subset[sample(nrow(edges_subset), 3000), ]
            }
        }

        node_data <- data.frame(
            x = layout_subset[, 1],
            y = layout_subset[, 2],
            dataset = cell_info_subset[["dataset"]],
            cell_type = cell_info_subset[["cell_type"]],
            community = cell_info_subset[["community"]],
            stringsAsFactors = FALSE
        )

        # Create edge data
        if (nrow(edges_subset) > 0) {
            edge_data <- data.frame(
                x = layout_subset[edges_subset[, 1], 1],
                y = layout_subset[edges_subset[, 1], 2],
                xend = layout_subset[edges_subset[, 2], 1],
                yend = layout_subset[edges_subset[, 2], 2]
            )
        } else {
            edge_data <- data.frame(x = numeric(0), y = numeric(0),
                                    xend = numeric(0), yend = numeric(0))
        }

        # Create plot based on color_by parameter
        if (color_by == "cell_type") {
            # Generate colors for cell types
            all_cell_types <- sort(unique(node_data[["cell_type"]]))
            cell_types_cases <- c()
            for(cell_type_case in all_cell_types){
                cell_types_cases <- c(cell_types_cases,
                                      paste(c("Reference", "Query"), cell_type_case))
            }
            cell_type_colors <- generateColors(cell_types_cases, paired = TRUE)

            # Set up data for node colors
            node_data[["cell_type_case"]] <- paste(node_data[["dataset"]],
                                                   node_data[["cell_type"]])

            p <- ggplot2::ggplot() +
                {if (nrow(edge_data) > 0) {
                    ggplot2::geom_segment(
                        data = edge_data,
                        ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                     xend = .data[["xend"]], yend = .data[["yend"]]),
                        alpha = 0.2, color = "gray80", linewidth = 0.2)
                }} +
                ggplot2::geom_point(
                    data = node_data,
                    ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                 color = .data[["cell_type_case"]]),
                    size = point_size, alpha = 0.7) +
                ggplot2::scale_color_manual(
                    values = cell_type_colors,
                    name = "Cell Type"
                )

            plot_title <- "Cell Network: Cell Type Distribution"

        } else if (color_by == "community_type") {
            # Add issue type classification
            node_data[["issue_type"]] <- "Well Integrated"

            # Identify reference-only communities for remaining nodes
            reference_only_comms_remaining <- reference_only_comms[
                reference_only_comms %in% unique(node_data[["community"]])
            ]

            node_data[["issue_type"]][node_data[["community"]] %in%
                                          reference_only_comms_remaining] <-
                "High Reference Proportion"
            node_data[["issue_type"]][node_data[["community"]] %in%
                                          x[["cross_type_mixing"]][["community"]]] <-
                "Cross-Type Mixing"
            node_data[["issue_type"]][node_data[["community"]] %in%
                                          x[["high_query_prop_analysis"]][["community"]]] <-
                "High Query Proportion"

            if (nrow(x[["local_annotation_inconsistencies"]]) > 0) {
                local_inconsistent_cells <- x[["local_annotation_inconsistencies"]][["query_cell_idx"]]

                if (exclude_reference_only) {
                    # Map original indices to current subset
                    original_keep_indices <- which(keep_cells)
                    subset_positions <- match(local_inconsistent_cells, original_keep_indices)
                    subset_positions <- subset_positions[!is.na(subset_positions)]

                    if (n_nodes > max_nodes) {
                        # Further map to final keep_nodes
                        final_positions <- match(subset_positions, keep_nodes)
                        valid_positions <- final_positions[!is.na(final_positions)]
                    } else {
                        valid_positions <- subset_positions
                    }
                } else {
                    if (n_nodes > max_nodes) {
                        subset_positions <- match(local_inconsistent_cells, keep_nodes)
                        valid_positions <- subset_positions[!is.na(subset_positions)]
                    } else {
                        valid_positions <- local_inconsistent_cells
                    }
                }

                if (length(valid_positions) > 0) {
                    # Only override if not already High Query Proportion
                    local_mask <- node_data[["issue_type"]][valid_positions] != "High Query Proportion"
                    node_data[["issue_type"]][valid_positions[local_mask]] <- "Local Inconsistent"
                }
            }

            p <- ggplot2::ggplot() +
                {if (nrow(edge_data) > 0) {
                    ggplot2::geom_segment(
                        data = edge_data,
                        ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                     xend = .data[["xend"]], yend = .data[["yend"]]),
                        alpha = 0.2, color = "gray80", linewidth = 0.2)
                }} +
                ggplot2::geom_point(data = node_data,
                                    ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                                 color = .data[["issue_type"]]),
                                    size = point_size, alpha = 0.7) +
                ggplot2::scale_color_manual(
                    values = c("High Query Proportion" = colors[["high_query_prop"]],
                               "Cross-Type Mixing" = colors[["cross_mixing"]],
                               "Local Inconsistent" = colors[["local_inconsistent"]],
                               "Well Integrated" = colors[["well_integrated"]],
                               "High Reference Proportion" = colors[["reference_only"]]),
                    name = "Annotation Status"
                )

            plot_title <- "Cell Network: Annotation Consistency Diagnostics"
        }

        # Create subtitle
        high_query_prop_comms <- nrow(x[["high_query_prop_analysis"]])
        cross_mixing_comms <- nrow(x[["cross_type_mixing"]])
        local_inconsistent_cells <- nrow(x[["local_annotation_inconsistencies"]])

        subtitle_parts <- c()
        if (high_query_prop_comms > 0) {
            subtitle_parts <- c(subtitle_parts,
                                paste0("Query-Only (",
                                       pluralize(high_query_prop_comms,
                                                 "community", "communities"), ")"))
        }
        if (cross_mixing_comms > 0) {
            subtitle_parts <- c(subtitle_parts,
                                paste0("Cross-Type Mixing (",
                                       pluralize(cross_mixing_comms,
                                                 "community", "communities"), ")"))
        }
        if (local_inconsistent_cells > 0) {
            subtitle_parts <- c(subtitle_parts,
                                paste0("Local Inconsistent (",
                                       pluralize(local_inconsistent_cells,
                                                 "cell"), ")"))
        }

        subtitle_text <- if (length(subtitle_parts) > 0) {
            paste(subtitle_parts, collapse = " | ")
        } else {
            "No annotation issues detected"
        }

        # Add titles and theme
        p <- p +
            ggplot2::labs(
                title = plot_title,
                subtitle = subtitle_text,
                caption = if(nrow(x[["graph_info"]][["layout"]]) >
                             max_nodes) paste0(
                                 "Showing ",
                                 nrow(node_data), " of ",
                                 nrow(x[["graph_info"]][["layout"]]),
                                 " cells (prioritizing problematic cells)") else ""
            ) +
            ggplot2::theme_void() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 14,
                                                   face = "bold",
                                                   hjust = 0.5),
                plot.subtitle = ggplot2::element_text(size = 10,
                                                      hjust = 0.5),
                legend.position = "right"
            )

    } else if (plot_type == "community_data") {
        community_composition <- x[["community_composition"]]

        if (nrow(community_composition) == 0) {
            p <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5,
                                  label = "No communities to display",
                                  size = 8, color = "gray50") +
                ggplot2::labs(title = "Community Overview") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 14,
                                                                  face = "bold",
                                                                  hjust = 0.5))
        } else {
            # Ensure n_cell_types is integer
            community_composition[["n_cell_types"]] <-
                as.integer(community_composition[["n_cell_types"]])

            # Filter out reference-only communities if requested
            if (exclude_reference_only) {
                # Identify reference-only communities (communities with very low query proportion)
                reference_only_mask <- community_composition[["query_proportion"]] < 0.1
                community_composition <- community_composition[!reference_only_mask, ]
            }

            if (nrow(community_composition) == 0) {
                p <- ggplot2::ggplot() +
                    ggplot2::annotate("text", x = 0.5, y = 0.5,
                                      label = "No communities to display after filtering",
                                      size = 8, color = "gray50") +
                    ggplot2::labs(title = "Community Overview") +
                    ggplot2::theme_void() +
                    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
            } else {
                # Classify communities
                community_composition[["issue_type"]] <-
                    "Well Integrated"
                community_composition[["issue_type"]][community_composition[["community"]] %in%
                                                          x[["cross_type_mixing"]][["community"]]] <-
                    "Cross-Type Mixing"
                community_composition[["issue_type"]][community_composition[["community"]] %in%
                                                          x[["high_query_prop_analysis"]][["community"]]] <-
                    "High Query Proportion"

                # Add Reference Only communities (communities with very low query proportion)
                community_composition[["issue_type"]][community_composition[["query_proportion"]] < 0.1 &
                                                          community_composition[["issue_type"]] == "Well Integrated"] <-
                    "High Reference Proportion"

                p <- ggplot2::ggplot(community_composition, ggplot2::aes(x = .data[["query_proportion"]],
                                                                         y = .data[["total_cells"]])) +
                    ggplot2::geom_point(ggplot2::aes(color = .data[["issue_type"]],
                                                     size = .data[["n_cell_types"]]),
                                        alpha = 0.7) +
                    ggplot2::scale_color_manual(
                        values = c("High Query Proportion" = colors[["high_query_prop"]],
                                   "Cross-Type Mixing" = colors[["cross_mixing"]],
                                   "Well Integrated" = colors[["well_integrated"]],
                                   "High Reference Proportion" = colors[["reference_only"]]),
                        name = "Community Type"
                    ) +
                    ggplot2::scale_size_continuous(
                        range = c(2, 8), name = "# Cell Types",
                        breaks = function(x) {
                            pretty_breaks <- pretty(x, n = 5)
                            return(pretty_breaks[pretty_breaks == round(pretty_breaks)])
                        },
                        labels = function(x) as.character(as.integer(x))) +
                    ggplot2::geom_text(ggplot2::aes(
                        label = .data[["community"]],
                        vjust = -0.75 - (.data[["n_cell_types"]] - min(.data[["n_cell_types"]])) /
                            (max(.data[["n_cell_types"]]) - min(.data[["n_cell_types"]]) + 1) * 1.5),
                        size = 3) +
                    ggplot2::labs(
                        title = "Community Data Overview",
                        x = "Query Proportion",
                        y = "Total Cells in Community"
                    ) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                        legend.position = "right"
                    )
            }
        }

    } else if (plot_type == "summary") {
        summary_data <- x[["annotation_consistency"]]
        summary_data[["total_issues"]] <- summary_data[["high_query_prop_cells"]] +
            summary_data[["true_cross_mixing_cells"]]

        local_summary <- x[["local_inconsistency_summary"]]
        summary_data <- merge(summary_data,
                              local_summary[, c("cell_type",
                                                "locally_inconsistent_cells")],
                              by = "cell_type")
        summary_data[["total_issues"]] <- summary_data[["total_issues"]] +
            summary_data[["locally_inconsistent_cells"]]
        summary_data[["total_issue_rate"]] <- summary_data[["total_issues"]] /
            summary_data[["total_query_cells"]]

        summary_data[["issue_category"]] <- cut(summary_data[["total_issue_rate"]],
                                                breaks = c(-Inf, 0.05, 0.15, Inf),
                                                labels = c("Excellent",
                                                           "Concerning",
                                                           "Problematic"))

        p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = reorder(.data[["cell_type"]],
                                                                    .data[["total_issue_rate"]]))) +
            ggplot2::geom_col(ggplot2::aes(y = .data[["total_issue_rate"]],
                                           fill = .data[["issue_category"]]),
                              alpha = 0.8, color = "black", linewidth = 0.3) +
            ggplot2::scale_fill_manual(
                values = c("Excellent" = colors[["excellent"]],
                           "Concerning" = colors[["concerning"]],
                           "Problematic" = colors[["problematic"]]),
                name = "Issue Level"
            ) +
            ggplot2::geom_text(ggplot2::aes(y = .data[["total_issue_rate"]],
                                            label = paste0(.data[["total_issues"]], "/",
                                                           .data[["total_query_cells"]])),
                               vjust = -0.5, size = 3) +
            ggplot2::labs(
                title = "Total Annotation Issues by Cell Type",
                subtitle = "Includes query-only, cross-mixing, and local inconsistencies",
                x = "Cell Type",
                y = "Proportion of Query Cells with Issues",
                caption = paste0("Mean issue rate: ",
                                 round(mean(summary_data[["total_issue_rate"]]), 3))
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                legend.position = "bottom"
            ) +
            ggplot2::coord_cartesian(ylim = c(0, max(1, max(summary_data[["total_issue_rate"]]) * 1.1)))

    } else if (plot_type == "local_issues") {
        if (nrow(x[["local_annotation_inconsistencies"]]) == 0) {
            p <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.6,
                                  label = "No local annotation inconsistencies detected!",
                                  size = 8, color = colors[["excellent"]], fontface = "bold") +
                ggplot2::annotate("text", x = 0.5, y = 0.4,
                                  label = "All query annotations supported by local neighborhoods",
                                  size = 6, color = colors[["excellent"]]) +
                ggplot2::labs(title = "Local Annotation Consistency") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
        } else {
            local_data <- x[["local_inconsistency_summary"]]
            local_data <- local_data[local_data[["locally_inconsistent_cells"]] > 0, ]

            if (nrow(local_data) > 0) {
                local_data[["inconsistency_category"]] <- cut(
                    local_data[["local_inconsistency_rate"]],
                    breaks = c(-Inf, 0.1, 0.3, Inf),
                    labels = c("Minor", "Moderate", "Severe"))

                p <- ggplot2::ggplot(local_data, ggplot2::aes(x = reorder(
                    .data[["cell_type"]],
                    .data[["local_inconsistency_rate"]]))) +
                    ggplot2::geom_col(ggplot2::aes(
                        y = .data[["local_inconsistency_rate"]],
                        fill = .data[["inconsistency_category"]]),
                        alpha = 0.8, color = "black", linewidth = 0.3) +
                    ggplot2::scale_fill_manual(
                        values = c("Minor" = colors[["excellent"]],
                                   "Moderate" = colors[["concerning"]],
                                   "Severe" = colors[["problematic"]]),
                        name = "Severity"
                    ) +
                    ggplot2::geom_text(ggplot2::aes(
                        y = .data[["local_inconsistency_rate"]],
                        label = paste0(.data[["locally_inconsistent_cells"]], "/",
                                       .data[["total_query_cells"]])),
                        vjust = -0.5, size = 3) +
                    ggplot2::labs(
                        title = "Local Annotation Inconsistencies by Cell Type",
                        subtitle = "Query cells not supported by reference neighbors",
                        x = "Cell Type",
                        y = "Local Inconsistency Rate",
                        caption = paste0(
                            "Mean rate: ",
                            round(x[["overall_metrics"]][["mean_local_inconsistency_rate"]], 3))
                    ) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                        legend.position = "bottom"
                    )
            } else {
                p <- ggplot2::ggplot() +
                    ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data to display", size = 6) +
                    ggplot2::theme_void()
            }
        }

    } else if (plot_type == "annotation_issues") {
        total_issues <- nrow(x[["high_query_prop_analysis"]]) + nrow(x[["cross_type_mixing"]]) +
            nrow(x[["local_annotation_inconsistencies"]])

        if (total_issues == 0) {
            p <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.6,
                                  label = "No annotation issues detected!",
                                  size = 8, color = colors[["excellent"]], fontface = "bold") +
                ggplot2::annotate("text", x = 0.5, y = 0.4,
                                  label = "Excellent annotation consistency across all metrics",
                                  size = 6, color = colors[["excellent"]]) +
                ggplot2::labs(title = "Comprehensive Annotation Assessment") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
        } else {
            high_query_prop_comms <- nrow(x[["high_query_prop_analysis"]])
            cross_mixing_comms <- nrow(x[["cross_type_mixing"]])
            local_inconsistent_cells <- nrow(x[["local_annotation_inconsistencies"]])

            issue_summary <- data.frame(
                issue_type = character(0),
                count = numeric(0),
                affected_cells = numeric(0),
                label_text = character(0),
                stringsAsFactors = FALSE
            )

            if (high_query_prop_comms > 0) {
                issue_summary <- rbind(issue_summary, data.frame(
                    issue_type = "Query-Only Communities",
                    count = high_query_prop_comms,
                    affected_cells = sum(x[["high_query_prop_analysis"]][["total_query_cells"]]),
                    label_text = paste0(
                        "Query-Only Communities\n(", pluralize(high_query_prop_comms,
                                                               "community", "communities"), ")"),
                    stringsAsFactors = FALSE
                ))
            }

            if (cross_mixing_comms > 0) {
                issue_summary <- rbind(issue_summary, data.frame(
                    issue_type = "Cross-Type Mixing",
                    count = cross_mixing_comms,
                    affected_cells = sum(x[["cross_type_mixing"]][["n_query"]]),
                    label_text = paste0(
                        "Cross-Type Mixing\n(", pluralize(cross_mixing_comms,
                                                          "community", "communities"), ")"),
                    stringsAsFactors = FALSE
                ))
            }

            if (local_inconsistent_cells > 0) {
                issue_summary <- rbind(issue_summary, data.frame(
                    issue_type = "Local Inconsistencies",
                    count = local_inconsistent_cells,
                    affected_cells = local_inconsistent_cells,
                    label_text = paste0(
                        "Local Inconsistencies\n(",
                        pluralize(local_inconsistent_cells, "cell"), ")"),
                    stringsAsFactors = FALSE
                ))
            }

            p <- ggplot2::ggplot(issue_summary, ggplot2::aes(
                x = reorder(.data[["label_text"]], .data[["affected_cells"]]),
                y = .data[["affected_cells"]])) +
                ggplot2::geom_col(ggplot2::aes(fill = .data[["issue_type"]]),
                                  alpha = 0.8, color = "black", linewidth = 0.3) +
                ggplot2::scale_fill_manual(
                    values = c("Query-Only Communities" = colors[["high_query_prop"]],
                               "Cross-Type Mixing" = colors[["cross_mixing"]],
                               "Local Inconsistencies" = colors[["local_inconsistent"]]),
                    name = "Issue Type"
                ) +
                ggplot2::labs(
                    title = "Annotation Issues Summary",
                    subtitle = "Number of affected query cells by issue type",
                    x = "Issue Type",
                    y = "Number of Affected Query Cells"
                ) +
                ggplot2::theme_bw() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 10),
                    axis.text.y = ggplot2::element_text(size = 10),
                    plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                    legend.position = "none",
                    plot.margin = ggplot2::margin(20, 40, 20, 40)
                ) +
                ggplot2::coord_flip() +
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))
        }
    }

    return(p)
}
