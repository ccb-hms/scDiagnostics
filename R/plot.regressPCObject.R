#' @title Plot Regression Results on Principal Components
#'
#' @description
#' The S3 plot method generates plots to visualize the results of regression analyses
#' performed on principal components (PCs) against cell types, datasets, or their interactions.
#'
#' @param x An object of class \code{regressPCObject} containing the output of the \code{regressPC} function
#' @param plot_type Type of plot to generate. Available options depend on analysis type.
#' @param alpha Significance threshold for p-values. Default is 0.05.
#' @param ... Additional arguments to be passed to the plotting functions.
#'
#' @return The S3 plot method returns a \code{ggplot} object representing the specified plot type.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{regressPC}}
#'
#' @rdname regressPC
#'
# Plot the results of regressPC
plot.regressPCObject <- function(x,
                                 plot_type = NULL,
                                 alpha = 0.05,
                                 ...){

    # Helper function to get variance explained percentages
    .getVarianceLabels <- function(pc_names, var_explained = NULL) {
        if(is.null(var_explained)) {
            return(pc_names)
        } else {
            pc_nums <- as.numeric(gsub("PC", "", pc_names))
            percentages <- round(var_explained[pc_nums], 1)
            return(paste0(pc_names, " (", percentages, "%)"))
        }
    }

    # Helper function to clean coefficient names
    .cleanCoeffNames <- function(coeff_names) {
        cleaned_names <- coeff_names

        # Process interaction terms (contain ":")
        interaction_idx <- grepl(":", coeff_names)
        if(any(interaction_idx)) {
            interaction_terms <- coeff_names[interaction_idx]

            for(i in seq_along(interaction_terms)) {
                term <- interaction_terms[i]
                parts <- strsplit(term, ":")[[1]]
                cleaned_parts <- character(length(parts))
                for(j in seq_along(parts)) {
                    part <- parts[j]
                    part <- gsub("^cell_type", "", part)
                    part <- gsub("^dataset", "", part)
                    part <- gsub("^batch", "", part)
                    cleaned_parts[j] <- part
                }

                if(length(cleaned_parts) == 2) {
                    if(grepl("Query|Reference", cleaned_parts[1])) {
                        cleaned_names[interaction_idx][i] <-
                            paste(cleaned_parts[1], cleaned_parts[2], sep = ":")
                    } else if(grepl("Query|Reference", cleaned_parts[2])) {
                        cleaned_names[interaction_idx][i] <-
                            paste(cleaned_parts[2], cleaned_parts[1], sep = ":")
                    } else {
                        cleaned_names[interaction_idx][i] <-
                            paste(cleaned_parts, collapse = ":")
                    }
                } else {
                    cleaned_names[interaction_idx][i] <-
                        paste(cleaned_parts, collapse = ":")
                }
            }
        }

        # Process non-interaction terms
        non_interaction_idx <- !interaction_idx
        if(any(non_interaction_idx)) {
            non_interaction_terms <- coeff_names[non_interaction_idx]
            cleaned_terms <- gsub("^cell_type", "", non_interaction_terms)
            cleaned_terms <- gsub("^dataset", "", cleaned_terms)
            cleaned_terms <- gsub("^batch", "", cleaned_terms)
            cleaned_names[non_interaction_idx] <- cleaned_terms
        }

        return(cleaned_names)
    }

    # Get variance explained - use reference if available, otherwise query
    var_explained <- if("reference_pca_var" %in% names(x)) {
        x[["reference_pca_var"]]
    } else if("query_pca_var" %in% names(x)) {
        x[["query_pca_var"]]
    } else {
        NULL
    }

    # =========================================================================
    # Case 1: Query only - PC ~ cell_type
    # =========================================================================
    if(x[["indep_var"]] == "cell_type"){

        if(is.null(plot_type)) plot_type <- "r_squared"
        valid_types <- c("r_squared", "heatmap", "variance_contribution")
        if(!plot_type %in% valid_types) {
            stop("For cell_type analysis, plot_type must be one of: ",
                 paste(valid_types, collapse = ", "))
        }

        if(plot_type == "r_squared"){
            pc_names <- names(x[["regression_summaries"]])
            pc_labels <- .getVarianceLabels(pc_names, var_explained)

            plot_data <- data.frame(
                PC = factor(pc_labels, levels = pc_labels),
                r_squared = x[["r_squared"]]
            )

            plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]],
                                                                y = .data[["r_squared"]])) +
                ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
                ggplot2::labs(title = expression(R^2 ~ "of PC ~ Cell Types"),
                              subtitle = "How well do cell types separate on each PC?",
                              x = "", y = expression(R^2)) +
                ggplot2::theme_bw()

        } else if(plot_type == "heatmap"){
            coeff_data <- data.frame()

            for(pc_name in names(x[["regression_summaries"]])){
                coeffs <- x[["regression_summaries"]][[pc_name]][["coefficients"]]

                # Remove intercept row
                coeff_rows <- rownames(coeffs)
                keep_rows <- !grepl("^\\(Intercept\\)$", coeff_rows)

                if(sum(keep_rows) > 0){
                    filtered_coeffs <- coeffs[keep_rows, , drop = FALSE]
                    clean_names <- .cleanCoeffNames(rownames(filtered_coeffs))

                    pc_data <- data.frame(
                        PC = pc_name,
                        cell_type = clean_names,
                        coefficient = filtered_coeffs[, "coef"],
                        p_value = filtered_coeffs[, "p.adjusted"],
                        stringsAsFactors = FALSE
                    )
                    coeff_data <- rbind(coeff_data, pc_data)
                }
            }

            if(nrow(coeff_data) == 0) stop("No coefficients found")

            # Fix PC ordering and add variance labels
            pc_nums <- as.numeric(gsub("PC", "", coeff_data[["PC"]]))
            unique_pc_names <- paste0("PC", sort(unique(pc_nums)))
            pc_labels <- .getVarianceLabels(unique_pc_names, var_explained)

            # Map PC names to labels
            pc_mapping <- setNames(pc_labels, unique_pc_names)
            coeff_data[["PC_labeled"]] <- pc_mapping[coeff_data[["PC"]]]
            coeff_data[["PC_labeled"]] <- factor(coeff_data[["PC_labeled"]], levels = pc_labels)
            coeff_data[["cell_type"]] <- factor(coeff_data[["cell_type"]],
                                                levels = unique(coeff_data[["cell_type"]]))
            coeff_data[["significant"]] <- coeff_data[["p_value"]] < alpha

            plot_obj <- ggplot2::ggplot(coeff_data, ggplot2::aes(x = .data[["PC_labeled"]],
                                                                 y = .data[["cell_type"]],
                                                                 fill = .data[["coefficient"]])) +
                ggplot2::geom_tile(color = "white") +
                ggplot2::geom_point(data = coeff_data[coeff_data[["significant"]], ],
                                    shape = 8, size = 2, color = "black") +
                ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                              name = "Coefficient") +
                ggplot2::scale_x_discrete(drop = TRUE) +  # This should fix extra x-axis ticks
                ggplot2::scale_y_discrete(drop = TRUE) +  # This should fix extra y-axis ticks
                ggplot2::labs(title = "Cell Type Coefficients by PC",
                              subtitle = "Asterisks indicate significant coefficients",
                              x = "", y = "Cell Type") +
                ggplot2::theme_bw()

        } else if(plot_type == "variance_contribution"){
            if(!"var_contributions" %in% names(x)) stop("Variance contributions not available")

            pc_names <- names(x[["regression_summaries"]])
            pc_labels <- .getVarianceLabels(pc_names, var_explained)

            plot_data <- data.frame(
                PC = factor(pc_labels, levels = pc_labels),
                var_contribution = x[["var_contributions"]]
            )

            plot_obj <- ggplot2::ggplot(plot_data,
                                        ggplot2::aes(x = .data[["PC"]],
                                                     y = .data[["var_contribution"]])) +
                ggplot2::geom_col(fill = "darkgreen", alpha = 0.7) +
                ggplot2::labs(title = "Variance Contribution by PC",
                              x = "", y = "Variance Contribution") +
                ggplot2::theme_bw()
        }

        # =========================================================================
        # Case 2: Query + Reference - PC ~ cell_type * dataset (unified model)
        # =========================================================================
    } else if(x[["indep_var"]] == "cell_type_dataset_interaction"){

        if(is.null(plot_type)) plot_type <- "r_squared"
        valid_types <- c("r_squared", "heatmap", "interaction_effects")
        if(!plot_type %in% valid_types) {
            stop("For cell_type_dataset_interaction analysis, plot_type must be one of: ",
                 paste(valid_types, collapse = ", "))
        }

        if(plot_type == "r_squared"){
            pc_names <- names(x[["regression_summaries"]])
            pc_labels <- .getVarianceLabels(pc_names, var_explained)

            plot_data <- data.frame(
                PC = factor(pc_labels, levels = pc_labels),
                r_squared = x[["r_squared"]]
            )

            plot_obj <- ggplot2::ggplot(plot_data,
                                        ggplot2::aes(x = .data[["PC"]],
                                                     y = .data[["r_squared"]])) +
                ggplot2::geom_col(fill = "orange", alpha = 0.7) +
                ggplot2::labs(title = expression(R^2 ~ "of PC ~ Cell Type x Dataset Model"),
                              subtitle = "How much variance is explained by cell types and dataset differences?",
                              x = "", y = expression(R^2)) +
                ggplot2::theme_bw()

        } else if(plot_type == "heatmap"){
            coeff_data <- data.frame()

            for(pc_name in names(x[["regression_summaries"]])){
                coeffs <- x[["regression_summaries"]][[pc_name]][["coefficients"]]

                coeff_rows <- rownames(coeffs)
                keep_rows <- !grepl("^\\(Intercept\\)$", coeff_rows)

                if(sum(keep_rows) > 0){
                    filtered_coeffs <- coeffs[keep_rows, , drop = FALSE]
                    coeff_names <- rownames(filtered_coeffs)
                    clean_names <- .cleanCoeffNames(coeff_names)

                    # Categorize coefficient types
                    coeff_types <- ifelse(grepl(":", coeff_names), "Interaction",
                                          ifelse(grepl("dataset", coeff_names),
                                                 "Dataset", "Cell Type"))

                    pc_data <- data.frame(
                        PC = pc_name,
                        coefficient_name = clean_names,
                        coefficient_type = coeff_types,
                        coefficient = filtered_coeffs[, "coef"],
                        p_value = filtered_coeffs[, "p.adjusted"],
                        stringsAsFactors = FALSE
                    )
                    coeff_data <- rbind(coeff_data, pc_data)
                }
            }

            if(nrow(coeff_data) == 0) stop("No coefficients found")

            # Fix PC ordering and add variance labels
            pc_nums <- as.numeric(gsub("PC", "", coeff_data[["PC"]]))
            unique_pc_names <- paste0("PC", sort(unique(pc_nums)))
            pc_labels <- .getVarianceLabels(unique_pc_names, var_explained)

            # Map PC names to labels
            pc_mapping <- setNames(pc_labels, unique_pc_names)
            coeff_data[["PC_labeled"]] <- pc_mapping[coeff_data[["PC"]]]
            coeff_data[["PC_labeled"]] <- factor(coeff_data[["PC_labeled"]], levels = pc_labels)
            coeff_data[["coefficient_name"]] <- factor(coeff_data[["coefficient_name"]],
                                                       levels = unique(coeff_data[["coefficient_name"]]))

            # Clean coefficient_type and ensure no NAs
            coeff_data <- coeff_data[!is.na(coeff_data[["coefficient_type"]]) &
                                         coeff_data[["coefficient_type"]] != "", ]
            coeff_data[["coefficient_type"]] <- factor(coeff_data[["coefficient_type"]],
                                                       levels = unique(coeff_data[["coefficient_type"]]))
            coeff_data[["significant"]] <- coeff_data[["p_value"]] < alpha

            plot_obj <- ggplot2::ggplot(coeff_data,
                                        ggplot2::aes(x = .data[["PC_labeled"]],
                                                     y = .data[["coefficient_name"]],
                                                     fill = .data[["coefficient"]])) +
                ggplot2::geom_tile(color = "white") +
                ggplot2::geom_point(data = coeff_data[coeff_data[["significant"]], ],
                                    shape = 8, size = 2, color = "black") +
                ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                              name = "Coefficient") +
                ggplot2::scale_x_discrete(drop = TRUE) +  # Fix extra x-axis ticks
                ggplot2::scale_y_discrete(drop = TRUE) +  # Fix extra y-axis ticks
                ggplot2::facet_grid(coefficient_type ~ ., scales = "free_y", space = "free_y") +
                ggplot2::labs(title = "Cell Type x Dataset Interaction Model Coefficients",
                              subtitle = "Asterisks indicate significant coefficients",
                              x = "", y = "Model Terms") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

        } else if(plot_type == "interaction_effects"){
            interaction_data <- data.frame()

            for(pc_name in names(x[["regression_summaries"]])){
                coeffs <- x[["regression_summaries"]][[pc_name]][["coefficients"]]

                coeff_rows <- rownames(coeffs)
                interaction_rows <- grepl(":", coeff_rows)

                if(sum(interaction_rows) > 0){
                    interaction_coeffs <- coeffs[interaction_rows, , drop = FALSE]
                    interaction_names <- rownames(interaction_coeffs)
                    clean_names <- .cleanCoeffNames(interaction_names)

                    pc_data <- data.frame(
                        PC = pc_name,
                        interaction = clean_names,
                        coefficient = interaction_coeffs[, "coef"],
                        p_value = interaction_coeffs[, "p.adjusted"],
                        stringsAsFactors = FALSE
                    )
                    interaction_data <- rbind(interaction_data, pc_data)
                }
            }

            if(nrow(interaction_data) == 0) stop("No interaction terms found")

            # Fix PC ordering and add variance labels
            pc_nums <- as.numeric(gsub("PC", "", interaction_data[["PC"]]))
            unique_pc_names <- paste0("PC", sort(unique(pc_nums)))
            pc_labels <- .getVarianceLabels(unique_pc_names, var_explained)

            # Map PC names to labels
            pc_mapping <- setNames(pc_labels, unique_pc_names)
            interaction_data[["PC_labeled"]] <- pc_mapping[interaction_data[["PC"]]]
            interaction_data[["PC_labeled"]] <- factor(interaction_data[["PC_labeled"]], levels = pc_labels)
            interaction_data[["interaction"]] <- factor(interaction_data[["interaction"]],
                                                        levels = unique(interaction_data[["interaction"]]))
            interaction_data[["significant"]] <- interaction_data[["p_value"]] < alpha

            plot_obj <- ggplot2::ggplot(interaction_data,
                                        ggplot2::aes(x = .data[["PC_labeled"]],
                                                     y = .data[["interaction"]],
                                                     fill = .data[["coefficient"]])) +
                ggplot2::geom_tile(color = "white") +
                ggplot2::geom_point(data = interaction_data[interaction_data[["significant"]], ],
                                    shape = 8, size = 2, color = "black") +
                ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                              name = "Interaction\nCoefficient") +
                ggplot2::scale_x_discrete(drop = TRUE) +  # Fix extra x-axis ticks
                ggplot2::scale_y_discrete(drop = TRUE) +  # Fix extra y-axis ticks
                ggplot2::labs(title = "Cell Type x Dataset Interaction Effects",
                              subtitle = "Large interactions suggest annotation problems for specific cell types",
                              x = "", y = "Interaction Term") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
        }

        # =========================================================================
        # Case 3: Batch interaction models - PC ~ cell_type * batch
        # =========================================================================
    } else if(x[["indep_var"]] == "cell_type_batch_interaction"){

        if(is.null(plot_type)) plot_type <- "r_squared"
        valid_types <- c("r_squared", "heatmap", "interaction_effects")
        if(!plot_type %in% valid_types) {
            stop("For cell_type_batch_interaction analysis, plot_type must be one of: ",
                 paste(valid_types, collapse = ", "))
        }

        if(plot_type == "r_squared"){
            pc_names <- names(x[["regression_summaries"]])
            pc_labels <- .getVarianceLabels(pc_names, var_explained)

            plot_data <- data.frame(
                PC = factor(pc_labels, levels = pc_labels),
                r_squared = x[["r_squared"]]
            )

            plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]],
                                                                y = .data[["r_squared"]])) +
                ggplot2::geom_col(fill = "purple", alpha = 0.7) +
                ggplot2::labs(title = expression(R^2 ~ "of PC ~ Cell Type x Batch Model"),
                              subtitle = "How much variance is explained by cell types and batch effects?",
                              x = "", y = expression(R^2)) +
                ggplot2::theme_bw()

        } else if(plot_type == "heatmap"){
            coeff_data <- data.frame()

            for(pc_name in names(x[["regression_summaries"]])){
                coeffs <- x[["regression_summaries"]][[pc_name]][["coefficients"]]

                coeff_rows <- rownames(coeffs)
                keep_rows <- !grepl("^\\(Intercept\\)$", coeff_rows)

                if(sum(keep_rows) > 0){
                    filtered_coeffs <- coeffs[keep_rows, , drop = FALSE]
                    coeff_names <- rownames(filtered_coeffs)
                    clean_names <- .cleanCoeffNames(coeff_names)

                    coeff_types <- ifelse(grepl(":", coeff_names), "Interaction",
                                          ifelse(grepl("batch", coeff_names), "Batch", "Cell Type"))

                    pc_data <- data.frame(
                        PC = pc_name,
                        coefficient_name = clean_names,
                        coefficient_type = coeff_types,
                        coefficient = filtered_coeffs[, "coef"],
                        p_value = filtered_coeffs[, "p.adjusted"],
                        stringsAsFactors = FALSE
                    )
                    coeff_data <- rbind(coeff_data, pc_data)
                }
            }

            if(nrow(coeff_data) == 0) stop("No coefficients found")

            # Fix PC ordering and add variance labels
            pc_nums <- as.numeric(gsub("PC", "", coeff_data[["PC"]]))
            unique_pc_names <- paste0("PC", sort(unique(pc_nums)))
            pc_labels <- .getVarianceLabels(unique_pc_names, var_explained)

            # Map PC names to labels
            pc_mapping <- setNames(pc_labels, unique_pc_names)
            coeff_data[["PC_labeled"]] <- pc_mapping[coeff_data[["PC"]]]
            coeff_data[["PC_labeled"]] <- factor(coeff_data[["PC_labeled"]], levels = pc_labels)
            coeff_data[["coefficient_name"]] <- factor(coeff_data[["coefficient_name"]],
                                                       levels = unique(coeff_data[["coefficient_name"]]))

            # Clean coefficient_type and ensure no NAs
            coeff_data <- coeff_data[!is.na(coeff_data[["coefficient_type"]]) &
                                         coeff_data[["coefficient_type"]] != "", ]
            coeff_data[["coefficient_type"]] <- factor(coeff_data[["coefficient_type"]],
                                                       levels = unique(coeff_data[["coefficient_type"]]))
            coeff_data[["significant"]] <- coeff_data[["p_value"]] < alpha

            plot_obj <- ggplot2::ggplot(coeff_data, ggplot2::aes(x = .data[["PC_labeled"]],
                                                                 y = .data[["coefficient_name"]],
                                                                 fill = .data[["coefficient"]])) +
                ggplot2::geom_tile(color = "white") +
                ggplot2::geom_point(data = coeff_data[coeff_data[["significant"]], ],
                                    shape = 8, size = 2, color = "black") +
                ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                              name = "Coefficient") +
                ggplot2::scale_x_discrete(drop = TRUE) +  # Fix extra x-axis ticks
                ggplot2::scale_y_discrete(drop = TRUE) +  # Fix extra y-axis ticks
                ggplot2::facet_grid(coefficient_type ~ ., scales = "free_y", space = "free_y") +
                ggplot2::labs(title = "Cell Type x Batch Interaction Model Coefficients",
                              subtitle = "Asterisks indicate significant coefficients",
                              x = "", y = "Model Terms") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

        } else if(plot_type == "interaction_effects"){
            interaction_data <- data.frame()

            for(pc_name in names(x[["regression_summaries"]])){
                coeffs <- x[["regression_summaries"]][[pc_name]][["coefficients"]]

                coeff_rows <- rownames(coeffs)
                interaction_rows <- grepl(":", coeff_rows)

                if(sum(interaction_rows) > 0){
                    interaction_coeffs <- coeffs[interaction_rows, , drop = FALSE]
                    interaction_names <- rownames(interaction_coeffs)
                    clean_names <- .cleanCoeffNames(interaction_names)

                    pc_data <- data.frame(
                        PC = pc_name,
                        interaction = clean_names,
                        coefficient = interaction_coeffs[, "coef"],
                        p_value = interaction_coeffs[, "p.adjusted"],
                        stringsAsFactors = FALSE
                    )
                    interaction_data <- rbind(interaction_data, pc_data)
                }
            }

            if(nrow(interaction_data) == 0) stop("No interaction terms found")

            # Fix PC ordering and add variance labels
            pc_nums <- as.numeric(gsub("PC", "", interaction_data[["PC"]]))
            unique_pc_names <- paste0("PC", sort(unique(pc_nums)))
            pc_labels <- .getVarianceLabels(unique_pc_names, var_explained)

            # Map PC names to labels
            pc_mapping <- setNames(pc_labels, unique_pc_names)
            interaction_data[["PC_labeled"]] <- pc_mapping[interaction_data[["PC"]]]
            interaction_data[["PC_labeled"]] <- factor(interaction_data[["PC_labeled"]], levels = pc_labels)
            interaction_data[["interaction"]] <- factor(interaction_data[["interaction"]],
                                                        levels = unique(interaction_data[["interaction"]]))
            interaction_data[["significant"]] <- interaction_data[["p_value"]] < alpha

            plot_obj <- ggplot2::ggplot(interaction_data, ggplot2::aes(x = .data[["PC_labeled"]],
                                                                       y = .data[["interaction"]],
                                                                       fill = .data[["coefficient"]])) +
                ggplot2::geom_tile(color = "white") +
                ggplot2::geom_point(data = interaction_data[interaction_data[["significant"]], ],
                                    shape = 8, size = 2, color = "black") +
                ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                              name = "Interaction\nCoefficient") +
                ggplot2::scale_x_discrete(drop = TRUE) +  # Fix extra x-axis ticks
                ggplot2::scale_y_discrete(drop = TRUE) +  # Fix extra y-axis ticks
                ggplot2::labs(title = "Cell Type x Batch Interaction Effects",
                              subtitle = "Large interactions suggest batch-specific cell type problems",
                              x = "", y = "Interaction Term") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
        }

    } else {
        stop("Unsupported analysis type: ", x[["indep_var"]])
    }

    return(plot_obj)
}
