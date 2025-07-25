#' @title Plot Regression Results on Principal Components
#'
#' @description
#' The S3 plot method generates plots to visualize the results of regression analyses
#' performed on principal components (PCs) against cell types, datasets, or their interactions.
#'
#' @param x An object of class \code{regressPCObject} containing the output of the \code{regressPC} function
#' @param plot_type Type of plot to generate. Available options:
#'   "r_squared", "variance_contribution", "coefficient_heatmap"
#' @param alpha Significance threshold for p-values. Default is 0.05.
#' @param coefficients_include Character vector specifying which coefficient types to include
#'   in the coefficient heatmap. Options are \code{c("cell_type", "batch", "interaction")}.
#'   Default is \code{NULL}, which includes all available coefficient types. Only applies
#'   to \code{plot_type = "coefficient_heatmap"}.
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
#' @importFrom stats aggregate formula complete.cases
#'
# Plot the results of regressPC
plot.regressPCObject <- function(x,
                                 plot_type = c("r_squared", "variance_contribution",
                                               "coefficient_heatmap"),
                                 alpha = 0.05,
                                 coefficients_include = NULL,
                                 ...) {

    # Match plot_type argument
    plot_type <- match.arg(plot_type)

    # Validate coefficients_include parameter
    if (!is.null(coefficients_include)) {
        # Check if coefficients_include is valid
        valid_coefficients <- c("cell_type", "batch", "interaction")
        if (!all(coefficients_include %in% valid_coefficients)) {
            stop("coefficients_include must be one or more of: ",
                 paste(valid_coefficients, collapse = ", "))
        }

        # Check if at least one coefficient type is specified
        if (length(coefficients_include) == 0) {
            stop("coefficients_include must contain at least one coefficient type")
        }

        # Check if the specified coefficients are available in the object
        available_coefficients <- .getAvailableCoefficients(x)
        unavailable <- coefficients_include[!coefficients_include %in% available_coefficients]
        if (length(unavailable) > 0) {
            stop("The following coefficient types are not available in the regression object: ",
                 paste(unavailable, collapse = ", "),
                 ". Available types: ", paste(available_coefficients, collapse = ", "))
        }
    }

    # Generate plot based on type
    switch(plot_type,
           "r_squared" = plotRSquared(x, ...),
           "variance_contribution" = plotVarianceContribution(x, ...),
           "coefficient_heatmap" = plotCoefficientHeatmap(x, alpha, coefficients_include, ...)
    )
}

#' @title Determine Available Coefficient Types
#'
#' @description
#' Helper function to determine which coefficient types are available in a regressPCObject
#' based on the model type and data structure.
#'
#' @param x An object of class \code{regressPCObject} containing regression results
#'   from \code{regressPC} function.
#'
#' @keywords internal
#'
#' @return A character vector of available coefficient types.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
.getAvailableCoefficients <- function(x) {
    available <- c("cell_type")  # cell_type is always available

    # Check for batch coefficients
    if (x[["indep_var"]] %in% c("cell_type_batch_interaction", "cell_type_dataset_interaction")) {
        available <- c(available, "batch")
    }

    # Check for interaction coefficients
    if (x[["indep_var"]] %in% c("cell_type_batch_interaction", "cell_type_dataset_interaction")) {
        available <- c(available, "interaction")
    }

    return(available)
}

# [plotRSquared and plotVarianceContribution functions remain unchanged]
#' @title Generate R-squared Bar Plot with Component Breakdown
#'
#' @description
#' Creates a bar plot visualization of R-squared values for each principal component,
#' with optional stacked bars showing the contribution of individual model components
#' (cell type, batch/dataset, and interaction effects) when component decomposition
#' is available.
#'
#' @details
#' This function generates either a simple bar plot or a stacked bar plot depending
#' on the availability of component-wise R-squared decomposition in the input object.
#' When component breakdown is available, the bars are stacked to show:
#' \itemize{
#'   \item Cell type main effect (blue)
#'   \item Batch or dataset main effect (orange)
#'   \item Interaction effect (green)
#' }
#'
#' Principal component labels include the percentage of total variance explained by
#' each PC. The total R-squared value is displayed above each bar. For query-only
#' analyses, the function uses query PCA variance; for query+reference analyses,
#' it uses reference PCA variance.
#'
#' @param x An object of class \code{regressPCObject} containing regression results
#'   from \code{regressPC} function.
#' @param ... Additional arguments passed to the plotting function (currently unused).
#'
#' @keywords internal
#'
#' @return A \code{ggplot2} object representing the R-squared bar plot with optional
#'   component breakdown.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_fill_manual labs theme_minimal theme element_text element_blank element_line margin ylim
#' @importFrom tools toTitleCase
#'
# Helper function: R-squared barplot with component breakdown
plotRSquared <- function(x, ...) {

    # Helper function to get model text
    .getModelText <- function(indep_var) {
        switch(indep_var,
               "cell_type" = "PC ~ Cell Type",
               "cell_type_batch_interaction" = "PC ~ Cell Type * Batch",
               "cell_type_dataset_interaction" = "PC ~ Cell Type * Dataset",
               "Unknown Model")
    }

    # Determine which PCA variance to use based on data type
    if (!is.null(x[["reference_pca_var"]])) {
        # Reference + Query data
        pca_var <- x[["reference_pca_var"]]
    } else {
        # Query only data
        pca_var <- x[["query_pca_var"]]
    }

    # Get PC indices for variance lookup
    pc_indices <- as.numeric(gsub("PC", "", names(x[["r_squared"]])))

    # Create PC names with variance labels
    pc_names <- names(x[["r_squared"]])
    pc_labels <- paste0(pc_names, "\n(", sprintf("%.1f%%", pca_var[pc_indices]), ")")

    # Check if component breakdown is available
    if (!is.null(x[["r_squared_components"]]) && length(x[["r_squared_components"]]) > 1) {
        # Create stacked bar plot with components

        # Prepare component data
        component_data <- list()
        for (component_name in names(x[["r_squared_components"]])) {
            component_data[[component_name]] <- data.frame(
                PC = factor(pc_labels, levels = pc_labels),
                Component = component_name,
                R_squared = x[["r_squared_components"]][[component_name]],
                stringsAsFactors = FALSE
            )
        }

        # Combine all components
        plot_data <- do.call(rbind, component_data)

        # Create proper component labels and colors
        plot_data[["Component"]] <- factor(plot_data[["Component"]],
                                           levels = names(x[["r_squared_components"]]))

        # Set up colors based on available components
        component_colors <- c(
            "cell_type" = "#1f77b4",     # Blue
            "batch" = "#ff7f0e",         # Orange
            "dataset" = "#ff7f0e",       # Orange (same as batch)
            "interaction" = "#2ca02c"    # Green
        )

        # Get colors for available components
        available_colors <- component_colors[names(x[["r_squared_components"]])]

        # Create component labels for legend
        component_labels <- names(x[["r_squared_components"]])
        component_labels <- gsub("_", " ", component_labels)
        component_labels <- tools::toTitleCase(component_labels)

        # Add total R² labels on top of bars
        total_r_squared <- aggregate(R_squared ~ PC, data = plot_data, FUN = sum)

        # Create plot
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]],
                                                     y = .data[["R_squared"]],
                                                     fill = .data[["Component"]])) +
            ggplot2::geom_col(alpha = 0.8, width = 0.6) +
            ggplot2::geom_text(data = total_r_squared,
                               ggplot2::aes(x = .data[["PC"]],
                                            y = .data[["R_squared"]],
                                            label = sprintf("%.3f", .data[["R_squared"]])),
                               vjust = -0.5, size = 3.5, fontface = "bold",
                               inherit.aes = FALSE) +
            ggplot2::scale_fill_manual(values = available_colors,
                                       labels = component_labels,
                                       name = "Component") +
            ggplot2::labs(
                title = "R-squared Values by Principal Component",
                subtitle = paste("Model:", .getModelText(x[["indep_var"]])),
                x = NULL,
                y = expression(R^2)
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 16, face = "bold", margin = ggplot2::margin(b = 10)),
                plot.subtitle = ggplot2::element_text(size = 12, color = "gray40", margin = ggplot2::margin(b = 20)),
                axis.text.x = ggplot2::element_text(size = 10, color = "gray20"),
                axis.text.y = ggplot2::element_text(size = 10, color = "gray20"),
                axis.title.y = ggplot2::element_text(size = 12, face = "bold", color = "gray20"),
                legend.position = "right",
                legend.title = ggplot2::element_text(size = 11, face = "bold"),
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.5),
                plot.margin = ggplot2::margin(20, 20, 20, 20)
            ) +
            ggplot2::ylim(0, max(total_r_squared[["R_squared"]]) * 1.1)

    } else {
        # Fallback to simple bar plot (original behavior)
        plot_data <- data.frame(
            PC = factor(pc_labels, levels = pc_labels),
            R_squared = x[["r_squared"]]
        )

        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]], y = .data[["R_squared"]])) +
            ggplot2::geom_col(fill = "steelblue", alpha = 0.8, width = 0.6) +
            ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data[["R_squared"]])),
                               vjust = -0.5, size = 3.5, fontface = "bold") +
            ggplot2::labs(
                title = "R-squared Values by Principal Component",
                subtitle = paste("Model:", .getModelText(x[["indep_var"]])),
                x = NULL,
                y = expression(R^2)
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 16, face = "bold", margin = ggplot2::margin(b = 10)),
                plot.subtitle = ggplot2::element_text(size = 12, color = "gray40", margin = ggplot2::margin(b = 20)),
                axis.text.x = ggplot2::element_text(size = 10, color = "gray20"),
                axis.text.y = ggplot2::element_text(size = 10, color = "gray20"),
                axis.title.y = ggplot2::element_text(size = 12, face = "bold", color = "gray20"),
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.5),
                plot.margin = ggplot2::margin(20, 20, 20, 20)
            ) +
            ggplot2::ylim(0, max(plot_data[["R_squared"]]) * 1.1)
    }

    return(p)
}

#' @title Generate Variance Contribution Bar Plot with Component Breakdown
#'
#' @description
#' Creates a bar plot visualization of variance contributions for each principal component,
#' showing how much total dataset variance is explained by the regression model.
#' When available, displays stacked bars showing individual component contributions
#' (cell type, batch/dataset, and interaction effects).
#'
#' @details
#' This function visualizes the variance contribution of each principal component,
#' calculated as the product of PC variance and R-squared values. The variance
#' contribution represents the percentage of total dataset variance explained by
#' the regression model for each PC.
#'
#' When component decomposition is available, the function creates stacked bars with:
#' \itemize{
#'   \item Cell type main effect contribution (blue)
#'   \item Batch or dataset main effect contribution (orange)
#'   \item Interaction effect contribution (green)
#' }
#'
#' The plot subtitle includes the total variance explained across all principal
#' components. PC labels show the individual variance percentage for each component.
#' The function automatically selects appropriate PCA variance values based on
#' analysis type (query-only vs. query+reference).
#'
#' @param x An object of class \code{regressPCObject} containing regression results
#'   from \code{regressPC} function.
#' @param ... Additional arguments passed to the plotting function (currently unused).
#'
#' @keywords internal
#'
#' @return A \code{ggplot2} object representing the variance contribution bar plot
#'   with optional component breakdown.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_fill_manual labs theme_minimal theme element_text element_blank element_line margin ylim
#' @importFrom tools toTitleCase
#'
# Helper function: Variance contribution barplot with component breakdown
plotVarianceContribution <- function(x, ...) {

    # Helper function to get model text
    .getModelText <- function(indep_var) {
        switch(indep_var,
               "cell_type" = "PC ~ Cell Type",
               "cell_type_batch_interaction" = "PC ~ Cell Type * Batch",
               "cell_type_dataset_interaction" = "PC ~ Cell Type * Dataset",
               "Unknown Model")
    }

    # Determine which PCA variance to use based on data type
    if (!is.null(x[["reference_pca_var"]])) {
        # Reference + Query data
        pca_var <- x[["reference_pca_var"]]
    } else {
        # Query only data
        pca_var <- x[["query_pca_var"]]
    }

    # Get PC indices for variance lookup
    pc_indices <- as.numeric(gsub("PC", "", names(x[["var_contributions"]])))

    # Create PC names with variance labels
    pc_names <- names(x[["var_contributions"]])
    pc_labels <- paste0(pc_names, "\n(", sprintf("%.1f%%", pca_var[pc_indices]), ")")

    # Check if component breakdown is available
    if (!is.null(x[["var_contributions_components"]]) && length(x[["var_contributions_components"]]) > 1) {
        # Create stacked bar plot with components

        # Prepare component data
        component_data <- list()
        for (component_name in names(x[["var_contributions_components"]])) {
            component_data[[component_name]] <- data.frame(
                PC = factor(pc_labels, levels = pc_labels),
                Component = component_name,
                Variance_Contribution = x[["var_contributions_components"]][[component_name]],
                stringsAsFactors = FALSE
            )
        }

        # Combine all components
        plot_data <- do.call(rbind, component_data)

        # Create proper component labels and colors
        plot_data[["Component"]] <- factor(plot_data[["Component"]],
                                           levels = names(x[["var_contributions_components"]]))

        # Set up colors based on available components
        component_colors <- c(
            "cell_type" = "#1f77b4",     # Blue
            "batch" = "#ff7f0e",         # Orange
            "dataset" = "#ff7f0e",       # Orange (same as batch)
            "interaction" = "#2ca02c"    # Green
        )

        # Get colors for available components
        available_colors <- component_colors[names(x[["var_contributions_components"]])]

        # Create component labels for legend
        component_labels <- names(x[["var_contributions_components"]])
        component_labels <- gsub("_", " ", component_labels)
        component_labels <- tools::toTitleCase(component_labels)

        # Add total variance contribution labels on top of bars
        total_var_contrib <- aggregate(Variance_Contribution ~ PC, data = plot_data, FUN = sum)

        # Create plot
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]],
                                                     y = .data[["Variance_Contribution"]],
                                                     fill = .data[["Component"]])) +
            ggplot2::geom_col(alpha = 0.8, width = 0.6) +
            ggplot2::geom_text(data = total_var_contrib,
                               ggplot2::aes(x = .data[["PC"]],
                                            y = .data[["Variance_Contribution"]],
                                            label = sprintf("%.2f%%", .data[["Variance_Contribution"]])),
                               vjust = -0.5, size = 3.5, fontface = "bold",
                               inherit.aes = FALSE) +
            ggplot2::scale_fill_manual(values = available_colors,
                                       labels = component_labels,
                                       name = "Component") +
            ggplot2::labs(
                title = "Variance Contribution by Principal Component",
                subtitle = paste("Total variance explained:",
                                 sprintf("%.2f%%", x[["total_variance_explained"]]),
                                 "| Model:", .getModelText(x[["indep_var"]])),
                x = NULL,
                y = "Variance Contribution (%)"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 16, face = "bold", margin = ggplot2::margin(b = 10)),
                plot.subtitle = ggplot2::element_text(size = 12, color = "gray40", margin = ggplot2::margin(b = 20)),
                axis.text.x = ggplot2::element_text(size = 10, color = "gray20"),
                axis.text.y = ggplot2::element_text(size = 10, color = "gray20"),
                axis.title.y = ggplot2::element_text(size = 12, face = "bold", color = "gray20"),
                legend.position = "right",
                legend.title = ggplot2::element_text(size = 11, face = "bold"),
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.5),
                plot.margin = ggplot2::margin(20, 20, 20, 20)
            ) +
            ggplot2::ylim(0, max(total_var_contrib[["Variance_Contribution"]]) * 1.1)

    } else {
        # Fallback to simple bar plot (original behavior)
        plot_data <- data.frame(
            PC = factor(pc_labels, levels = pc_labels),
            Variance_Contribution = x[["var_contributions"]],
            PCA_Variance = pca_var[pc_indices],
            R_squared = x[["r_squared"]]
        )

        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]], y = .data[["Variance_Contribution"]])) +
            ggplot2::geom_col(fill = "darkgreen", alpha = 0.8, width = 0.6) +
            ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f%%", .data[["Variance_Contribution"]])),
                               vjust = -0.5, size = 3.5, fontface = "bold") +
            ggplot2::labs(
                title = "Variance Contribution by Principal Component",
                subtitle = paste("Total variance explained:",
                                 sprintf("%.2f%%", x[["total_variance_explained"]]),
                                 "| Model:", .getModelText(x[["indep_var"]])),
                x = NULL,
                y = "Variance Contribution (%)"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 16, face = "bold", margin = ggplot2::margin(b = 10)),
                plot.subtitle = ggplot2::element_text(size = 12, color = "gray40", margin = ggplot2::margin(b = 20)),
                axis.text.x = ggplot2::element_text(size = 10, color = "gray20"),
                axis.text.y = ggplot2::element_text(size = 10, color = "gray20"),
                axis.title.y = ggplot2::element_text(size = 12, face = "bold", color = "gray20"),
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.5),
                plot.margin = ggplot2::margin(20, 20, 20, 20)
            ) +
            ggplot2::ylim(0, max(plot_data[["Variance_Contribution"]]) * 1.1)
    }

    return(p)
}

#' @title Generate Regression Coefficients Heatmap
#'
#' @description
#' Creates a heatmap visualization of regression coefficients from principal component
#' regression analysis, organized by coefficient type (cell type, batch, interaction)
#' and annotated with significance indicators.
#'
#' @details
#' This function generates a comprehensive heatmap showing regression coefficients
#' for each principal component and model term. The visualization includes:
#' \itemize{
#'   \item Color-coded coefficient values (blue = negative, red = positive)
#'   \item Significance indicators (asterisks) for adjusted p-values below threshold
#'   \item Faceted organization by coefficient category (Cell Type, Batch, Interaction)
#'   \item Clean term labels with proper formatting and reference category information
#' }
#'
#' The function handles different model types automatically:
#' \itemize{
#'   \item Simple cell type models: \code{PC ~ cell_type}
#'   \item Batch interaction models: \code{PC ~ cell_type * batch}
#'   \item Dataset interaction models: \code{PC ~ cell_type * dataset}
#' }
#'
#' Term labels are cleaned and formatted for better readability, with batch/dataset
#' terms converted to consistent "Query Batch" terminology. The plot includes
#' comprehensive subtitle information showing model specification, significance
#' threshold, and reference categories.
#'
#' @param x An object of class \code{regressPCObject} containing regression results
#'   from \code{regressPC} function.
#' @param alpha Numeric value specifying the significance threshold for p-value
#'   adjustment. Default is 0.05.
#' @param coefficients_include Character vector specifying which coefficient types to include.
#'   Options are \code{c("cell_type", "batch", "interaction")}. Default is \code{NULL},
#'   which includes all available coefficient types.
#' @param ... Additional arguments passed to the plotting function (currently unused).
#'
#' @keywords internal
#'
#' @return A \code{ggplot2} object representing the regression coefficients heatmap
#'   with faceted organization and significance annotations.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_gradient2 facet_grid labs theme_minimal theme element_text element_blank element_rect unit scale_y_discrete
#'
# Helper function: Coefficient heatmap
plotCoefficientHeatmap <- function(x, alpha = 0.05, coefficients_include = NULL, ...) {

    # Determine which PCA variance to use based on data type
    if (!is.null(x[["reference_pca_var"]])) {
        # Reference + Query data
        pca_var <- x[["reference_pca_var"]]
    } else {
        # Query only data
        pca_var <- x[["query_pca_var"]]
    }

    # Get PC indices for variance lookup
    pc_names <- names(x[["regression_summaries"]])
    pc_indices <- as.numeric(gsub("PC", "", pc_names))

    # Create PC labels with variance percentages
    pc_labels <- paste0(pc_names, "\n(", sprintf("%.1f%%", pca_var[pc_indices]), ")")

    # Extract coefficient data
    coef_list <- lapply(pc_names, function(pc) {
        coef_data <- x[["regression_summaries"]][[pc]][["coefficients"]]
        coef_data[["PC"]] <- pc
        coef_data[["Term"]] <- rownames(coef_data)
        return(coef_data)
    })

    # Combine into single dataframe
    plot_data <- do.call(rbind, coef_list)

    # Remove intercept for cleaner visualization
    plot_data <- plot_data[!grepl("Intercept", plot_data[["Term"]]), ]

    # Remove any rows with NA values
    plot_data <- plot_data[complete.cases(plot_data), ]

    # Add significance indicator
    plot_data[["Significant"]] <- plot_data[["p.adjusted"]] < alpha

    # Initialize category and clean term names
    plot_data[["Category"]] <- "Cell Type"
    plot_data[["Term_clean"]] <- plot_data[["Term"]]

    # Handle different model types
    if (x[["indep_var"]] == "cell_type_batch_interaction") {

        # Identify main effect cell types (no "batch" and no ":" and no "dataset")
        celltype_main <- !grepl("batch", plot_data[["Term"]]) &
            !grepl(":", plot_data[["Term"]]) &
            !grepl("dataset", plot_data[["Term"]])
        plot_data[["Category"]][celltype_main] <- "Cell Type"

        # Identify main effect batch (starts with "batch" but no ":")
        batch_main <- grepl("^batch", plot_data[["Term"]]) & !grepl(":", plot_data[["Term"]])
        plot_data[["Category"]][batch_main] <- "Batch"

        # Identify interactions (contains ":")
        interactions <- grepl(":", plot_data[["Term"]])
        plot_data[["Category"]][interactions] <- "Cell Type : Batch"

        # Clean up term names
        plot_data[["Term_clean"]] <- gsub("cell_type", "", plot_data[["Term_clean"]])

        # Special handling for batch names
        batch_pattern <- grepl("^batch", plot_data[["Term_clean"]])
        if (any(batch_pattern)) {
            batch_names <- gsub("^batch", "", plot_data[["Term_clean"]][batch_pattern])
            # Remove any remaining "batch" or "sample" prefixes
            batch_names <- gsub("^batch", "", batch_names, ignore.case = TRUE)
            batch_names <- gsub("^sample", "", batch_names, ignore.case = TRUE)
            # Remove Query_ prefix if present (fixes the main effect double Query issue)
            batch_names <- gsub("^Query_", "", batch_names)
            # Remove leading non-alphanumeric characters
            batch_names <- gsub("^[^A-Za-z0-9]+", "", batch_names)
            plot_data[["Term_clean"]][batch_pattern] <- paste("Query Batch", batch_names)
        }

        # Clean up remaining terms
        plot_data[["Term_clean"]] <- gsub("^:", "", plot_data[["Term_clean"]])

        # Special handling for interactions
        interaction_indices <- grepl(":", plot_data[["Term_clean"]])
        if (any(interaction_indices)) {
            interaction_terms <- plot_data[["Term_clean"]][interaction_indices]
            interaction_parts <- strsplit(interaction_terms, ":")
            formatted_interactions <- sapply(interaction_parts, function(parts) {
                if (length(parts) == 2) {
                    cell_type_part <- trimws(parts[1])
                    batch_part <- trimws(parts[2])

                    # Clean batch part thoroughly
                    batch_part <- gsub("^batch", "", batch_part, ignore.case = TRUE)
                    batch_part <- gsub("^sample", "", batch_part, ignore.case = TRUE)
                    batch_part <- gsub("^[^A-Za-z0-9]+", "", batch_part)

                    # Remove Query_ prefix if present (fixes the double Query issue)
                    batch_part <- gsub("^Query_", "", batch_part)
                    batch_part <- trimws(batch_part)

                    # Add "Query Batch" prefix if not already there
                    if (!grepl("Query Batch", batch_part)) {
                        batch_part <- paste("Query Batch", batch_part)
                    }
                    return(paste("Cell Type", cell_type_part, ":", batch_part))
                }
                return(paste(parts, collapse = " : "))
            })
            plot_data[["Term_clean"]][interaction_indices] <- formatted_interactions
        }

    } else if (x[["indep_var"]] == "cell_type_dataset_interaction") {

        # Identify main effect cell types (no "dataset" and no ":")
        celltype_main <- !grepl("dataset", plot_data[["Term"]]) & !grepl(":", plot_data[["Term"]])
        plot_data[["Category"]][celltype_main] <- "Cell Type"

        # Identify main effect dataset (starts with "dataset" but no ":")
        dataset_main <- grepl("^dataset", plot_data[["Term"]]) & !grepl(":", plot_data[["Term"]])
        plot_data[["Category"]][dataset_main] <- "Batch"

        # Identify interactions (contains ":")
        interactions <- grepl(":", plot_data[["Term"]])
        plot_data[["Category"]][interactions] <- "Cell Type : Batch"

        # Clean up term names
        plot_data[["Term_clean"]] <- gsub("cell_type", "", plot_data[["Term_clean"]])

        # Special handling for dataset names - convert to batch terminology
        dataset_pattern <- grepl("^dataset", plot_data[["Term_clean"]])
        if (any(dataset_pattern)) {
            dataset_names <- gsub("^dataset", "", plot_data[["Term_clean"]][dataset_pattern])
            # For dataset interaction, "Query" becomes "Query Batch"
            dataset_names <- gsub("Query", "Query Batch", dataset_names)
            plot_data[["Term_clean"]][dataset_pattern] <- dataset_names
        }

        # Clean up remaining terms
        plot_data[["Term_clean"]] <- gsub("^:", "", plot_data[["Term_clean"]])

        # Special handling for interactions
        interaction_indices <- grepl(":", plot_data[["Term_clean"]])
        if (any(interaction_indices)) {
            interaction_terms <- plot_data[["Term_clean"]][interaction_indices]
            interaction_parts <- strsplit(interaction_terms, ":")
            formatted_interactions <- sapply(interaction_parts, function(parts) {
                if (length(parts) == 2) {
                    cell_type_part <- trimws(parts[1])
                    dataset_part <- trimws(parts[2])

                    # Clean dataset part
                    dataset_part <- gsub("^dataset", "", dataset_part, ignore.case = TRUE)
                    dataset_part <- trimws(dataset_part)

                    # Convert "Query" to "Query Batch"
                    if (dataset_part == "Query") {
                        dataset_part <- "Query Batch"
                    }

                    return(paste("Cell Type", cell_type_part, ":", dataset_part))
                }
                return(paste(parts, collapse = " : "))
            })
            plot_data[["Term_clean"]][interaction_indices] <- formatted_interactions
        }

    } else if (x[["indep_var"]] == "cell_type") {
        # For cell type only model
        plot_data[["Term_clean"]] <- gsub("cell_type", "", plot_data[["Term_clean"]])
    }

    # Filter coefficients based on coefficients_include parameter
    if (!is.null(coefficients_include)) {
        # Map coefficient types to categories
        category_mapping <- list(
            "cell_type" = "Cell Type",
            "batch" = "Batch",
            "interaction" = "Cell Type : Batch"
        )

        # Get categories to include
        categories_to_include <- unlist(category_mapping[coefficients_include])

        # Filter plot data
        plot_data <- plot_data[plot_data[["Category"]] %in% categories_to_include, ]

        # If no data remains after filtering, stop with informative message
        if (nrow(plot_data) == 0) {
            stop("No coefficients found for the specified types: ",
                 paste(coefficients_include, collapse = ", "))
        }
    }

    # Create proper factor ordering for facets (only include categories present in data)
    available_categories <- unique(plot_data[["Category"]])
    all_category_levels <- c("Cell Type", "Batch", "Cell Type : Batch")
    category_levels <- all_category_levels[all_category_levels %in% available_categories]

    plot_data[["Category"]] <- factor(plot_data[["Category"]], levels = category_levels)

    # Create PC factor with labels
    existing_pc_names <- unique(plot_data[["PC"]])
    existing_pc_indices <- match(existing_pc_names, pc_names)
    existing_pc_labels <- pc_labels[existing_pc_indices]

    plot_data[["PC_labeled"]] <- factor(plot_data[["PC"]],
                                        levels = existing_pc_names,
                                        labels = existing_pc_labels)

    # Create y-axis ordering - special handling for interactions
    if (x[["indep_var"]] %in% c("cell_type_batch_interaction", "cell_type_dataset_interaction")) {
        # For interactions, we want to order by batch first, then by cell type
        interaction_data <- plot_data[plot_data[["Category"]] == "Cell Type : Batch", ]

        if (nrow(interaction_data) > 0) {
            # Extract batch and cell type parts from interaction terms
            interaction_parts <- strsplit(interaction_data[["Term_clean"]], " : ")
            batch_parts <- sapply(interaction_parts, function(x) if(length(x) == 2) trimws(x[2]) else "")
            cell_type_parts <- sapply(interaction_parts, function(x) if(length(x) == 2) gsub("Cell Type ", "", trimws(x[1])) else "")

            # Create ordering data frame with unique terms only
            ordering_df <- data.frame(
                Term_clean = interaction_data[["Term_clean"]],
                batch_part = batch_parts,
                cell_type_part = cell_type_parts,
                stringsAsFactors = FALSE
            )

            # Remove duplicates
            ordering_df <- unique(ordering_df)

            # Order by batch first, then by cell type
            ordering_df <- ordering_df[order(ordering_df[["batch_part"]], ordering_df[["cell_type_part"]]), ]

            # Get the ordered interaction terms
            ordered_interaction_terms <- ordering_df[["Term_clean"]]

            # Get unique non-interaction terms
            non_interaction_terms <- unique(plot_data[plot_data[["Category"]] != "Cell Type : Batch", "Term_clean"])

            # Combine all unique terms in desired order
            all_terms <- unique(c(non_interaction_terms, ordered_interaction_terms))

            # Create factor with custom ordering
            plot_data[["Term_factor"]] <- factor(plot_data[["Term_clean"]], levels = all_terms)
        } else {
            # No interactions, use default ordering
            plot_data[["Term_factor"]] <- factor(plot_data[["Term_clean"]])
        }
    } else {
        # For non-interaction models, use default ordering
        plot_data[["Term_factor"]] <- factor(plot_data[["Term_clean"]])
    }

    # Check if Batch category has only one term and add padding if needed
    if (x[["indep_var"]] %in% c("cell_type_batch_interaction", "cell_type_dataset_interaction") &&
        "Batch" %in% available_categories) {
        batch_terms <- unique(plot_data[plot_data[["Category"]] == "Batch", "Term_clean"])
        if (length(batch_terms) == 1) {
            # Add a dummy row for better spacing
            dummy_row <- plot_data[plot_data[["Category"]] == "Batch", ][1, ]
            dummy_row[["Term_clean"]] <- ""
            dummy_row[["Term_factor"]] <- factor("", levels = c(levels(plot_data[["Term_factor"]]), ""))
            dummy_row[["coef"]] <- NA
            dummy_row[["Significant"]] <- FALSE
            plot_data <- rbind(plot_data, dummy_row)
        }
    }

    # Create main subtitle
    model_text <- switch(x[["indep_var"]],
                         "cell_type" = "PC ~ Cell Type",
                         "cell_type_batch_interaction" = "PC ~ Cell Type * Batch",
                         "cell_type_dataset_interaction" = "PC ~ Cell Type * Batch")

    subtitle_text1 <- paste("Model:", model_text, "| * indicates significance at alpha =", alpha)

    # Add information about coefficient filtering if applied
    if (!is.null(coefficients_include)) {
        coefficient_text <- paste("Showing coefficients:", paste(coefficients_include, collapse = ", "))
        subtitle_text1 <- paste(subtitle_text1, "|", coefficient_text)
    }

    # Create reference category subtitle
    reference_text_parts <- c()

    # Always include cell type reference
    if (!is.null(x[["reference_cell_type"]])) {
        reference_text_parts <- c(reference_text_parts, paste("Cell Type:", x[["reference_cell_type"]]))
    }

    # Include batch/dataset reference based on model type
    if (x[["indep_var"]] == "cell_type_batch_interaction") {
        # For query+reference with multiple batches OR query-only with batches
        if (!is.null(x[["reference_batch"]])) {
            reference_text_parts <- c(reference_text_parts, paste("Batch: Query Batch", x[["reference_batch"]]))
        } else if (!is.null(x[["reference_pca_var"]])) {
            # Query+reference with multiple batches (Reference is reference category)
            reference_text_parts <- c(reference_text_parts, "Batch: Reference")
        }
    } else if (x[["indep_var"]] == "cell_type_dataset_interaction") {
        # Query + Reference (Reference is always the reference category)
        reference_text_parts <- c(reference_text_parts, "Batch: Reference")
    }

    # Create subtitle
    if (length(reference_text_parts) > 0) {
        subtitle_text2 <- paste("Reference categories -", paste(reference_text_parts, collapse = "; "))
    } else {
        subtitle_text2 <- "Reference category is first alphabetically"
    }

    # Create plot with facets
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC_labeled"]],
                                                 y = .data[["Term_factor"]],
                                                 fill = .data[["coef"]])) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::geom_point(data = plot_data[plot_data[["Significant"]] & !is.na(plot_data[["coef"]]), ],
                            ggplot2::aes(x = .data[["PC_labeled"]],
                                         y = .data[["Term_factor"]]),
                            shape = 8, size = 2, color = "black") +
        ggplot2::scale_fill_gradient2(
            low = "blue", mid = "white", high = "red",
            midpoint = 0, name = "Coefficient",
            na.value = "white"
        ) +
        ggplot2::facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
        ggplot2::labs(
            title = "Regression Coefficients Heatmap",
            subtitle = paste(subtitle_text1, "\n", subtitle_text2),
            x = NULL,
            y = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            axis.text.x = ggplot2::element_text(size = 10),
            plot.title = ggplot2::element_text(size = 16, face = "bold", margin = ggplot2::margin(b = 10)),
            plot.subtitle = ggplot2::element_text(size = 11, color = "gray40", margin = ggplot2::margin(b = 20)),
            legend.position = "right",
            legend.title = ggplot2::element_text(size = 11, face = "bold"),
            panel.grid = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(20, 20, 20, 20),
            strip.text.y = ggplot2::element_text(size = 10, face = "bold", color = "black"),
            strip.background = ggplot2::element_rect(fill = "gray92", color = "black", linewidth = 0.5),
            panel.spacing = ggplot2::unit(0.8, "lines"),
            panel.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.5)
        ) +
        ggplot2::scale_y_discrete(labels = function(x) ifelse(x == "", "", x))

    return(p)
}
