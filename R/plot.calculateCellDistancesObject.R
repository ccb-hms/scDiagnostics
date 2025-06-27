#' @title Plot Distance Density Comparison for a Specific Cell Type and Selected Cells
#'
#' @description
#' The S3 plot method plots the density functions for the reference data and the distances from a specified query cells
#' to all reference cell within a specified cell type.
#'
#' @details
#' The S3 plot method first checks if the specified cell type and cell names are present in the object. If the
#' specified cell type or cell name is not found, an error is thrown. It then extracts the distances within the reference dataset
#' and the distances from the specified query cell to the reference cells The function creates a density plot using \code{ggplot2}
#' to compare the distance distributions. The density plot will show two distributions: one for the pairwise distances within the
#' reference dataset and one for the distances from the specified query cell to each reference cell. These distributions are
#' plotted in different colors to visually assess how similar the query cell is to the reference cells of the specified cell type.
#'
#' @param x A list containing the distance data computed by \code{calculatecellDistances}.
#' @param ref_cell_type A string specifying the reference cell type.
#' @param cell_names A string specifying the query cell name for which to plot the distances.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @keywords internal
#'
#' @return The S3 plot method returns a \code{ggplot} density plot comparing the reference distances and the distances from the specified cell to the reference cells.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateCellDistances}}
#'
#' @rdname calculateCellDistances
#'
# Function to plot density functions for the reference data and the specified cell
plot.calculateCellDistancesObject <- function(x, ref_cell_type, cell_names, ...) {

    # Check if cell type is available
    if(length(ref_cell_type) != 1 || !(ref_cell_type %in% names(x)))
        stop("The specified \'ref_cell_type\' is not available.")

    # Filter distance data for the specified cell type
    distance_data <- x[[ref_cell_type]]

    # Check if cells are available in data for that cell type
    if(!all(cell_names %in% rownames(
        distance_data[["query_to_ref_distances"]])))
        stop("One or more specified 'cell_names' are not available for that cell type.")

    # Extract distances within the reference dataset
    ref_distances <- distance_data[["ref_distances"]]

    # Initialize an empty list to store data frames for each cell
    plot_data_list <- vector("list", length = length(cell_names))
    names(plot_data_list) <- cell_names

    # Loop through each cell to create the combined data frame
    for(s in cell_names) {
        # Extract distances for the current cell
        cell_distances <- distance_data[["query_to_ref_distances"]][s, ]

        # Create a data frame for the current cell and reference distances
        cell_data <- data.frame(cell = s, Distance = cell_distances,
                                Distance_Type = "Query")
        ref_data <- data.frame(cell = s, Distance = ref_distances,
                               Distance_Type = "Reference")

        # Combine the reference and cell data frames
        combined_data <- rbind(ref_data, cell_data)

        # Append the combined data frame to the list
        plot_data_list[[s]] <- combined_data
    }

    # Combine all data frames into one data frame
    plot_data <- do.call(rbind, plot_data_list)

    # Keep order of cell names
    plot_data[["Query"]] <- factor(plot_data[["cell"]], levels = cell_names)

    # Plot density comparison with facets for each cell
    density_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(
        x = .data[["Distance"]], fill = .data[["Distance_Type"]])) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::labs(title = NULL,
                      x = "Distance", y = "Density",
                      fill = "Distance Type") +
        ggplot2::facet_wrap(~ .data[["Query"]], scales = "free_y") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = "white",
                                                     color = "black",
                                                     linewidth = 0.5),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(
                color = "gray",
                linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14,
                                               face = "bold",
                                               hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10),
            legend.position = "bottom")
    return(density_plot)
}







