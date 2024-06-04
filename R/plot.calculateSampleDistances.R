#' @title Plot Distance Density Comparison for a Specific Cell Type and Selected Samples
#'
#' @description This function plots the density functions for the reference data and the distances from a specified query samples 
#' to all reference samples within a specified cell type.
#'
#' @details The function first checks if the specified cell type and sample names are present in the \code{x}. If the 
#' specified cell type or sample name is not found, an error is thrown. It then extracts the distances within the reference dataset 
#' and the distances from the specified query sample to the reference samples. The function creates a density plot using \code{ggplot2} 
#' to compare the distance distributions. The density plot will show two distributions: one for the pairwise distances within the 
#' reference dataset and one for the distances from the specified query sample to each reference sample. These distributions are 
#' plotted in different colors to visually assess how similar the query sample is to the reference samples of the specified cell type.
#'
#' @param x A list containing the distance data computed by \code{calculateSampleDistances}.
#' @param ref_cell_type A string specifying the reference cell type.
#' @param sample_names A string specifying the query sample name for which to plot the distances.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot2 density plot comparing the reference distances and the distances from the specified sample to the reference samples.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{calculateSampleDistances}}
#' 
#' @rdname calculateSampleDistances
#'
# Function to plot density functions for the reference data and the specified sample
plot.calculateSampleDistances <- function(x, ref_cell_type, sample_names, ...) {
    
    # Check if cell type is available
    if(length(ref_cell_type) != 1 || !(ref_cell_type %in% names(x)))
        stop("The specified \'ref_cell_type\' is not available.")
    
    # Filter distance data for the specified cell type
    cell_distances <- x[[ref_cell_type]]
    
    # Check if samples are available in data for that cell type
    if(!all(sample_names %in% rownames(cell_distances$query_to_ref_distances)))
        stop("One or more specified 'sample_names' are not available for that cell type.")
    
    # Extract distances within the reference dataset
    ref_distances <- cell_distances$ref_distances
    
    # Initialize an empty list to store data frames for each sample
    plot_data_list <- list()
    
    # Loop through each sample to create the combined data frame
    for(s in sample_names) {
        # Extract distances for the current sample
        sample_distances <- cell_distances$query_to_ref_distances[s, ]
        
        # Create a data frame for the current sample and reference distances
        sample_data <- data.frame(Sample = s, Distance = sample_distances, Distance_Type = "Sample")
        ref_data <- data.frame(Sample = s, Distance = ref_distances, Distance_Type = "Reference")
        
        # Combine the reference and sample data frames
        combined_data <- rbind(ref_data, sample_data)
        
        # Append the combined data frame to the list
        plot_data_list[[s]] <- combined_data
    }
    
    # Combine all data frames into one data frame
    plot_data <- do.call(rbind, plot_data_list)
    
    # Keep order of sample names
    plot_data$Sample <- factor(plot_data$Sample, levels = sample_names)
    
    # Plot density comparison with facets for each sample
    density_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Distance, fill = Distance_Type)) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::labs(title = paste("Distance Density Comparison for Cell Type:", ref_cell_type),
                      x = "Distance", y = "Density") +
        ggplot2::scale_fill_manual(name = "Distance Type", values = c("Reference" = "blue", "Sample" = "red")) +
        ggplot2::facet_wrap(~ Sample, scales = "free_y", labeller = ggplot2::labeller(Sample = label_parsed)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = "lightgrey", color = "grey50"),
            strip.text = ggplot2::element_text(color = "grey20", size = 10, face = "bold"),
            panel.grid.major = ggplot2::element_line(color = "grey90", linetype = "dashed"),
            panel.grid.minor = ggplot2::element_line(color = "grey95", linetype = "dashed")
        )
    
    # Print the density plot
    print(density_plot)
}







