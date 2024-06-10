#' @title Plot Cosine Similarities Between Samples and PCs
#'
#' @description 
#' This function creates a heatmap plot to visualize the cosine similarities between samples and principal components (PCs).
#'
#' @details 
#' This function reshapes the input data frame to create a long format suitable for plotting as a heatmap. It then
#' creates a heatmap plot using ggplot2, where the x-axis represents the PCs, the y-axis represents the samples, and the
#' color intensity represents the cosine similarity values.
#'
#' @param x An object of class 'calculateSampleSimilarityPCA' containing a dataframe of cosine similarity values 
#' between samples and PCs.
#' @param pc_subset A numeric vector specifying the subset of principal components to include in the plot (default: c(1:5)).
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object representing the cosine similarity heatmap.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{calculateSampleSimilarityPCA}}
#' 
#' @rdname calculateSampleSimilarityPCA
#'
# Function to plot cosine similarities between samples and PCs
plot.calculateSampleSimilarityPCA <- function(x, pc_subset = c(1:5), ...){
    
    # Subset data
    x <- x[, paste0("PC", pc_subset)]
    
    # Initialize empty vectors for reshaped data
    sample_names <- c()
    pc_names <- c()
    cosine_values <- c()
    
    # Loop through the data frame to manually reshape it
    for (sample in rownames(x)) {
        for (pc in colnames(x)) {
            sample_names <- c(sample_names, sample)
            pc_names <- c(pc_names, pc)
            cosine_values <- c(cosine_values, x[sample, pc])
        }
    }
    
    # Create a data frame with the reshaped data
    cosine_long <- data.frame(Sample = factor(sample_names, levels = rev(rownames(x))), 
                              PC = pc_names, CosineSimilarity = cosine_values)
    
    # Create the heatmap plot
    plot <- ggplot2::ggplot(cosine_long, ggplot2::aes(x = PC, y = Sample, fill = CosineSimilarity)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", CosineSimilarity)), size = 3) +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                             limits = c(-1, 1), space = "Lab", name = "Cosine Similarity") +
        ggplot2::labs(title = "Cosine Similarity Heatmap", x = "", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       plot.title = ggplot2::element_text(hjust = 0.5))
    return(plot)
}

