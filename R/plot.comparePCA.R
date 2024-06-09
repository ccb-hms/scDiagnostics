#' @title Plot Heatmap of Cosine Similarities Between Principal Components
#' 
#' @description This function generates a heatmap to visualize the cosine similarities between 
#' principal components from the output of the `comparePCA` function.
#' 
#' @details The function converts the input matrix into a long-format data frame 
#' suitable for plotting with `ggplot2`. The rows in the heatmap are ordered in 
#' reverse to match the conventional display format. The heatmap uses a blue-white-red 
#' color gradient to represent cosine similarity values, where blue indicates negative 
#' similarity, white indicates zero similarity, and red indicates positive similarity.
#' 
#' @param x A numeric matrix output from the `comparePCA` function, representing 
#' cosine similarities between query and reference principal components.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object representing the heatmap of cosine similarities.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{comparePCA}}
#' 
#' @rdname comparePCA
#' 
# Function to produce the heatmap from output of comparePCA function
plot.comparePCA <- function(x, ...){
    
    # Convert the matrix to a data frame
    similarity_df <- data.frame(
        Ref = factor(rep(rownames(x), each = ncol(x)), levels = rev(rownames(x))),
        Query = rep(colnames(x), times = nrow(x)),
        value = as.vector(x))
    
    # Create the heatmap
    pc_plot <- ggplot2::ggplot(similarity_df, ggplot2::aes(x = Query, y = Ref, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(aes(label = sprintf("%.2f", value)), size = 3) +  # Add text with rounded values
        ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                      midpoint = 0, limit = c(min(x, -0.5), max(x, 0.5)), space = "Lab", 
                                      name = "Cosine Similarity") +
        ggplot2::theme_minimal() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, 
                                                           size = 12, hjust = 1)) +
        ggplot2::labs(x = "", y = "", 
                      title = "Heatmap of Cosine Similarities Between PCs")
    return(pc_plot)
}
