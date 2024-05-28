#' @title Plot Density of Probabilities for Cell Type Classification
#'
#' @description This function generates a density plot showing the distribution of probabilities for each sample of belonging to 
#' either the reference or query dataset for each cell type.
#'
#' @details This function creates a density plot to visualize the distribution of probabilities for each sample belonging to the 
#' reference or query dataset for each cell type. It utilizes the ggplot2 package for plotting.
#'
#' @param x An object of class \code{nearestNeighbotDiagnostics} containing the probabilities calculated by the \code{\link{nearestNeighborDiagnostics}} function.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types in \code{x} will be plotted. Default is NULL.
#' @param prob_type A character string specifying the type of probability to plot. Must be either "query" or "reference". Default is "query".
#' @param ... Additional arguments to be passed to \code{\link[ggplot2]{geom_density}}.
#'
#' @return A ggplot2 density plot.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{nearestNeighborDiagnostics}}
#' 
#' @examples
#' # Load necessary library
#' library(scRNAseq)
#' library(scuttle)
#' library(scran)
#' library(SingleR)
#' library(scater)
#'
#' # Load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#' 
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#' 
#' # log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#' 
#' # Get cell type scores using SingleR (or any other cell type annotation method)
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' 
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#' 
#' # Selecting highly variable genes (can be customized by the user)
#' ref_var <- getTopHVGs(ref_data, n = 500)
#' query_var <- getTopHVGs(query_data, n = 500)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#' 
#' # Run PCA on the reference data
#' ref_data_subset <- runPCA(ref_data_subset)
#'
#' # Project the query data onto PCA space of reference
#' nn_output <- nearestNeighborDiagnostics(query_data_subset, ref_data_subset,
#'                                         n_neighbor = 15, 
#'                                         n_components = 10,
#'                                         pc_subset = c(1:10),
#'                                         query_cell_type_col = "labels", 
#'                                         ref_cell_type_col = "reclustered.broad")
#' 
#' # Plot output
#' plot(nn_output, cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'      prob_type = "query")
#' 
# Function to plot probabilities of each sample of belonging to reference or query dataset for each cell type
plot.nearestNeighborDiagnostics <- function(x, cell_types = NULL,
                                            prob_type = c("query", "reference")[1], ...) {
    
    # Check input for probability type
    if(!(prob_type %in% c("query", "reference")))
        stop("\'prob_type\' must be one of \'query\' or \'reference\'.")
    
    # Convert probabilities to data frame
    probabilities_df <- do.call(rbind, lapply(names(x), function(ct) {
        data.frame(cell_types = ct, 
                   probability = x[[ct]][[ifelse(prob_type == "reference", "prob_ref", "prob_query")]])
    }))
    
    if(!is.null(cell_types)){
        
        if(!all(cell_types %in% names(x)))
            stop("One or more of the \'cell_types'\ is not available.")
        
        # Subset cell types
        probabilities_df <- probabilities_df[probabilities_df$cell_types %in% cell_types,]
    }
    
    # Create density plot
    density_plot <- ggplot2::ggplot(probabilities_df, ggplot2::aes(x = probability, fill = cell_types)) +
        ggplot2::geom_density(alpha = 0.7) +
        ggplot2::labs(x = "Probability", y = "Density", title = "Density Plot of Probabilities") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = "none",
            strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
            strip.text = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::facet_wrap(~cell_types, scales = "free", labeller = ggplot2::labeller(cell_types = label_value))
    if(length(unique(probabilities_df$cell_types)) > 2)
        density_plot <- density_plot + 
        ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(n = nlevels(as.factor(probabilities_df$cell_types)), 
                                                                     name = "Set1")) 
    
    return(density_plot)
}
