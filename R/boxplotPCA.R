#' @title Plot Principal Components for Different Cell Types
#'
#' @description This function generates a \code{ggplot2} boxplot visualization of principal components (PCs) for different 
#' cell types across two datasets (query and reference).
#'
#' @details
#' The function \code{boxplotPCA} is designed to provide a visualization of principal component analysis (PCA) results. It projects 
#' the query dataset onto the principal components obtained from the reference dataset. The results are then visualized 
#' as boxplots, grouped by cell types and datasets (query and reference). This allows for a comparative analysis of the 
#' distributions of the principal components across different cell types and datasets. The function internally calls \code{projectPCA} 
#' to perform the PCA projection. It then reshapes the output data into a long format suitable for ggplot2 plotting. 
#' The color scheme is automatically determined using the \code{RColorBrewer} package, ensuring a visually distinct and appealing plot.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param n_components An integer specifying the number of principal components to use for projection. Defaults to 10. 
#' Must be less than or equal to the number of components available in the reference PCA.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} 
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} 
#' that identifies the cell types.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#'
#' @return A ggplot object representing the boxplots of specified principal components for the given cell types and datasets.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load required libraries
#' library(scRNAseq)
#' library(scuttle)
#' library(SingleR)
#' library(scran)
#' library(scater)
#'
#' # Load data (replace with your data loading)
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#' 
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#' 
#' # log transform datasets
#' ref_data <- scuttle::logNormCounts(ref_data)
#' query_data <- scuttle::logNormCounts(query_data)
#' 
#' # Get cell type scores using SingleR (or any other cell type annotation method)
#' scores <- SingleR::SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' 
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#' 
#' # Selecting highly variable genes (can be customized by the user)
#' ref_var <- scran::getTopHVGs(ref_data, n = 2000)
#' query_var <- scran::getTopHVGs(query_data, n = 2000)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#'
#' # Run PCA on the reference data (assumed to be prepared)
#' ref_data_subset <- runPCA(ref_data_subset)
#'
#' pc_plot <- boxplotPCA(query_data_subset, ref_data_subset,
#'                       n_components = 10,
#'                       cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                       query_cell_type_col = "labels", 
#'                       ref_cell_type_col = "reclustered.broad", 
#'                       pc_subset = c(1:5))
#' pc_plot
#' 
#' @importFrom stats approxfun cancor density setNames
#' @importFrom utils combn
#'                          
# Function to plot PC for different cell types
boxplotPCA <- function(query_data, reference_data, 
                       n_components = 10, 
                       cell_types = NULL,
                       query_cell_type_col = NULL, 
                       ref_cell_type_col = NULL, 
                       pc_subset = c(1:5)){
    
    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data, reference_data = reference_data, 
                             n_components = n_components, 
                             query_cell_type_col = query_cell_type_col, 
                             ref_cell_type_col = ref_cell_type_col)
    
    # Create the long format data frame manually
    pca_output <- pca_output[!is.na(pca_output$cell_type),]
    if(!is.null(cell_types)){
        if(all(cell_types %in% pca_output$cell_type)){
            pca_output <- pca_output[which(pca_output$cell_type %in% cell_types),]
        } else{
            stop("One or more of the specified \'cell_types\' are not available.")
        }
    }
    pca_long <- data.frame(PC = rep(paste0("pc", pc_subset), each = nrow(pca_output)),
                           Value = unlist(c(pca_output[, pc_subset])),
                           dataset = rep(pca_output$dataset, length(pc_subset)),
                           cell_type = rep(pca_output$cell_type, length(pc_subset)))
    pca_long$PC <- toupper(pca_long$PC)
    
    # Create a new variable representing the combination of cell type and dataset
    pca_long$cell_type_dataset <- paste(pca_long$dataset, pca_long$cell_type, sep = " ")
    
    # Define the order of cell type and dataset combinations
    order_combinations <- paste(rep(c("Reference", "Query"), length(unique(pca_long$cell_type))),
                                rep(sort(unique(pca_long$cell_type)), each = 2))
    
    # Reorder the levels of cell type and dataset factor
    pca_long$cell_type_dataset <- factor(pca_long$cell_type_dataset, levels = order_combinations)
    
    # Define the colors for cell types
    color_palette <- RColorBrewer::brewer.pal(length(order_combinations), "Paired")
    color_mapping <- setNames(color_palette, order_combinations)
    cell_type_colors <- color_mapping[order_combinations]
    
    # Create the ggplot
    plot <- ggplot2::ggplot(pca_long, aes(x = cell_type, y = Value, fill = cell_type_dataset)) +
        ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7) + 
        ggplot2::facet_wrap(~ PC, scales = "free") +
        ggplot2::scale_fill_manual(values = cell_type_colors, name = "Cell Types") + 
        ggplot2::labs(x = "", y = "Value") +  
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "right",  
                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),  
                       axis.title = ggplot2::element_text(size = 14), 
                       strip.text = ggplot2::element_text(size = 12, face = "bold"), 
                       panel.grid.major = ggplot2::element_line(color = "grey", linetype = "dotted", linewidth = 0.7),  
                       panel.grid.minor = ggplot2::element_blank(),  
                       panel.border = ggplot2::element_blank(),  
                       strip.background = ggplot2::element_rect(fill = "lightgrey", color = "grey", linewidth = 0.5),  
                       plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5))
    
    # Return the plot
    return(plot)
}


