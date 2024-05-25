#' @title Create Faceted Scatter Plots for Specified PC Combinations From \code{detectAnomaly} Object
#'
#' @description This function generates faceted scatter plots for specified principal component (PC) combinations
#' within an anomaly detection object. It allows visualization of the relationship between specified
#' PCs and highlights anomalies detected by the Isolation Forest algorithm.
#'
#' @details The function extracts the specified PCs from the given anomaly detection object and generates
#' scatter plots for each pair of PCs. It uses \code{ggplot2} to create a faceted plot where each facet represents
#' a pair of PCs. Anomalies are highlighted in red, while normal points are shown in black.
#'
#' @param anomaly_object A list object containing the anomaly detection results from the \code{detectAnomaly} function. 
#' Each element of the list should correspond to a cell type and contain \code{reference_pca_subset}, \code{query_pca_subset}, 
#' \code{var_explained}, and \code{anomaly}.
#' @param cell_type A character string specifying the cell type for which the plots should be generated. This should
#' be a name present in \code{anomaly_object}.
#' @param pc_subset A numeric vector specifying the indices of the PCs to be included in the plots. If NULL, all PCs
#' in \code{reference_pca_subset} will be included.
#' @param ... Additional arguments.
#' 
#' @return A ggplot object displaying the faceted scatter plots for the specified PC combinations.
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
#' ref_var <- getTopHVGs(ref_data, n = 2000)
#' query_var <- getTopHVGs(query_data, n = 2000)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#' 
#' # Run PCA on the reference data
#' ref_data_subset <- runPCA(ref_data_subset, ncomponents = 50)
#' 
#' # Store PCA anomaly data and plots
#' anomaly_output <- detectAnomaly(reference_data = ref_data_subset, 
#'                                 query_data = query_data_subset, 
#'                                 reference_labels = ref_data$reclustered.broad, 
#'                                 query_labels = query_data$labels,
#'                                 n_components = 10,
#'                                 n_tree = 500,
#'                                 anomaly_treshold = 0.5,
#'                                 store_plots = TRUE) 
#' 
#' # Plot the output for a cell type
#' plot(anomaly_object = anomaly_output, cell_type = "CD8", pc_subset = c(1:5))
#' 
# Function to create faceted scatter plots for specified PC combinations
plot.detectAnomaly <- function(anomaly_object, cell_type, pc_subset = NULL, ...) {
    
    # Check input for cell type
    if(!(cell_type %in% names(anomaly_object)))
        stop("\'cell_type\' is not available in \'anomaly_object\'.")
    
    # Check input for pc_subset
    if(!is.null(pc_subset)){
        if(!all(pc_subset %in% 1:ncol(anomaly_object[[cell_type]]$reference_pca_subset)))
            stop("\'pc_subset\' is out of range.")
    } else{
        pc_subset <- 1:ncol(anomaly_object[[cell_type]]$reference_pca_subset)
    }
    
    # Filter data to include only specified PCs
    data_subset <- anomaly_object[[cell_type]]$query_pca_subset[, pc_subset, drop = FALSE]
    
    # Modify column names to include percentage of variance explained
    colnames(data_subset) <- paste0("PC", pc_subset, 
                                    " (", sprintf("%.1f%%", anomaly_object[[cell_type]]$var_explained[pc_subset] * 100), ")")
    
    # Create all possible pairs of specified PCs
    pc_names <- colnames(data_subset)
    pairs <- expand.grid(x = pc_names, y = pc_names) |>
        dplyr::filter(x != y)
    
    # Create a new data frame with all possible pairs of specified PCs
    data_pairs <- pairs |>
        dplyr::mutate(pair_data = purrr::map2(x, y, ~ data_subset |>
                                                  dplyr::select(dplyr::all_of(c(.x, .y))) |>
                                                  dplyr::rename(x_value = .x, y_value = .y) |>
                                                  dplyr::mutate(x = .x, y = .y))) |>
        dplyr::select(-x, -y) |>
        tidyr::unnest(cols = pair_data)
    
    # Remove points on the diagonal
    data_pairs <- data_pairs |>
        dplyr::filter(x_value != y_value)
    
    # Remove redundant data (to avoid duplicated plots)
    data_pairs <- data_pairs[as.numeric(substr(data_pairs$x, 3, 3)) < as.numeric(substr(data_pairs$y, 3, 3)),]
    
    # Add anomalies vector to data_pairs dataframe
    data_pairs$anomaly <- rep(anomaly_object[[cell_type]]$anomaly, choose(length(pc_subset), 2))
    
    # Create the ggplot object with facets
    plot <- ggplot2::ggplot(data_pairs, ggplot2::aes(x = x_value, y = y_value, color = factor(anomaly))) +
        ggplot2::geom_point(size = 2) + 
        ggplot2::scale_color_manual(values = c("black", "red"), labels = c("Normal", "Anomaly")) + 
        ggplot2::facet_grid(rows = ggplot2::vars(y), cols = ggplot2::vars(x), scales = "free") +
        ggplot2::theme_minimal() +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey85", color = "grey70"),   
                       strip.text = ggplot2::element_text(size = 10, face = "bold", color = "black"), 
                       axis.title = ggplot2::element_blank(),        
                       axis.text = ggplot2::element_text(size = 10), 
                       panel.grid = ggplot2::element_blank(),        
                       panel.background = ggplot2::element_rect(fill = "white", color = "black"), 
                       legend.position = "right",          
                       plot.title = ggplot2::element_text(size = 14, hjust = 0.5), 
                       plot.background = ggplot2::element_rect(fill = "white")) + 
        ggplot2::labs(title = paste0("Isolation Forest Anomaly Plot: ", cell_type), color = "iForest Type")
    print(plot)
}
