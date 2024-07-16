#' @title Density Plot of Wasserstein Distances Under Null Distribution
#'
#' @description 
#' This function generates a density plot of Wasserstein distances under the null hypothesis that the two datasets
#' come from the same distribution. It computes the null distribution of Wasserstein distances and compares 
#' it to the distances between reference and query data.
#'
#' @details
#' This function first projects the query data onto the PCA space defined by the reference data. It then computes the 
#' Wasserstein distances between randomly sampled subsets of the reference data to generate a null distribution. The 
#' function also computes the Wasserstein distances between the reference and query data. Finally, it visualizes these 
#' distances using a density plot with vertical lines indicating the significance threshold and the reference-query distance.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is 1:5.
#' @param n_resamples An integer specifying the number of resamples to use for generating the null distribution. Default is 300.
#' @param alpha A numeric value specifying the significance level for thresholding. Default is 0.05.
#'
#' @return A ggplot2 object representing the density plot of Wasserstein distances under the null distribution.
#' 
#' @references
#' Schuhmacher, D., Bernhard, S., & Book, M. (2019). "A Review of Approximate Transport in Machine Learning". 
#' In *Journal of Machine Learning Research* (Vol. 20, No. 117, pp. 1-61).
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#' 
#' # Extract CD4 cells
#' ref_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
#' query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]
#' 
#' # Selecting highly variable genes (can be customized by the user)
#' ref_top_genes <- scran::getTopHVGs(ref_data_subset, n = 500)
#' query_top_genes <- scran::getTopHVGs(query_data_subset, n = 500)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_top_genes, query_top_genes)
#' ref_data_subset <- ref_data_subset[common_genes,]
#' query_data_subset <- query_data_subset[common_genes,]
#' 
#' # Run PCA on reference data
#' ref_data_subset <- scater::runPCA(ref_data_subset)
#' 
#' # Compute Wasserstein null distribution using reference data and observed distance with query data
#' wasserstein_plot <- plotWassersteinDistance(query_data = query_data_subset,
#'                                             reference_data = ref_data_subset, 
#'                                             query_cell_type_col = "expert_annotation", 
#'                                             ref_cell_type_col = "expert_annotation", 
#'                                             pc_subset = 1:5,
#'                                             n_resamples = 300,
#'                                             alpha = 0.05)
#' wasserstein_plot
#'  
#' @importFrom stats quantile
#'  
# Function to generate density of Wasserstein distances under null distribution
plotWassersteinDistance <- function(query_data,
                                    reference_data, 
                                    ref_cell_type_col, 
                                    query_cell_type_col, 
                                    pc_subset = 1:5,
                                    n_resamples = 300,
                                    alpha = 0.05){
    
    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  unique_cell_type = TRUE,
                  pc_subset_ref = pc_subset)
    
    # Check if n_resamples is a positive integer
    if (!inherits(n_resamples, "numeric")) {
        stop("\'n_resamples\' should be numeric.")
    } else if (any(!n_resamples == floor(n_resamples), n_resamples < 1)) {
        stop("\'n_resamples\' should be an integer, greater than zero.")
    }
    
    # Input check for alpha
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
        stop("\'alpha\' must be a positive number greater than 0 and less than 1.")
    }

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data, 
                             reference_data = reference_data, 
                             pc_subset = pc_subset, 
                             query_cell_type_col = query_cell_type_col, 
                             ref_cell_type_col = ref_cell_type_col)
    
    # Get sample size for Wasserstein null distribution
    n_null <- min(floor(ncol(reference_data)/2), ncol(query_data), 500)
    
    # Extract variance explained
    weights <- attributes(reducedDim(reference_data, "PCA"))[["varExplained"]][pc_subset] / 
        sum(attributes(reducedDim(reference_data, "PCA"))[["varExplained"]][pc_subset])
    
    # Compute reference-reference PCA weighted distances
    pca_ref <- pca_output[pca_output$dataset == "Reference", paste0("PC", pc_subset)]
    pca_ref_weighted <- t(apply(pca_ref, 1, function(x, weights) return(x * weights), weights = sqrt(weights)))
    weighted_dist_ref <- as.matrix(dist(pca_ref_weighted))
    
    # Computing Wasserstein distances of null distribution
    null_dist <- numeric(n_resamples)
    prob_masses <- rep(1/n_null, n_null)
    for(iter in seq_len(n_resamples)){
        
        sample_ref_1 <- sample(seq_len(nrow(pca_ref)), n_null, replace = FALSE)
        sample_ref_2 <- sample(seq_len(nrow(pca_ref))[-sample_ref_1], n_null, replace = FALSE)
        cost_mat <- weighted_dist_ref[sample_ref_1, sample_ref_2]
        opt_plan <- transport::transport(prob_masses, prob_masses, costm = cost_mat)
        null_dist[iter] <- transport::wasserstein(prob_masses, prob_masses, tplan = opt_plan, costm = cost_mat)
    }
    
    # Compute reference-query PCA weighted distances
    pca_query <- pca_output[pca_output$dataset == "Query", paste0("PC", pc_subset)]
    pca_query_weighted <- t(apply(pca_query, 1, function(x, weights) return(x * weights), weights = sqrt(weights)))
    weighted_dist_query <- outer(rowSums(pca_ref_weighted^2), rowSums(pca_query_weighted^2), "+") - 
        2 * pca_ref_weighted %*% t(pca_query_weighted)

    # Computing Wasserstein distances for query data
    query_dist <- numeric(n_resamples)
    for(iter in seq_len(n_resamples)){
        
        sample_ref <- sample(seq_len(nrow(pca_ref)), n_null, replace = FALSE)
        sample_query <- sample(seq_len(nrow(pca_query)), n_null, replace = FALSE)
        cost_mat <- weighted_dist_query[sample_ref, sample_query]
        opt_plan <- transport::transport(prob_masses, prob_masses, costm = cost_mat)
        query_dist[iter] <- transport::wasserstein(prob_masses, prob_masses, tplan = opt_plan, costm = cost_mat)
    }
    
    # Visualize results
    threshold_text <- bquote(paste("Signifiance Threshold (", alpha, " = ", .(alpha), ")"))
    vline_data <- data.frame(xintercept = c(quantile(null_dist, 1 - alpha), mean(query_dist)),
                             line_type = c("Signifiance Threshold", "Reference-Query Distance"))
    density_plot <- ggplot2::ggplot(data.frame(null_dist), ggplot2::aes(x = null_dist)) +
        ggplot2::geom_density(alpha = 0.7, fill = "#00BBC4") + 
        ggplot2::labs(title = paste0("Density of Wasserstein Distances For Reference Distribution of ", 
                                     unique(reference_data[[ref_cell_type_col]])), 
                      x = "Wasserstein Distances", y = "Density") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "right", panel.grid.minor = ggplot2::element_blank(), 
                       panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                       plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10)) + 
        ggplot2::geom_vline(data = vline_data, ggplot2::aes(xintercept = .data[["xintercept"]], 
                                                            linetype = .data[["line_type"]]), 
                            color = "black", linewidth = c(1, 1)) +
        ggplot2::scale_linetype_manual(name = NULL, 
                                       values = c("Signifiance Threshold" = "solid", "Reference-Query Distance" = "dashed"),
                                       labels = c("Reference-Query Distance", threshold_text)) +
        ggplot2::guides(linetype = ggplot2::guide_legend(nrow = 2, override.aes = list(color = "black", size = 0.5), 
                                                         direction = "horizontal", 
                                                         keywidth = ggplot2::unit(1, "line"), 
                                                         keyheight = ggplot2::unit(1.5, "line")))
    return(density_plot)
}







