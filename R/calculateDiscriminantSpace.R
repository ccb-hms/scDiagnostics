#' @title Project Query Data onto a Unified Discriminant Space of Reference Data
#'
#' @description
#' This function projects query single-cell RNA-seq data onto a unified discriminant space defined by reference data. The reference data
#' is used to identify important variables across all cell types and compute discriminant vectors, which are then used to project both reference and query
#' data. Similarity between the query and reference projections can be assessed using cosine similarity and Mahalanobis distance.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Identifies the top important variables to distinguish cell types from the reference data by taking the union of important variables from pairwise comparisons.
#'   \item Computes the Ledoit-Wolf shrinkage estimate of the covariance matrix for each cell type using these important genes.
#'   \item Constructs within-class and between-class scatter matrices.
#'   \item Solves the generalized eigenvalue problem to obtain discriminant vectors.
#'   \item Projects both reference and query data onto the unified discriminant space.
#'   \item Assesses similarity of the query data projection to the reference data using cosine similarity and Mahalanobis distance.
#' }
#'
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' If NULL, only the projected reference data is returned. Default is NULL.
#' @param ref_cell_type_col The column name in \code{reference_data} indicating cell type labels.
#' @param query_cell_type_col The column name in \code{query_data} indicating cell type labels.
#' @param cell_types A character vector specifying the cell types to include in the analysis. If NULL, all cell types are included.
#' @param n_tree An integer specifying the number of trees for the random forest used in variable importance calculation.
#' @param n_top An integer specifying the number of top variables to select based on importance scores from each pairwise comparison.
#' @param eigen_threshold A numeric value specifying the threshold for retaining eigenvalues in discriminant analysis.
#' @param calculate_metrics Parameter to determine if cosine similarity and Mahalanobis distance metrics should be computed. Default is FALSE.
#' @param alpha A numeric value specifying the significance level for Mahalanobis distance cutoff.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A list with the following components:
#' \item{discriminant_eigenvalues}{Eigenvalues from the discriminant analysis.}
#' \item{discriminant_eigenvectors}{Eigenvectors from the discriminant analysis.}
#' \item{ref_proj}{Reference data projected onto the discriminant space.}
#' \item{query_proj}{Query data projected onto the discriminant space (if query_data is provided).}
#' \item{query_mahalanobis_dist}{Mahalanobis distances of query projections (if calculate_metrics is TRUE).}
#' \item{mahalanobis_crit}{Cutoff value for Mahalanobis distance significance (if calculate_metrics is TRUE).}
#' \item{query_cosine_similarity}{Cosine similarity scores of query projections (if calculate_metrics is TRUE).}
#'
#' @references
#' \itemize{
#' \item Fisher, R. A. (1936). "The Use of Multiple Measurements in Taxonomic Problems". *Annals of Eugenics*. 7 (2): 179–188. doi:10.1111/j.1469-1809.1936.tb02137.x.
#' \item Hastie, T., Tibshirani, R., & Friedman, J. (2009). *The Elements of Statistical Learning: Data Mining, Inference, and Prediction*. Springer. Chapter 4: Linear Methods for Classification.
#' \item Ledoit, O., & Wolf, M. (2004). "A well-conditioned estimator for large-dimensional covariance matrices". *Journal of Multivariate Analysis*. 88 (2): 365–411. doi:10.1016/S0047-259X(03)00096-4.
#' \item De Maesschalck, R., Jouan-Rimbaud, D., & Massart, D. L. (2000). "The Mahalanobis distance". *Chemometrics and Intelligent Laboratory Systems*. 50 (1): 1–18. doi:10.1016/S0169-7439(99)00047-7.
#' \item Breiman, L. (2001). "Random Forests". *Machine Learning*. 45 (1): 5–32. doi:10.1023/A:1010933404324.
#' }
#'
#' @export
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateDiscriminantSpaceObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Compute discriminant space using unified model across all cell types
#' disc_output <- calculateDiscriminantSpace(reference_data = reference_data,
#'                                           query_data = query_data,
#'                                           query_cell_type_col = "SingleR_annotation",
#'                                           ref_cell_type_col = "expert_annotation",
#'                                           n_tree = 500,
#'                                           n_top = 50,
#'                                           eigen_threshold  = 1e-1,
#'                                           calculate_metrics = FALSE,
#'                                           alpha = 0.01)
#'
#' # Generate scatter and boxplot
#' plot(disc_output, plot_type = "scatterplot")
#' plot(disc_output, cell_types = c("CD4", "CD8"), plot_type = "boxplot")
#'
#' # Check comparison
#' table(Expert_Annotation = query_data$expert_annotation,
#'       SingleR = query_data$SingleR_annotation)
#'
#' @importFrom stats cov qchisq mahalanobis
#'
# Function to get the discriminant spaces and the projected reference/query data on the discriminant space.
# Similarity measures (cosine similarity/Mahalanobis distance) for the projected query data are also available.
calculateDiscriminantSpace <- function(reference_data,
                                       query_data = NULL,
                                       ref_cell_type_col,
                                       query_cell_type_col = NULL,
                                       cell_types = NULL,
                                       n_tree = 500,
                                       n_top = 20,
                                       eigen_threshold  = 1e-1,
                                       calculate_metrics = FALSE,
                                       alpha = 0.01,
                                       assay_name = "logcounts",
                                       max_cells = 2500){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  assay_name = assay_name)

    # Downsample reference and query data
    reference_data <- downsampleSCE(sce = reference_data,
                                    max_cells = max_cells)
    if(!is.null(query_data)){
        query_data <- downsampleSCE(sce = query_data,
                                    max_cells = max_cells)
    }

    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        if(is.null(query_data)){
            cell_types <- na.omit(
                unique(c(reference_data[[ref_cell_type_col]])))
        } else{
            cell_types <- na.omit(
                unique(c(reference_data[[ref_cell_type_col]],
                         query_data[[query_cell_type_col]])))
        }
    }

    # Check if n_tree is a positive integer
    if (!is.numeric(n_tree) || n_tree <= 0 || n_tree != as.integer(n_tree)) {
        stop("\'n_tree\' must be a positive integer.")
    }

    # Check if n_top is a positive integer
    if (!is.numeric(n_top) || n_top <= 0 || n_top != as.integer(n_top)) {
        stop("\'n_top\' must be a positive integer.")
    }

    # Input check for eigen_threshold
    if (!is.numeric(eigen_threshold) || eigen_threshold <= 0) {
        stop("\'eigen_threshold\' must be a positive number greater than 0.")
    }

    # Input check for alpha
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
        stop("\'alpha\' must be a positive number greater than 0 and less than 1.")
    }

    # Getting top variables
    var_imp <- calculateVarImpOverlap(reference_data = reference_data,
                                      ref_cell_type_col = ref_cell_type_col,
                                      n_tree = n_tree,
                                      n_top = n_top)

    # Get union of top n_top genes from all pairwise combinations
    all_top_genes <- unique(unlist(
        lapply(var_imp[["var_imp_ref"]], function(x) {
            x[["Gene"]][seq_len(n_top)]
        })))

    # Create a single discriminant model using all_top_genes
    # Extract reference matrix using top genes
    ref_mat <- t(as.matrix(assay(reference_data, assay_name)))[, all_top_genes]

    # Compute within-class and between-class scatter matrix
    sw <- sb <- matrix(0, length(all_top_genes), length(all_top_genes))
    overall_mean <- colMeans(ref_mat)

    for(cell_type in cell_types){
        # Extract matrix for cell type
        ref_mat_subset <- ref_mat[which(
            reference_data[[ref_cell_type_col]] == cell_type),]

        # Ledoit-Wolf estimation for the current class
        lw_cov <- ledoitWolf(ref_mat_subset)
        # Update within-class scatter matrix
        sw <- sw + (nrow(ref_mat_subset) - 1) * lw_cov

        # Update between-class scatter matrix
        class_mean <- colMeans(ref_mat_subset)
        sb <- sb + nrow(ref_mat_subset) * (class_mean - overall_mean) %*%
            t(class_mean - overall_mean)
    }

    # Ensure symmetry of sw
    sw <- (sw + t(sw)) / 2

    # Solve generalized eigenvalue problem
    eig <- eigen(solve(sw) %*% sb)

    # Sort eigenvectors by eigenvalues
    discriminant_eigenvalues <- Re(
        eig[["values"]][which(Re(eig[["values"]]) > eigen_threshold),
                        drop = FALSE])
    discriminant_eigenvectors <- Re(
        eig[["vectors"]][, which(Re(eig[["values"]]) > eigen_threshold),
                         drop = FALSE])
    rownames(discriminant_eigenvectors) <- all_top_genes
    colnames(discriminant_eigenvectors) <- paste0(
        "DV", seq_len(ncol(discriminant_eigenvectors)))

    # Compute projected data for reference
    ref_proj <- data.frame(
        ref_mat[which(reference_data[[ref_cell_type_col]] %in% cell_types),] %*%
            discriminant_eigenvectors,
        reference_data[[ref_cell_type_col]][reference_data[[ref_cell_type_col]] %in%
                                                cell_types])
    colnames(ref_proj) <- c(
        paste0("DV", seq_len(ncol(discriminant_eigenvectors))),
        "cell_type")

    # Create a single entry in discriminant_output
    discriminant_output <- list()
    discriminant_output <- list()
    discriminant_output[["discriminant_eigenvalues"]] <- discriminant_eigenvalues
    discriminant_output[["discriminant_eigenvectors"]] <- discriminant_eigenvectors
    discriminant_output[["ref_proj"]] <- ref_proj

    # Computations for query data
    if(!is.null(query_data)){

        # Projection on discriminant space
        query_mat <- t(as.matrix(assay(query_data, assay_name)))[, all_top_genes]
        query_proj <- data.frame(
            query_mat[which(query_data[[query_cell_type_col]] %in%
                                cell_types),] %*% discriminant_eigenvectors,
            query_data[[query_cell_type_col]][query_data[[query_cell_type_col]] %in%
                                                  cell_types])
        colnames(query_proj) <- c(
            paste0("DV", seq_len(ncol(discriminant_eigenvectors))), "cell_type")
        discriminant_output[["query_proj"]] <- query_proj

        if(calculate_metrics){

            # Cosine similarity between mean vector of reference projection
            # Mahalanobis distance between each query cell projected on reference discriminant space
            # and reference data projected on reference discriminant space
            cosine_similarity <- mahalanobis_dist <- numeric(nrow(query_proj))
            mahalanobis_crit <- numeric(length(cell_types))
            for(type_idx in seq_along(cell_types)){
                type <- cell_types[type_idx]
                query_cells_of_type <- query_proj[
                    query_proj[, "cell_type"] == type,
                    paste0("DV",
                           seq_len(length(discriminant_eigenvalues)))]
                ref_cells_of_type <- ref_proj[
                    ref_proj[, "cell_type"] == type,
                    paste0("DV",
                           seq_len(length(discriminant_eigenvalues)))]

                # Skip if we have no cells of this type
                if(nrow(query_cells_of_type) == 0 || nrow(ref_cells_of_type) == 0) {
                    next
                }

                # Calculate mean of reference projection for this cell type
                ref_mean <- colMeans(ref_cells_of_type)

                # Calculate covariance of reference projection for this cell type
                ref_cov <- cov(ref_cells_of_type)

                # Check if covariance matrix is invertible
                if(any(is.na(ref_cov)) || determinant(ref_cov)[["modulus"]][1] <= 0) {
                    # If not invertible, use a regularized version
                    ref_cov <- ledoitWolf(ref_cells_of_type)
                }

                # Calculate Mahalanobis distance
                mahalanobis_dist[query_proj[, "cell_type"] == type] <- mahalanobis(
                    query_cells_of_type, ref_mean, ref_cov)

                # Calculate cosine similarity
                cosine_similarity[query_proj[, "cell_type"] == type] <-
                    apply(query_cells_of_type, 1,
                          function(x, y) return(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))),
                          y = ref_mean)
            }
            discriminant_output[["query_mahalanobis_dist"]] <-
                mahalanobis_dist
            discriminant_output[["mahalanobis_crit"]] <-
                qchisq(1 - alpha, df = length(cell_types))
            discriminant_output[["query_cosine_similarity"]] <-
                cosine_similarity
        }
    }

    # Return data projected onto (reference) discriminant space for single combined model
    class(discriminant_output) <- c(class(discriminant_output),
                                    "calculateDiscriminantSpaceObject")
    return(discriminant_output)
}

#' @title Ledoit-Wolf Covariance Matrix Estimation
#'
#' @description
#' Estimate the covariance matrix using the Ledoit-Wolf shrinkage method.
#'
#' @details
#' This function computes the Ledoit-Wolf shrinkage covariance matrix estimator,
#' which improves the accuracy of the sample covariance matrix by shrinking it
#' towards a structured estimator, typically the diagonal matrix with the mean
#' of variances as its diagonal elements.
#'
#' @param class_data A numeric matrix or data frame containing the data for covariance estimation,
#' where rows represent observations and columns represent variables.
#'
#' @keywords internal
#'
#' @return A numeric matrix representing the Ledoit-Wolf estimated covariance matrix.
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to compute Ledoit-Wolf covariance matrix
ledoitWolf <- function(class_data) {

    # Sample covariance matrix
    sample_cov <- cov(class_data)

    # Check for zero column means and replace them with a small constant if necessary
    col_means <- colMeans(class_data)
    col_means[col_means == 0] <- 1e-10

    # Calculate the shrinkage target (identity matrix scaled by the average variance)
    mean_variance <- mean(diag(sample_cov))
    shrinkage_target <- diag(mean_variance, ncol(class_data),
                             ncol(class_data))

    # Calculate the shrinkage intensity
    phi_hat <- sum((class_data - col_means) ^ 2) /
        (nrow(class_data) - 1)
    shrinkage_intensity <- (1 / nrow(class_data)) *
        min(phi_hat, mean_variance ^ 2)

    # Ledoit-Wolf estimated covariance matrix
    lw_cov <- (1 - shrinkage_intensity) *
        sample_cov + shrinkage_intensity * shrinkage_target

    # Ensure symmetry
    lw_cov <- (lw_cov + t(lw_cov)) / 2

    # Return Ledoit-Wolf estimated covariance matrix
    return(lw_cov)
}
