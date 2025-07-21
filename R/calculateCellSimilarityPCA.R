#' @title Calculate Cell Similarity Using PCA Loadings
#'
#' @description
#' This function calculates the cosine similarity between cells based on the principal components (PCs)
#' obtained from PCA (Principal Component Analysis) loadings.
#'
#' @details
#' This function calculates the cosine similarity between cells based on the loadings of the selected
#' principal components obtained from PCA. It extracts the rotation matrix from the PCA results of the
#' \code{\linkS4class{SingleCellExperiment}} object and identifies the high-loading variables for each selected PC.
#' Then, it computes the cosine similarity between cells using the high-loading variables for each PC.
#'
#' @param se_object A \code{\linkS4class{SingleCellExperiment}} object containing expression data.
#' @param cell_names A character vector specifying the cell names for which to compute the similarity.
#' @param pc_subset A numeric vector specifying the subset of principal components to consider. Default is 1:5..
#' @param n_top_vars An integer indicating the number of top loading variables to consider for each PC. Default is 50.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A data frame containing cosine similarity values between cells for each selected principal component.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateCellSimilarityPCAObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Store PCA anomaly data and plots
#' anomaly_output <- detectAnomaly(reference_data = reference_data,
#'                                 query_data = query_data,
#'                                 ref_cell_type_col = "expert_annotation",
#'                                 query_cell_type_col = "SingleR_annotation",
#'                                 pc_subset = 1:10,
#'                                 n_tree = 500,
#'                                 anomaly_threshold = 0.5)
#' top6_anomalies <- names(sort(anomaly_output$Combined$reference_anomaly_scores,
#'                              decreasing = TRUE)[1:6])
#'
#' # Compute cosine similarity between anomalies and top PCs
#' cosine_similarities <- calculateCellSimilarityPCA(reference_data,
#'                                                   cell_names = top6_anomalies,
#'                                                   pc_subset = 1:25,
#'                                                   n_top_vars = 50)
#' cosine_similarities
#'
#' # Plot similarities
#' plot(cosine_similarities, pc_subset = 15:25)
#'
# Function to calculate cosine similarities between cells and PCs
calculateCellSimilarityPCA <- function(se_object,
                                       cell_names,
                                       pc_subset = 1:5,
                                       n_top_vars = 50,
                                       assay_name = "logcounts",
                                       max_cells = 2500){

    # Check standard input arguments
    argumentCheck(query_data = se_object,
                  cell_names_query = cell_names,
                  pc_subset_query = pc_subset,
                  assay_name = assay_name)

    # Downsample SCE object
    se_object <- downsampleSCE(sce = se_object,
                               max_cells = max_cells)

    # Check if n_top_vars is a positive integer
    if (!is.numeric(n_top_vars) || n_top_vars <= 0 ||
        n_top_vars != as.integer(n_top_vars)) {
        stop("\'n_top_vars\' must be a positive integer.")
    }
    if(is.null(n_top_vars)){
        n_top_vars <- nrow(se_object)
    }

    # Extract rotation matrix for SingleCellExperiment object
    rotation_mat <- attributes(
        reducedDim(se_object, "PCA"))$rotation[, pc_subset]

    # Function to identify high-loading variables for each PC
    .getHighLoadingVars <- function(rotation_mat, n_top_vars) {
        high_loading_vars <- lapply(
            seq_len(ncol(rotation_mat)), function(pc) {
                abs_loadings <- abs(rotation_mat[, pc])
                top_vars <-
                    names(sort(abs_loadings,
                               decreasing = TRUE))[seq_len(n_top_vars)]
                return(top_vars)
            })
        return(high_loading_vars)
    }

    # Extract high-loading variables
    high_loading_vars <- .getHighLoadingVars(rotation_mat, n_top_vars)

    # Function to compute cosine similarity
    .cosine_similarity <- function(vector1, vector2) {
        sum(vector1 * vector2) / (sqrt(sum(vector1^2)) *
                                      sqrt(sum(vector2^2)))
    }

    # Function to compute cosine similarity for each PC using high-loading variables
    .computeCosineSimilarity <- function(cell_names,
                                         rotation_mat,
                                         high_loading_vars) {
        similarities <- lapply(
            seq_len(length(high_loading_vars)), function(pc) {
                vars <- high_loading_vars[[pc]]
                cell_subset <- cell_names[, vars, drop = FALSE]
                pc_vector <- rotation_mat[vars, pc]
                apply(cell_subset, 1, .cosine_similarity,
                      vector2 = pc_vector)
            })
        return(similarities)
    }

    # Calculate similarities
    assay_mat <- t(as.matrix(assay(se_object[, cell_names, drop = FALSE],
                                   assay_name)))
    similarities <- .computeCosineSimilarity(assay_mat, rotation_mat,
                                             high_loading_vars)

    # Format the result into a data frame for easy interpretation
    similarity_df <- do.call(cbind, similarities)
    colnames(similarity_df) <- paste0("PC", seq_len(ncol(rotation_mat)))

    # Update class of output
    class(similarity_df) <- c(class(similarity_df),
                              "calculateCellSimilarityPCAObject")
    return(similarity_df)
}
