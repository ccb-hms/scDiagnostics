#' @title Project Query Data Onto PCA Space of Reference Data
#'
#' @description 
#' This function projects a query singleCellExperiment object onto the PCA space of a reference 
#' singleCellExperiment object. The PCA analysis on the reference data is assumed to be pre-computed and stored within the object.
#'
#' @details 
#' This function assumes that the "PCA" element exists within the \code{reducedDims} of the reference data 
#' (obtained using \code{reducedDim(reference_data)}) and that the genes used for PCA are present in both the reference and query data. 
#' It performs centering and scaling of the query data based on the reference data before projection.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param n_components An integer specifying the number of principal components to use for projection. Defaults to 10. 
#' Must be less than or equal to the number of components available in the reference PCA.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} 
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} 
#' that identifies the cell types.
#' @param return_value A character string specifying the format of the returned data. Can be \code{data.frame} (combined reference 
#' and query projections) or \code{list} (separate lists for reference and query projections) (default = \code{data.frame}).
#'
#' @return A \code{data.frame} containing the projected data in rows (reference and query data combined) or a \code{list} containing 
#' separate matrices for reference and query projections, depending on the \code{return_value} argument.
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
#' library(RColorBrewer)
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
#' # Project the query data onto PCA space of reference
#' pca_output <- projectPCA(query_data_subset, ref_data_subset,
#'                          n_components = 10,
#'                          query_cell_type_col = "labels",
#'                          ref_cell_type_col = "reclustered.broad",
#'                          return_value = c("data.frame", "list")[1])
#'
#' # Compute t-SNE and UMAP using first 10 PCs
#' tsne_data <- data.frame(calculateTSNE(t(pca_output[, paste0("PC", 1:10)])))
#' umap_data <- data.frame(calculateUMAP(t(pca_output[, paste0("PC", 1:10)])))
#'
#' # Combine the cell type labels from both datasets
#' tsne_data$Type <- paste(pca_output$dataset, pca_output$cell_type)
#'
#' # Define the cell types and legend order
#' legend_order <- c("Query CD8",
#'                   "Reference CD8",
#'                   "Query CD4",
#'                   "Reference CD4",
#'                   "Query B_and_plasma",
#'                   "Reference B_and_plasma")
#'
#' # Define the colors for cell types
#' color_palette <- brewer.pal(length(legend_order), "Paired")
#' color_mapping <- setNames(color_palette, legend_order)
#' cell_type_colors <- color_mapping[legend_order]
#'
#' # Visualize t-SNE output
#' tsne_plot <- ggplot(tsne_data[tsne_data$Type %in% legend_order,],
#'                     aes(x = TSNE1, y = TSNE2, color = factor(Type, levels = legend_order))) +
#'     geom_point(alpha = 0.5, size = 1) +
#'     scale_color_manual(values = cell_type_colors) +
#'     theme_bw() +
#'     guides(color = guide_legend(title = "Cell Types"))
#' tsne_plot     
#' 
#'
# Function to project query data onto PCA space of reference data
projectPCA <- function(query_data, reference_data, 
                       n_components = 10, 
                       query_cell_type_col = NULL, 
                       ref_cell_type_col = NULL, 
                       return_value = c("data.frame", "list")[1]){
    
    # Check if query_data is a SingleCellExperiment object
    if (!is(query_data, "SingleCellExperiment")) {
        stop("query_data must be a SingleCellExperiment object.")
    }
    
    # Check if reference_data is a SingleCellExperiment object
    if (!is(reference_data, "SingleCellExperiment")) {
        stop("reference_data must be a SingleCellExperiment object.")
    }
    
    # Check if "PCA" is present in reference's reduced dimensions
    if (!"PCA" %in% names(reducedDims(reference_data))) {
        stop("Reference data must have pre-computed PCA in \'reducedDims\'.")
    }
    
    # Check if n_components is a positive integer
    if (!inherits(n_components, "numeric")) {
        stop("n_components should be numeric")
    } else if (any(!n_components == floor(n_components), n_components < 1)) {
        stop("n_components should be an integer, greater than zero.")
    }
    
    # Check if requested number of components is within available components
    if (ncol(reducedDim(reference_data, "PCA")) < n_components) {
        stop("\'n_components\' is larger than number of available components in reference PCA.")
    }
    
    # Returning output as single matrix or a list
    if (!return_value %in% c("data.frame", "list")) {
        stop("Invalid \'return_value\'. Must be 'data.frame' or \'list\'.")
    }
    
    # Extract reference PCA components and rotation matrix
    ref_mat <- reducedDim(reference_data, "PCA")[, 1:n_components]
    rotation_mat <- attributes(reducedDim(reference_data, "PCA"))$rotation[, 1:n_components]
    PCA_genes <- rownames(rotation_mat)
    
    # Check if genes used for PCA are available in query data
    if (!all(PCA_genes %in% rownames(assay(query_data)))) {
        stop("Genes in reference PCA are not found in query data.")
    }
    
    # Center and scale query data based on reference for projection
    centering_vec <- apply(t(as.matrix(assay(reference_data, "logcounts"))), 2, mean)[PCA_genes]
    query_mat <- scale(t(as.matrix(assay(query_data, "logcounts")))[, PCA_genes], center = centering_vec, scale = FALSE) %*% 
        rotation_mat
    
    # Returning output as single matrix or a list
    if (return_value == "data.frame") {
        return(data.frame(rbind(ref_mat, query_mat), 
                          dataset = c(rep("Reference", nrow(ref_mat)), rep("Query", nrow(query_mat))),
                          cell_type = c(ifelse(rep(is.null(ref_cell_type_col), nrow(ref_mat)), 
                                               rep(NA, nrow(ref_mat)), 
                                               colData(reference_data)[[ref_cell_type_col]]),
                                        ifelse(rep(is.null(query_cell_type_col), nrow(query_mat)), 
                                               rep(NA, nrow(query_mat)), 
                                               colData(query_data)[[query_cell_type_col]]))))
    } else if (return_value == "list") {
        return(list(ref = data.frame(ref_mat, 
                                     cell_type = ifelse(rep(is.null(ref_cell_type_col), nrow(ref_mat)), 
                                                        rep(NA, nrow(ref_mat)), 
                                                        colData(reference_data)[[ref_cell_type_col]])), 
                    query = data.frame(query_mat,
                                       cell_type = ifelse(rep(is.null(query_cell_type_col), nrow(query_mat)), 
                                                          rep(NA, nrow(query_mat)), 
                                                          colData(query_data)[[query_cell_type_col]]))))
    }
}