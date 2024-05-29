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
#' @examples
#' \donttest{
#' # Load necessary library
#' library(scRNAseq)
#' library(scuttle)
#' library(scran)
#' library(SingleR)
#' library(ComplexHeatmap)
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
#' # Log transform datasets
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
#' # Subset reference and query data for a specific cell type
#' ref_data_subset <- ref_data_subset[, which(ref_data_subset$reclustered.broad == "CD8")]
#' query_data_subset <- query_data_subset[, which(colData(query_data_subset)$labels == "CD8")]
#'
#' # Run PCA on the reference and query datasets
#' ref_data_subset <- runPCA(ref_data_subset)
#' query_data_subset <- runPCA(query_data_subset)
#'
#' # Call the PCA comparison function
#' similarity_mat <- comparePCA(query_data_subset, ref_data_subset, 
#'                              n_components = 5, 
#'                              metric = c("cosine", "correlation")[1], 
#'                              correlation_method = c("spearman", "pearson")[1])
#'
#' # Create the heatmap
#' plot(similarity_mat)
#' }
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
        ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                      midpoint = 0, limit = c(min(x, -0.5), max(x, 0.5)), space = "Lab", 
                                      name = "Cosine Similarity") +
        ggplot2::theme_minimal() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, 
                                                           size = 12, hjust = 1)) +
        ggplot2::labs(x = "", y = "", 
                      title = "Heatmap of Cosine Similarities Between PCs")
    print(pc_plot)
}
