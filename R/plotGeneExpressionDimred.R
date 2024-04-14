#' Visualize gene expression on a dimensional reduction plot
#'
#' This function plots gene expression on a dimensional reduction plot
#' using methods like t-SNE, UMAP, or PCA. Each single cell is
#' color-coded based on the expression of a specific gene or feature.
#'
#' @param se_object An object of class "SingleCellExperiment"
#'     containing log-transformed expression matrix and other
#'     metadata.  It can be either a reference or query dataset.
#'
#' @param method The reduction method to use for visualization. It
#'     should be one of the supported methods: "tSNE", "UMAP", or
#'     "PCA".
#'
#' @param n_components A numeric vector of length 2 indicating the
#'     first two dimensions to be used for plotting.
#'
#' @param feature A character string representing the name of the gene
#'     or feature to be visualized.
#'
#' @return A ggplot object representing the dimensional reduction plot
#'     with gene expression.
#'
#' @examples
#'
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#'
#' # Load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(
#'     ncol(assay(sce)),
#'     size = floor(0.7 * ncol(assay(sce))),
#'     replace = FALSE
#' )
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' # Log transform datasets
#' query_data <- logNormCounts(query_data)
#'
#' # Run PCA
#' query_data <- runPCA(query_data)
#'
#' # Plot gene expression on PCA plot
#' plotGeneExpressionDimred(se_object = query_data,
#'                          method = "PCA",
#'                          n_components = c(1, 2),
#'                          feature = "VPREB3")
#'
#'
#' @importFrom ggplot2 ggplot xlab ylab scale_color_gradient theme_bw
#' @importFrom rlang .data
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
#'
#' @export
#'
plotGeneExpressionDimred <- 
    function(se_object,
             method,
             n_components = c(1, 2),
             feature) 
{
    ## Error handling and validation
    supported_methods <- c("tSNE", "UMAP", "PCA")
    if (!(method %in% supported_methods)) {
        stop("Unsupported method. Please choose one of: ",
             paste(supported_methods, collapse = ", "))
    }

    if (length(n_components) != 2) {
        stop("n_components should be a numeric vector of length 2.")
    }

    if (!feature %in% rownames(assay(se_object, "logcounts"))) {
        stop("Specified feature does not exist in the expression matrix.")
    }

    ## Extract dimension reduction coordinates from SingleCellExperiment object
    reduction <- reducedDim(se_object, method)[, n_components]

    ## Extract gene expression vector
    expression <- assay(se_object, "logcounts")[feature, ]

    ## Prepare data for plotting
    df <- data.frame(Dim1 = reduction[, 1],
                     Dim2 = reduction[, 2],
                     Expression = expression)

    ## Create the plot object
    plot <- ggplot(df, aes(x = .data$Dim1, y = .data$Dim2)) +
        geom_point(aes(color = .data$Expression)) +
        scale_color_gradient(low = "grey90", high = "blue") +
        xlab("Dimension 1") +
        ylab("Dimension 2") +
        theme_bw()

    plot
}
