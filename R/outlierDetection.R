#' Calculate cell outlier scores using isolation forests
#'
#' @description This function calculates outlier scores for each cell using
#'   \link[isotree]{isolation.forest}. Isolation forests explicitly try to
#'   estimate the isolation of datapoints, in contrast to estimating the full
#'   distribution of the data. To quote the isotree documentation: "Isolation
#'   Forest is an algorithm originally developed for outlier detection that
#'   consists in splitting sub-samples of the data according to some
#'   attribute/feature/column at random." Basically, it calculates how quickly
#'   each cell gets split off from the rest of the group when the dataset is
#'   randomly partitioned by a hyperplane. There are many adjustable parameters
#'   (e.g. number of trees, dimensionality of the splits, etc), but the defaults
#'   usually work pretty well.
#' @param sce a SingleCellExperiment object containing a logcounts matrix
#' @param group.var if specified, run the outlier scoring by group (e.g. cell
#'   cluster) instead of across all cells globally.
#' @param use_pcs if TRUE, run the outlier detection in principal component
#'   space. Otherwise, use logcounts.
#' @param prediction_thresh threshold to use for binary outlier calling
#' @param ... arguments to pass to \link[isotree]{isolation.forest}
#' @details If \code{use_pcs} is turned on, it will use all available PCs in the
#'   PCA reduced dimension slot if it exists, and run \link[scater]{runPCA} with
#'   default settings and use that if it does not.
#'
#'   \code{prediction_thresh} is the threshold to use on the isoforest score
#'   predictions
#' @returns the SingleCellExperiment object with \code{outlier_score} and
#'   \code{is_outlier} added to the colData. Also, the resulting isotree model
#'   is attached to the object's metadata.
#'
#' @examples
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#'
#' # Load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' query_data <- logNormCounts(sce)
#'
#' query_data <- runPCA(query_data)
#'
#' query_data <- calculateOutlierScore(query_data, use_pcs = TRUE)
#' @export
calculateOutlierScore <- function(sce,
                                  group.var = NULL,
                                  use_pcs = FALSE,
                                  prediction_thresh = 0.5,
                                  ...) {
  if (!is.null(group.var)) {
    if (is.null(colnames(sce))) stop("This function requires the input SCE to have unique column names when calculating group-wise outlier scores.")
    
    # There's probably a better way to split an sce by a factor...
    # Anyway, manually passing it through factor() like this way retains NAs.
    cluster_list <- lapply(
      split(
        colData(sce),
        factor(colData(sce)[[group.var]],
               exclude = NULL
        )
      ),
      \(x) sce[, rownames(x)]
    )
    
    # TODO: assess if it would be worth making the threshold adaptive. Different
    # clusters can have different distributions of outlier scores.
    score_list <- lapply(cluster_list,
                         calculateOutlierScore,
                         group.var = NULL,
                         use_pcs = use_pcs,
                         prediction_thresh = prediction_thresh
    )
    
    score_df <- Reduce(
      lapply(
        score_list,
        \(x) colData(x)[, c("outlier_score", "is_outlier")]
      ),
      f = rbind
    )
    
    # If we were using an ID column instead of column names of the input we could
    # do a join.
    score_df <- score_df[rownames(colData(sce)), ]
    colData(sce) <- cbind(colData(sce), score_df)
    
    S4Vectors::metadata(sce)$cluster_isotree_model <- lapply(score_list, S4Vectors::metadata)
    
    return(sce)
  }
  
  if (use_pcs) {
    if (!("PCA" %in% SingleCellExperiment::reducedDimNames(sce))) {
      sce <- scater::runPCA(sce)
    }
    
    X <- sce |>
      SingleCellExperiment::reducedDim("PCA")
  } else {
    X <- sce |>
      SingleCellExperiment::logcounts() |>
      BiocGenerics::t()
  }
  
  isotree_res <- X |>
    isotree::isolation.forest(
      output_score = TRUE,
      ...
    )
  
  # V These scores are the same as running predict(isotree_res$model, X)
  sce$outlier_score <- isotree_res$scores
  sce$is_outlier <- isotree_res$scores >= prediction_thresh
  
  S4Vectors::metadata(sce)$isotree_model <- isotree_res$model
  
  return(sce)
}

