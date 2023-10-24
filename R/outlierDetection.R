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
#' @param dimred if specified, plot the specified reduced dimension with
#'   \link[scater]{plotReducedDim}, colored by outlier score
#' @param use_pcs if TRUE, run the outlier detection in principal component
#'   space
#' @param prediction_thresh threshold to use for outlier calling
#' @param ... arguments to pass to \link[isotree]{isolation.forest}
#' @details If \code{use_pcs} is turned on, it will use all available PCs in the
#'   PCA reduced dimension slot if it exists, and run \link[scater]{runPCA} with
#'   default settings and use that if it does not.
#'
#'   \code{prediction_thresh} is the threshold to use on the isoforest score
#'   predictions
#' @returns the SingleCellExperiment object with \code{outlier_score} and
#'   \code{is_outlier} added to the colData.
#'
#' @export
calculateOutlierScore <- function(sce,
    dimred = NULL,
    use_pcs = FALSE,
    prediction_thresh = 0.5,
    ...) {
    if (use_pcs) {
        if (!("PCA" %in% SingleCellExperiment::reducedDimNames(sce))) {
            sce <- scater::runPCA(sce)
        }

        X <- sce |>
            reducedDim("PCA")

        isotree_res <- X |>
            isotree::isolation.forest(
                output_score = TRUE,
                ...
            )
    } else {
        X <- sce |>
            logcounts() |>
            t()

        isotree_res <- X |>
            isotree::isolation.forest(
                output_score = TRUE,
                ...
            )
    }

    # V These scores are the same as running predict(isotree_res$model, X)
    colData(sce)$outlier_score <- isotree_res$scores
    colData(sce)$is_outlier <- isotree_res$scores >= prediction_thresh

    if (!is.null(dimred)) {
        scater::plotReducedDim(sce, dimred, color_by = "outlier_score")
    }

    return(sce)
}
