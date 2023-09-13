#' Calculate Categorization Entropy
#' @description This function takes a matrix of category scores (cell type by
#'   cells) and calculates the entropy of the category probabilities for each
#'   cell. This gives a sense of how confident the cell type assignments are.
#'   High entropy = lots of plausible category assignments = low confidence. Low
#'   entropy = only one or two plausible categories = high confidence. This is
#'   confidence in the vernacular sense, not in the "confidence interval"
#'   statistical sense. Also note that the entropy tells you nothing about
#'   whether or not the assignments are correct -- see the other functionality
#'   in the package for that. This functionality can be used for assessing how
#'   comparatively confident different sets of assignments are (given that the
#'   number of categories is the same).
#' @param X a matrix of category scores
#' @param plot if TRUE, plot a histogram of the entropies
#' @returns A vector of entropy values for each column in X.
#' @details The function checks if X is already on the probability scale.
#'   Otherwise, it divides out the sum columnwise.
#'
#'   You can think about entropies on a scale from 0 to a maximum that depends
#'   on the number of categories. This is the function for entropy (minus input
#'   checking): \code{entropy(p) = -sum(p*log(p))} . If that input vector p is a
#'   uniform distribution over the \code{length(p)} categories, the entropy will
#'   be a high as possible.
#' @export
#' @examples
#' # Prepare the type assigment matrix X using SingleR like in the readme. Skip to the bottom two commands if you already have your type assignment matrix.
#' library(scRNAseq)
#' library(scuttle)
#' library(SingleR)
#' library(scater)
#'
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' # Log-transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Run PCA
#' ref_data <- runPCA(ref_data)
#' query_data <- runPCA(query_data)
#'
#' # Get cell type scores using SingleR
#' pred <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' pred <- as.data.frame(pred)
#'
#' X <- pred[, grepl("scores", names(pred))] |>
#'     as.matrix() |>
#'     t()
#' entropy_scores <- calculateCategorizationEntropy(X)
calculateCategorizationEntropy <- function(X, plot = TRUE, verbose = TRUE) {
    if (any(X < 0)) {
        stop("All entries in X must be >= 0.")
    }

    colSumsX <- colSums(X)

    X_is_probabilities <- all(X >= 0 & X <= 1) &
        all((colSumsX - 1) <= 1e-8)

    if (!X_is_probabilities) {
        if (verbose) message("X doesn't seem to be on the probability scale, dividing out the colSums.")
        X <- sweep(X, MARGIN = 2, STATS = colSumsX, FUN = "/")
    }

    if (verbose) {
        ncat <- nrow(X)
        
        message(
            "Max possible entropy given ", ncat, " categories: ",
            round(calculate_entropy(rep(1 / ncat, ncat)),
                  digits = 2
            )
        )
    }

    entropies <- apply(X, 2, calculate_entropy)

    if (plot) {
        p <- data.frame(entropies = entropies) |>
            ggplot(aes(entropies)) +
            geom_histogram(color = "black", fill = "white",
                           bins = 30) +
            theme_bw()
        print(p)
    }

    return(entropies)
}

calculate_entropy <- function(p) {
    # p is one column of X, a vector of probabilities summing to 1.

    nonzeros <- p != 0

    -sum(p[nonzeros] * log(p[nonzeros]))
}
