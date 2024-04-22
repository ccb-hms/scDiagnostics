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
#' @param inverse_normal_transform if TRUE, apply
#' @param verbose if TRUE, display messages about the calculations
#' @param plot if TRUE, plot a histogram of the entropies
#' @returns A vector of entropy values for each column in X.
#' @details The function checks if X is already on the probability scale.
#'   Otherwise, it applies softmax columnwise.
#'
#'   You can think about entropies on a scale from 0 to a maximum that depends
#'   on the number of categories. This is the function for entropy (minus input
#'   checking): \code{entropy(p) = -sum(p*log(p))} . If that input vector p is a
#'   uniform distribution over the \code{length(p)} categories, the entropy will
#'   be a high as possible.
#' @export
#' @examples
#' # Simulate 500 cells with scores on 4 possible cell types
#' X <- rnorm(500 * 4) |> matrix(nrow = 4)
#' X[1, 1:250] <- X[1, 1:250] + 5 # Make the first category highly scored in the first 250 cells
#'
#' # The function will issue a message about softmaxing the scores, and the entropy histogram will be
#' # bimodal since we made half of the cells clearly category 1 while the other half are roughly even.
#' # entropy_scores <- calculateCategorizationEntropy(X)
calculateCategorizationEntropy <- function(X,
    inverse_normal_transform = FALSE,
    plot = TRUE,
    verbose = TRUE) {
    if (inverse_normal_transform) {
        # https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html#inverse-normal-transformation
        if (verbose) message("Applying global inverse normal transformation.")
        # You can't do the INT column-wise (by cell) because it will set a
        # constant "range" to the probabilities, eliminating the differences in
        # confidence across methods we're trying to quantify.

        # You can't do the INT row-wise (by cell-type) because even though
        # different cell types exhibit different marginal distributions of
        # scores (in SingleR at least), doing the transformation row-wise would
        # eliminate any differences in which cell types are "hard to predict".
        # You don't want a score of .5 for cytotoxic T cells (hard to predict
        # type) to overwhelm a score of .62 from erythroid type 2 (easy to
        # predict), even though the first would be extraordinary within its cell
        # type and the latter unexceptional within its cell type.

        X <- inverse_normal_trans(X)
    }

    colSumsX <- colSums(X)

    X_is_probabilities <- all(X >= 0 & X <= 1) &
        all((colSumsX - 1) <= 1e-8)

    if (!X_is_probabilities) {
        if (verbose) message("X doesn't seem to be on the probability scale, applying column-wise softmax.")
        expX <- exp(X)

        X <- sweep(expX, MARGIN = 2, STATS = colSums(expX), FUN = "/")
    }

    ncat <- nrow(X)

    max_ent <- calculate_entropy(rep(1 / ncat, ncat))

    if (verbose) {
        message(
            "Max possible entropy given ", ncat, " categories: ",
            round(max_ent,
                digits = 2
            )
        )
    }

    entropies <- apply(X, 2, calculate_entropy)

    if (plot) {
        p <- data.frame(entropies = entropies) |>
            ggplot(aes(entropies)) +
            geom_histogram(
                color = "black", fill = "white",
                bins = 30,
                boundary = 0
            ) +
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

n_elements <- function(X) ifelse(is.matrix(X), prod(dim(X)), length(X))

inverse_normal_trans <- function(X, constant = 3 / 8) {
    n <- n_elements(X)

    rankX <- rank(X)

    intX <- qnorm((rankX - constant) / (n - 2 * constant + 1))

    if (is.matrix(X)) {
        intX <- matrix(intX, nrow = nrow(X))
    }

    return(intX)
}
