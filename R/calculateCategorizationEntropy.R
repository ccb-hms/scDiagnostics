#' @title Calculate Categorization Entropy
#'
#' @description
#' This function takes a matrix of category scores (cell type by
#' cells) and calculates the entropy of the category probabilities for each
#' cell. This gives a sense of how confident the cell type assignments are.
#' High entropy = lots of plausible category assignments = low confidence. Low
#' entropy = only one or two plausible categories = high confidence. This is
#' confidence in the vernacular sense, not in the "confidence interval"
#' statistical sense. Also note that the entropy tells you nothing about
#' whether or not the assignments are correct -- see the other functionality
#' in the package for that. This functionality can be used for assessing how
#' comparatively confident different sets of assignments are (given that the
#' number of categories is the same).
#'
#' @details
#' The function checks if X is already on the probability scale.
#' Otherwise, it applies softmax columnwise.
#'
#' You can think about entropies on a scale from 0 to a maximum that depends
#' on the number of categories. This is the function for entropy (minus input
#' checking): \code{entropy(p) = -sum(p*log(p))} . If that input vector p is a
#' uniform distribution over the \code{length(p)} categories, the entropy will
#' be a high as possible.
#
#' @param X A matrix of category scores.
#' @param inverseNormalTransformationform If TRUE, apply inverse normal transformation to X. Default is FALSE.
#' @param verbose If TRUE, display messages about the calculations. Default is TRUE.
#' @param plot If TRUE, plot a histogram of the entropies. Default is TRUE.
#'
#' @returns A vector of entropy values for each column in X.
#'
#' @export
#'
#' @author Andrew Ghazi, \email{andrew_ghazi@hms.harvard.edu}
#'
#' @examples
#' # Simulate 500 cells with scores on 4 possible cell types
#' X <- rnorm(500 * 4) |> matrix(nrow = 4)
#' X[1, 1:250] <- X[1, 1:250] + 5 # Make the first category highly scored in the first 250 cells
#'
#' # The function will issue a message about softmaxing the scores, and the entropy histogram will be
#' # bimodal since we made half of the cells clearly category 1 while the other half are roughly even.
#' entropy_scores <- calculateCategorizationEntropy(X)
#'
# Function to calculate categorization entropy
calculateCategorizationEntropy <- function(X,
                                           inverseNormalTransformationform = FALSE,
                                           plot = TRUE,
                                           verbose = TRUE) {

    if (inverseNormalTransformationform) {

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

        X <- inverseNormalTransformation(X)
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

    max_ent <- calculateEntropy(rep(1 / ncat, ncat))

    if (verbose) {
        message(
            "Max possible entropy given ", ncat, " categories: ",
            round(max_ent,
                  digits = 2
            )
        )
    }

    entropies <- apply(X, 2, calculateEntropy)

    if (plot) {
        p <- data.frame(entropies = entropies) |>
            ggplot2::ggplot(ggplot2::aes(entropies)) +
            ggplot2::geom_histogram(
                color = "black", fill = "white",
                bins = 30,
                boundary = 0
            ) +
            ggplot2::theme_bw()
        p
    }

    return(entropies)
}

#' @title Calculate Entropy
#'
#' @description
#' This function calculates the entropy of a probability distribution.
#'
#' @details
#' The entropy is calculated using the formula \eqn{-\sum p \log(p)}, where the sum is over all non-zero elements of \code{p}.
#'
#' @param p A numeric vector representing a probability distribution. The elements should sum to 1.
#'
#' @keywords internal
#'
#' @return A numeric value representing the entropy of the probability distribution.
#'
# Function to calculate entropy
calculateEntropy <- function(p) {

    nonzeros <- p != 0

    -sum(p[nonzeros] * log(p[nonzeros]))
}

#' @title Number of Elements
#'
#' @description
#' This function returns the number of elements in a matrix or vector.
#'
#' @details If \code{X} is a matrix, the function returns the product of its dimensions. If \code{X} is a vector, the function returns
#' its length.
#'
#' @param X A matrix or vector.
#'
#' @keywords internal
#'
#' @return An integer representing the number of elements in \code{X}.
#'
#' @author Andrew Ghazi, \email{andrew_ghazi@hms.harvard.edu}
#'
# Function to return the number of elements
nElements <- function(X){

    return(ifelse(is.matrix(X), prod(dim(X)), length(X)))
}

#' @title Inverse Normal Transformation
#'
#' @description
#' This function performs an inverse normal transformation on a matrix or vector.
#'
#' @details
#' The function ranks the elements of \code{X} and then applies the inverse normal transformation using the formula \eqn{qnorm((rank - constant) / (n - 2 * constant + 1))}.
#'
#' @param X A numeric matrix or vector.
#' @param constant A numeric value used in the transformation. Default is \code{3 / 8}.
#'
#' @keywords internal
#'
#' @author Andrew Ghazi, \email{andrew_ghazi@hms.harvard.edu}
#'
#' @return A matrix or vector with the same dimensions as \code{X}, with values transformed using the inverse normal transformation.
#'
# Function to compute the inverse normal rank transformation
inverseNormalTransformation <- function(X,
                                        constant = 3 / 8) {

    n <- nElements(X)

    rankX <- rank(X)

    intX <- qnorm((rankX - constant) / (n - 2 * constant + 1))

    if (is.matrix(X)) {
        intX <- matrix(intX, nrow = nrow(X))
    }

    return(intX)
}
