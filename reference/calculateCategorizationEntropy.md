# Calculate Categorization Entropy

This function takes a matrix of category scores (cell type by cells) and
calculates the entropy of the category probabilities for each cell. This
gives a sense of how confident the cell type assignments are. High
entropy = lots of plausible category assignments = low confidence. Low
entropy = only one or two plausible categories = high confidence. This
is confidence in the vernacular sense, not in the "confidence interval"
statistical sense. Also note that the entropy tells you nothing about
whether or not the assignments are correct – see the other functionality
in the package for that. This functionality can be used for assessing
how comparatively confident different sets of assignments are (given
that the number of categories is the same).

## Usage

``` r
calculateCategorizationEntropy(
  X,
  inverseNormalTransformationform = FALSE,
  plot = TRUE,
  verbose = TRUE
)
```

## Arguments

- X:

  A matrix of category scores.

- inverseNormalTransformationform:

  If TRUE, apply inverse normal transformation to X. Default is FALSE.

- plot:

  If TRUE, plot a histogram of the entropies. Default is TRUE.

- verbose:

  If TRUE, display messages about the calculations. Default is TRUE.

## Value

A vector of entropy values for each column in X.

## Details

The function checks if X is already on the probability scale. Otherwise,
it applies softmax columnwise.

You can think about entropies on a scale from 0 to a maximum that
depends on the number of categories. This is the function for entropy
(minus input checking): `entropy(p) = -sum(p*log(p))` . If that input
vector p is a uniform distribution over the `length(p)` categories, the
entropy will be a high as possible.

## Author

Andrew Ghazi, <andrew_ghazi@hms.harvard.edu>

## Examples

``` r
# Simulate 500 cells with scores on 4 possible cell types
X <- rnorm(500 * 4) |> matrix(nrow = 4)
X[1, 1:250] <- X[1, 1:250] + 5 # Make the first category highly scored in the first 250 cells

# The function will issue a message about softmaxing the scores, and the entropy histogram will be
# bimodal since we made half of the cells clearly category 1 while the other half are roughly even.
entropy_scores <- calculateCategorizationEntropy(X)
#> X doesn't seem to be on the probability scale, applying column-wise softmax.
#> Max possible entropy given 4 categories: 1.39
```
