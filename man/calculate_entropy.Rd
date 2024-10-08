% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateCategorizationEntropy.R
\name{calculate_entropy}
\alias{calculate_entropy}
\title{Calculate Entropy}
\usage{
calculate_entropy(p)
}
\arguments{
\item{p}{A numeric vector representing a probability distribution. The elements should sum to 1.}
}
\value{
A numeric value representing the entropy of the probability distribution.
}
\description{
This function calculates the entropy of a probability distribution.
}
\details{
The entropy is calculated using the formula \eqn{-\sum p \log(p)}, where the sum is over all non-zero elements of \code{p}.
}
\keyword{internal}
