# Ledoit-Wolf Covariance Matrix Estimation

Estimate the covariance matrix using the Ledoit-Wolf shrinkage method.

## Usage

``` r
ledoitWolf(class_data)
```

## Arguments

- class_data:

  A numeric matrix or data frame containing the data for covariance
  estimation, where rows represent observations and columns represent
  variables.

## Value

A numeric matrix representing the Ledoit-Wolf estimated covariance
matrix.

## Details

This function computes the Ledoit-Wolf shrinkage covariance matrix
estimator, which improves the accuracy of the sample covariance matrix
by shrinking it towards a structured estimator, typically the diagonal
matrix with the mean of variances as its diagonal elements.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
