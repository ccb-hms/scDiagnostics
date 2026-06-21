# Compute Maximum Mean Discrepancy Statistic

Compute the Maximum Mean Discrepancy (MMD) statistic between two
datasets using either Gaussian or linear kernels for distribution
comparison.

## Usage

``` r
computeMMDStatistic(X, Y, kernel_type = "gaussian", sigma = NULL)
```

## Arguments

- X:

  A numeric matrix representing the first dataset, where rows are
  observations and columns are features.

- Y:

  A numeric matrix representing the second dataset, where rows are
  observations and columns are features.

- kernel_type:

  A character string specifying the kernel type. Options are "gaussian"
  for RBF kernel or "linear" for linear kernel. Default is "gaussian".

- sigma:

  A numeric value specifying the bandwidth parameter for the Gaussian
  kernel. If NULL, it is estimated using the median heuristic. Default
  is NULL.

## Value

A numeric value representing the MMD^2 statistic between the two
datasets.

## Details

This function calculates the MMD statistic, which measures the distance
between two probability distributions by comparing their embeddings in a
reproducing kernel Hilbert space (RKHS). For the Gaussian kernel, an
optimized median heuristic is used to estimate the bandwidth parameter
sigma when not provided. The linear kernel provides a computationally
faster alternative.

The MMD statistic is computed as: MMD^2 = E\[k(X,X')\] + E\[k(Y,Y')\] -
2\*E\[k(X,Y)\] where k is the chosen kernel function.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
