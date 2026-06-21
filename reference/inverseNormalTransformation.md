# Inverse Normal Transformation

This function performs an inverse normal transformation on a matrix or
vector.

## Usage

``` r
inverseNormalTransformation(X, constant = 3/8)
```

## Arguments

- X:

  A numeric matrix or vector.

- constant:

  A numeric value used in the transformation. Default is `3 / 8`.

## Value

A matrix or vector with the same dimensions as `X`, with values
transformed using the inverse normal transformation.

## Details

The function ranks the elements of `X` and then applies the inverse
normal transformation using the formula \\qnorm((rank - constant) / (n -
2 \* constant + 1))\\.

## Author

Andrew Ghazi, <andrew_ghazi@hms.harvard.edu>
