# Calculate Entropy

This function calculates the entropy of a probability distribution.

## Usage

``` r
calculateEntropy(p)
```

## Arguments

- p:

  A numeric vector representing a probability distribution. The elements
  should sum to 1.

## Value

A numeric value representing the entropy of the probability
distribution.

## Details

The entropy is calculated using the formula \\-\sum p \log(p)\\, where
the sum is over all non-zero elements of `p`.
