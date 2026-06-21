# Fast Custom Linear Regression for Principal Components

Performs efficient linear regression of principal component scores
against categorical predictors (cell types, batches, datasets, or their
interactions) using QR decomposition for numerical stability and
computational efficiency.

## Usage

``` r
regressFastCustom(pc, indep_var, df)
```

## Arguments

- pc:

  A character string specifying the principal component column name in
  the data frame.

- indep_var:

  A character string specifying the independent variable specification.
  Options include "cell_type", "cell_type \* batch", "cell_type \*
  dataset", or other interaction specifications.

- df:

  A data frame containing the principal component scores and categorical
  predictors. Must include columns for the specified PC and predictor
  variables.

## Value

A list containing:

- coefficients:

  A data frame with columns `coef`, `se`, `t`, and `p.value` containing
  regression coefficients, standard errors, t-statistics, and p-values.

- r_squared:

  A numeric value representing the R-squared of the model.

## Details

This function implements a custom linear regression optimized for
categorical predictors commonly used in single-cell RNA sequencing
analysis. It uses QR decomposition instead of normal equations for
improved numerical stability and handles rank-deficient design matrices
gracefully. The function supports various model specifications
including:

- Simple cell type effects: `PC ~ cell_type`

- Cell type and batch interactions: `PC ~ cell_type * batch`

- Cell type and dataset interactions: `PC ~ cell_type * dataset`

The output format is compatible with `speedglm` results to maintain
consistency with existing plotting and analysis workflows.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
