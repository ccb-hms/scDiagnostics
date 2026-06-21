# Decompose R-squared by Model Components

Decomposes the total R-squared from a linear model into individual
components representing the variance explained by main effects (cell
type, batch/dataset) and their interaction using sequential sum of
squares.

## Usage

``` r
decomposeR2(pc, indep_var, df)
```

## Arguments

- pc:

  A character string specifying the principal component column name in
  the data frame.

- indep_var:

  A character string specifying the independent variable specification.
  Options are "cell_type", "cell_type \* batch", or "cell_type \*
  dataset".

- df:

  A data frame containing the principal component scores and categorical
  predictors. Must include columns for the specified PC and predictor
  variables.

## Value

A named list containing R-squared components:

- cell_type:

  Numeric value representing the R-squared explained by cell type main
  effect.

- batch/dataset:

  Numeric value representing the R-squared explained by batch or dataset
  main effect (only for interaction models).

- interaction:

  Numeric value representing the R-squared explained by the interaction
  term (only for interaction models).

## Details

This function performs R-squared decomposition using a sequential sum of
squares approach, which partitions the total explained variance into
additive components. The decomposition follows the hierarchical
structure:

1.  Cell type main effect

2.  Batch/dataset main effect (after accounting for cell type)

3.  Cell type × batch/dataset interaction (after accounting for main
    effects)

For simple cell type models, only the cell type component is returned.
For interaction models, all three components are computed. The method
uses group means and residual analysis to avoid computationally
expensive matrix operations while maintaining mathematical accuracy
equivalent to ANOVA decomposition.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>
