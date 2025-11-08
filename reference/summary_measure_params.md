# Shared parameters for summary functions

Shared parameters for summary functions

## Arguments

- y:

  vector of observed values (n)

- ypred:

  matrix of posterior draws (S x n) of posterior predictive draws

- ylp:

  matrix of posterior draws (S x n) of pointwise LOO log predictive
  densities

- mupred:

  matrix of posterior draws (S x n) of point predictions

- log_weights:

  matrix of loo weights (S x n) on the log scale

- pointwise:

  optional precomputed pointwise squared errors (n)
