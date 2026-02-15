# Shared parameters for summary functions

Shared parameters for summary functions

## Arguments

- y:

  vector of observed values (n)

- ypred:

  matrix of posterior draws (S x n) of posterior predictive draws

- ylp:

  matrix of posterior draws (S x n) of pointwise log predictive
  densities

- mupred:

  matrix of posterior draws (S x n) of point predictions

- log_weights:

  matrix of standardized loo weights (S x n) on the log scale

- pointwise:

  optional precomputed pointwise squared errors (n)

## Assumptions

`log_weights` are on the log scale and standardized. `y` is a vector of
length `n` and any relevant amongst `mupred`, `ypred`, `ylp`, and
`log_weights` are matrices of size `S x n`.
