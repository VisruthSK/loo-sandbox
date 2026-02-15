# Shared parameters for pointwise functions

Shared parameters for pointwise functions

## Arguments

- y:

  scalar, leave one out value

- ypred:

  vector (S) of posterior predictive draws

- ylp:

  vector (S) of pointwise log predictive densities

- mupred:

  vector (S) of point predictions

- log_weights:

  vector of standardized loo weights (S) on the log scale

## Assumptions

`log_weights` are on the log scale and standardized. `y` is a scalar and
any relevant amongst `mupred`, `ypred`, `ylp`, and `log_weights` are
vectors of length `S`.
