# Root mean squared error

Root mean squared error

## Usage

``` r
rmse(y, mupred, log_weights)
```

## Arguments

- y:

  vector of observed values (n)

- mupred:

  matrix of posterior draws (S x n) of point predictions

- log_weights:

  matrix of standardized loo weights (S x n) on the log scale
