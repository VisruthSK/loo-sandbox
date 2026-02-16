# In-Sample Predictive Measure

In-Sample Predictive Measure

## Usage

``` r
insample_pred_measure(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  predperf = NULL,
  measure = .pred_measure_choices(),
  group_ids = NULL,
  model_name = NULL
)
```

## Arguments

- y:

  vector of observed values (n)

- ypred:

  matrix of posterior predictive draws (S x n)

- mupred:

  matrix of posterior point predictions (S x n)

- ylp:

  matrix of pointwise log predictive densities (S x n)

- predperf:

  existing predictive measure object to update/accumulate additional
  metrics

- measure:

  predictive metric to compute

- group_ids:

  optional grouping ids

- model_name:

  optional model label stored in metadata

## Value

Placeholder
