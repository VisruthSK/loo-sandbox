# Shared parameters for predictive measure wrappers

Shared parameters for predictive measure wrappers

## Arguments

- y:

  vector of observed values (n)

- ypred:

  matrix of posterior predictive draws (S x n)

- mupred:

  matrix of posterior point predictions (S x n)

- ylp:

  matrix of pointwise log predictive densities (S x n)

- ylp_insample:

  optional matrix of in-sample pointwise log predictive densities (S x
  n), used for effective number of parameters in out-of-sample settings

- predperf:

  existing predictive measure object to update/accumulate additional
  metrics

- measure:

  predictive metric to compute

- group_ids:

  optional grouping ids

- model_name:

  optional model label stored in metadata

- fold_id:

  vector of fold ids (required for k-fold wrappers)

- loo:

  `loo` object, preferably fit with `save_psis = TRUE`

- psis_object:

  optional `psis` object used to derive LOO weights

- save_psis:

  logical; if `TRUE`, store the `psis` object in output for LOO paths
