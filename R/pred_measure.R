# pred_measure.R: primary user-facing constructors for predictive
# performance objects across data sources.


#' Shared parameters for predictive measure wrappers
#'
#' @param y vector of observed values (n)
#' @param ypred matrix of posterior predictive draws (S x n)
#' @param mupred matrix of posterior point predictions (S x n)
#' @param ylp matrix of pointwise log predictive densities (S x n)
#' @param ylp_insample optional matrix of in-sample pointwise log predictive
#'   densities (S x n), used for effective number of parameters in out-of-sample
#'   settings
#' @param predperf existing predictive measure object to update/accumulate
#'   additional metrics
#' @param measure predictive metric to compute
#' @param group_ids optional grouping ids
#' @param model_name optional model label stored in metadata
#' @param fold_id vector of fold ids (required for k-fold wrappers)
#' @param loo `loo` object, preferably fit with `save_psis = TRUE`
#' @param psis_object optional `psis` object used to derive LOO weights
#' @param save_psis logical; if `TRUE`, store the `psis` object in output for
#'   LOO paths
#'
#' @keywords internal
#' @name pred_measure_params
NULL

#' In-Sample Predictive Measure
#'
#' @inheritParams pred_measure_params
#'
#' @return Placeholder
#'
#' @export
insample_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  predperf = NULL,
  measure = "logscore",
  group_ids = NULL,
  model_name = NULL
) {
  measure <- .validate_measure_names(measure, .pred_measure_choices)

  .pred_measure_engine(
    source = "insample",
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_insample = NULL,
    measure = measure,
    predperf = predperf,
    fold_id = NULL,
    group_ids = group_ids,
    model_name = model_name,
    psis_object = NULL,
    save_psis = FALSE
  )
}

#' K-fold Predictive Measure
#'
#' @inheritParams pred_measure_params
#'
#' @return Placeholder
#'
#' @export
kfold_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_insample = NULL,
  fold_id = NULL,
  predperf = NULL,
  measure = "logscore",
  group_ids = NULL,
  model_name = NULL
) {
  measure <- .validate_measure_names(measure, .pred_measure_choices)

  .pred_measure_engine(
    source = "kfold",
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_insample = ylp_insample,
    measure = measure,
    predperf = predperf,
    fold_id = fold_id,
    group_ids = group_ids,
    model_name = model_name,
    psis_object = NULL,
    save_psis = FALSE
  )
}

#' Test-Set Predictive Measure
#'
#' @inheritParams pred_measure_params
#'
#' @return Placeholder
#'
#' @export
test_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_insample = NULL,
  predperf = NULL,
  measure = "logscore",
  group_ids = NULL,
  model_name = NULL
) {
  measure <- .validate_measure_names(measure, .pred_measure_choices)

  .pred_measure_engine(
    source = "test",
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_insample = ylp_insample,
    measure = measure,
    predperf = predperf,
    fold_id = NULL,
    group_ids = group_ids,
    model_name = model_name,
    psis_object = NULL,
    save_psis = FALSE
  )
}

# TODO: document function--look at loo #281
#' LOO Predictive Measure
#'
#' @inheritParams pred_measure_params
#'
#' @return Placeholder
#'
#' @export
loo_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measure = "logscore",
  group_ids = NULL,
  model_name = NULL,
  loo = NULL, # This `loo` object needs to be fit with `save_psis = TRUE`
  psis_object = NULL,
  save_psis = FALSE
) {

  measure <- .validate_measure_names(measure, .pred_measure_choices)

  .pred_measure_engine(
    source = "loo",
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_insample = NULL,
    measure = measure,
    predperf = loo,
    fold_id = NULL,
    group_ids = group_ids,
    model_name = model_name,
    psis_object = psis_object,
    save_psis = save_psis
  )
}
