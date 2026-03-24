# pred_measure_engine.R: internal orchestration engine for validation,
# weighting, computation, and object assembly.

.pred_measure_engine <- function(
  source = c("loo", "insample", "kfold", "test"),
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_insample = NULL,
  measure,
  predperf = NULL,
  fold_id = NULL,
  group_ids = NULL,
  model_name = NULL,
  psis_object = NULL,
  save_psis = FALSE
) {
  state <- .prepare_pred_measure_engine_state(
    source = source,
    measure = measure,
    predperf = predperf,
    ypred = ypred,
    ylp = ylp,
    mupred = mupred,
    ylp_insample = ylp_insample,
    fold_id = fold_id
  )

  weight_info <- .resolve_log_weights(
    source = state$source,
    predperf = state$predperf,
    psis_object = psis_object,
    ylp = state$ylp,
    n_draws = state$n_draws,
    n_obs = state$n_obs
  )

  existing <- .extract_predperf_components(state$predperf)

  base_measure <- .compute_base_measure_if_needed(
    source = state$source,
    ylp = state$ylp,
    ylp_insample = state$ylp_insample,
    normalized_log_weights = weight_info$normalized_log_weights,
    existing_estimates = existing$estimates,
    existing_pointwise = existing$pointwise
  )

  merged <- .accumulate_requested_measures(
    source = state$source,
    measure_names = state$measure,
    y = y,
    ypred = state$ypred,
    mupred = state$mupred,
    ylp = state$ylp,
    normalized_log_weights = weight_info$normalized_log_weights,
    n_draws = state$n_draws,
    n_obs = state$n_obs,
    existing_estimates = existing$estimates,
    existing_pointwise = existing$pointwise,
    base_measure = base_measure
  )

  diagnostics <- .resolve_diagnostics(
    predperf = state$predperf,
    psis_used = weight_info$psis_used
  )
  metadata <- .build_pred_measure_metadata(
    source = state$source,
    predperf = state$predperf,
    group_ids = group_ids,
    model_name = model_name,
    fold_id = state$fold_id
  )

  .new_pred_measure(
    source = state$source,
    estimates = merged$estimates,
    pointwise = merged$pointwise,
    diagnostics = diagnostics,
    metadata = metadata,
    log_weights = weight_info$log_weights,
    psis_used = weight_info$psis_used,
    save_psis = save_psis
  )
}

# internal helper functions ---------------------------------------------------

.prepare_pred_measure_engine_state <- function(
  source,
  measure,
  predperf,
  ypred,
  ylp,
  mupred,
  ylp_insample,
  fold_id
) {
  source <- match.arg(source, c("loo", "insample", "kfold", "test"))
  measure <- .validate_measure_names(measure, .pred_measure_choices)

  if (source == "kfold" && is.null(fold_id)) {
    stop("`fold_id` is required for `kfold_pred_measure()`.")
  }

  .assert_predperf_compatibility(source, predperf)

  if (source != "loo") {
    coerce <- .coerce_nonloo_matrix_input
    ypred <- coerce(ypred)
    ylp <- coerce(ylp)
    mupred <- coerce(mupred)
    ylp_insample <- coerce(ylp_insample)
  }

  dims <- .infer_measure_dims(
    ypred = ypred,
    ylp = ylp,
    mupred = mupred,
    fallback_log_weights = predperf$log_weights,
    fallback_predperf = predperf
  )
  n_draws <- dims[1]
  n_obs <- dims[2]

  checkmate::assert_int(n_draws, lower = 1)
  checkmate::assert_int(n_obs, lower = 1)

  if (source == "kfold") {
    checkmate::assert_vector(fold_id, len = n_obs, any.missing = FALSE)
  }

  list(
    source = source,
    measure = measure,
    predperf = predperf,
    ypred = ypred,
    ylp = ylp,
    mupred = mupred,
    ylp_insample = ylp_insample,
    fold_id = fold_id,
    n_draws = n_draws,
    n_obs = n_obs
  )
}

# internal function to build log weights for all supported sources
.resolve_log_weights <- function(
  source,
  predperf,
  psis_object,
  ylp,
  n_draws,
  n_obs
) {
  if (source == "loo") {
    weight_info <- .compute_loo_log_weights(
      loo = predperf,
      psis_object = psis_object,
      ylp = ylp,
      n_draws = n_draws,
      n_obs = n_obs
    )
  } else {
    log_weights <- predperf$log_weights
    if (is.null(log_weights)) {
      log_weights <- matrix(
        -log(n_draws),
        nrow = n_draws,
        ncol = n_obs
      )
    } else {
      checkmate::assert_matrix(
        log_weights, nrows = n_draws, ncols = n_obs
      )
    }
    weight_info <- list(log_weights = log_weights, psis_used = NULL)
  }

  .assert_finite_log_weight_columns(weight_info$log_weights)

  weight_info$normalized_log_weights <- .normalize_log_weights(
    weight_info$log_weights
  )
  weight_info
}

# internal function to assert each weight column can be normalized
.assert_finite_log_weight_columns <- function(log_weights) {
  col_log_sums <- matrixStats::colLogSumExps(log_weights)
  if (any(!is.finite(col_log_sums))) {
    stop(
      "At least one `log_weights` column has no finite values and cannot be normalized."
    )
  }
  invisible(NULL)
}

# internal function to get components from existing predperf objects
.extract_predperf_components <- function(predperf) {
  if (is.null(predperf)) {
    return(list(estimates = NULL, pointwise = NULL))
  }
  list(
    estimates = predperf$estimates,
    pointwise = predperf$pointwise
  )
}

# internal function to compute measure summaries with reuse of existing columns
.compute_measure_values <- function(
  source,
  measure,
  spec,
  y,
  ypred,
  mupred,
  ylp,
  normalized_log_weights,
  n_draws,
  n_obs,
  existing_estimates,
  existing_pointwise
) {
  pointwise_col <- spec$pointwise_col
  results_col <- spec$results_col

  use_existing_elpd <- source == "loo" &&
    measure %in% c("elpd", "logscore") &&
    !is.null(existing_estimates) &&
    !is.null(existing_pointwise) &&
    results_col %in% rownames(existing_estimates) &&
    pointwise_col %in% colnames(existing_pointwise)

  if (use_existing_elpd) {
    return(list(
      estimate = existing_estimates[results_col, "Estimate"],
      se = existing_estimates[results_col, "SE"],
      pointwise = existing_pointwise[, pointwise_col]
    ))
  }

  args <- .validate_measure_args(
    measure = measure,
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    n_draws = n_draws,
    n_obs = n_obs
  )
  args$log_weights <- normalized_log_weights

  if (!is.null(existing_pointwise) && pointwise_col %in% colnames(existing_pointwise)) {
    args$pointwise <- existing_pointwise[, pointwise_col]
  }

  do.call(spec$summary_fun, args)
}

# internal function to optionally construct base (elpd/p_eff/ic) summaries
.compute_base_measure_if_needed <- function(
  source,
  ylp,
  ylp_insample,
  normalized_log_weights,
  existing_estimates,
  existing_pointwise
) {
  needs_base_measure <- !is.null(ylp) &&
    (is.null(existing_estimates) || is.null(existing_pointwise))
  if (!needs_base_measure) {
    return(NULL)
  }

  is_out_of_sample <- source %in% c("kfold", "test")
  .build_base_measure(
    ylp = ylp,
    normalized_log_weights = normalized_log_weights,
    include_p_eff = source == "loo" || (is_out_of_sample && !is.null(ylp_insample)),
    ylp_insample = if (is_out_of_sample) ylp_insample else NULL
  )
}

# internal function to compute and merge one or multiple requested measures
.accumulate_requested_measures <- function(
  source,
  measure_names,
  y,
  ypred,
  mupred,
  ylp,
  normalized_log_weights,
  n_draws,
  n_obs,
  existing_estimates,
  existing_pointwise,
  base_measure
) {
  merged <- list(estimates = existing_estimates, pointwise = existing_pointwise)

  if (is.null(merged$estimates) && !is.null(base_measure)) {
    merged$estimates <- base_measure$estimates
  }
  if (is.null(merged$pointwise) && !is.null(base_measure)) {
    merged$pointwise <- base_measure$pointwise
  }

  for (measure_name in measure_names) {
    spec <- .get_measure_spec(measure_name)
    output_names <- .resolve_measure_output_names(
      measure = measure_name,
      spec = spec
    )

    measure_values <- .compute_measure_values(
      source = source,
      measure = measure_name,
      spec = spec,
      y = y,
      ypred = ypred,
      mupred = mupred,
      ylp = ylp,
      normalized_log_weights = normalized_log_weights,
      n_draws = n_draws,
      n_obs = n_obs,
      existing_estimates = merged$estimates,
      existing_pointwise = merged$pointwise
    )

    merged <- .merge_measure_outputs(
      output_names = output_names,
      measure_values = measure_values,
      existing_estimates = merged$estimates,
      existing_pointwise = merged$pointwise,
      base_measure = NULL
    )
  }

  merged
}

# internal function to resolve output row/column names for each measure
.resolve_measure_output_names <- function(measure, spec) {
  if (identical(measure, "logscore")) {
    return(list(results_col = "logscore", pointwise_col = "logscore"))
  }
  list(results_col = spec$results_col, pointwise_col = spec$pointwise_col)
}

# internal function to merge previous and current summaries
.merge_measure_outputs <- function(
  output_names,
  measure_values,
  existing_estimates,
  existing_pointwise,
  base_measure
) {
  estimates <- existing_estimates
  pointwise <- existing_pointwise

  if (is.null(estimates) && !is.null(base_measure)) {
    estimates <- base_measure$estimates
  }
  if (is.null(pointwise) && !is.null(base_measure)) {
    pointwise <- base_measure$pointwise
  }

  estimates <- .merge_estimate_row(
    estimates = estimates,
    row_name = output_names$results_col,
    estimate = measure_values$estimate,
    se = measure_values$se
  )
  pointwise <- .merge_pointwise_col(
    pointwise = pointwise,
    col_name = output_names$pointwise_col,
    values = measure_values$pointwise
  )

  list(estimates = estimates, pointwise = pointwise)
}

# internal function to preserve diagnostics with PSIS fallback
.resolve_diagnostics <- function(predperf, psis_used) {
  if (!is.null(predperf) && !is.null(predperf$diagnostics)) {
    return(predperf$diagnostics)
  }
  if (!is.null(psis_used)) {
    return(psis_used$diagnostics)
  }
  NULL
}

# internal function to build metadata consistently
.build_pred_measure_metadata <- function(
  source, predperf, group_ids, model_name, fold_id
) {
  metadata <- if (!is.null(predperf) && !is.null(predperf$metadata)) {
    predperf$metadata
  } else {
    list()
  }

  metadata$source <- source
  if (!is.null(group_ids)) {
    metadata$group_ids <- group_ids
  }
  if (!is.null(model_name)) {
    metadata$model_name <- model_name
  }
  if (source == "kfold") {
    metadata$fold_id <- fold_id
  }

  metadata
}

# internal function to construct a pred_measure object
.new_pred_measure <- function(
  source,
  estimates,
  pointwise,
  diagnostics,
  metadata,
  log_weights,
  psis_used,
  save_psis
) {
  structure(
    list(
      estimates = estimates,
      pointwise = pointwise,
      diagnostics = diagnostics,
      metadata = metadata,
      log_weights = log_weights,
      psis_object = if (source == "loo" && save_psis) psis_used else NULL
    ),
    class = c(
      switch(
        source,
        loo = "loo_pred_measure",
        insample = "insample_pred_measure",
        kfold = "kfold_pred_measure",
        test = "test_pred_measure"
      ),
      "pred_measure"
    ),
    dims = dim(log_weights)
  )
}

# internal function to infer the dimensions of the measure
.infer_measure_dims <- function(
  ypred,
  ylp,
  mupred,
  fallback_log_weights = NULL,
  fallback_predperf = NULL
) {
  matrix_candidates <- list(ypred, ylp, mupred, fallback_log_weights)

  for (x in matrix_candidates) {
    d <- dim(x)

    if (!is.null(d)) {
      checkmate::assert_vector(d, len = 2)
      return(d)
    }
  }

  psis_dims <- attr(fallback_predperf$psis_object, "dims")
  if (!is.null(psis_dims)) {
    return(psis_dims)
  }

  stop(
    "Could not infer dimensions from inputs. Supply at least one matrix input."
  )
}

.validate_measure_args <- function(
  measure,
  y,
  ypred,
  mupred,
  ylp,
  n_draws,
  n_obs
) {
  switch(
    measure,
    mae = ,
    mse = ,
    rmse = ,
    r2 = ,
    acc = ,
    balanced_acc = {
      .validate_measure_arg_not_null(measure, list(y = y, mupred = mupred))
      checkmate::assert_vector(y, len = n_obs, any.missing = FALSE)
      checkmate::assert_matrix(mupred, nrows = n_draws, ncols = n_obs)
      list(y = y, mupred = mupred)
    },
    rps = ,
    crps = ,
    srps = ,
    scrps = {
      .validate_measure_arg_not_null(measure, list(y = y, ypred = ypred))
      checkmate::assert_vector(y, len = n_obs, any.missing = FALSE)
      checkmate::assert_matrix(ypred, nrows = n_draws, ncols = n_obs)
      list(y = y, ypred = ypred)
    },
    elpd = ,
    logscore = ,
    mlpd = {
      .validate_measure_arg_not_null(measure, list(ylp = ylp))
      checkmate::assert_matrix(ylp, nrows = n_draws, ncols = n_obs)
      list(ylp = ylp)
    },
    stop("Unsupported measure: ", measure)
  )
}

# internal function to resolve the log-weights for the LOO case
.compute_loo_log_weights <- function(
  loo, psis_object, ylp, n_draws, n_obs
) {
  log_weights <- loo$log_weights
  psis_loo <- loo$psis_object
  has_psis_arg <- !missing(psis_object) && !is.null(psis_object)
  has_loo_psis <- !is.null(psis_loo)

  if (has_psis_arg && has_loo_psis) {
    warning(
      "Passing both PSIS and loo objects is not advised--defaulting to getting log-weights from loo object."
    )
  }

  psis_used <- if (has_loo_psis) psis_loo else if (has_psis_arg) psis_object else NULL

  if (!is.null(psis_used) && !loo::is.psis(psis_used)) {
    if (has_loo_psis) {
      stop(
        "No valid `psis` object found in provided `loo` object. Make sure you rerun `loo()` with the argument `save_psis = TRUE`.\nIf you want to use equal weighting, call `pred_measure()` instead."
      )
    }
    stop("Provided `psis_object` is not a valid `psis` object.")
  }

  if (is.null(log_weights)) {
    if (is.null(psis_used)) {
      if (is.null(ylp)) {
        stop(
          "No possible way to obtain log-weights. Pass psis object, loo object with the argument `save_psis = TRUE`, or ylp."
        )
      }
      psis_used <- loo::psis(-ylp)
    }
    log_weights <- weights(psis_used, normalize = TRUE, log = TRUE)
  }
  checkmate::assert_matrix(log_weights, nrows = n_draws, ncols = n_obs)
  list(log_weights = log_weights, psis_used = psis_used)
}

.build_base_measure <- function(
  ylp,
  normalized_log_weights,
  include_p_eff = TRUE,
  ylp_insample = NULL
) {
  base_elpd <- .elpd_summary(ylp, normalized_log_weights)
  # TODO: make a measure--use ELPD ptwise
  ic_pointwise <- -2 * base_elpd$pointwise
  ic <- list(
    estimate = sum(ic_pointwise),
    se = ncol(ylp) * .se_helper(ic_pointwise, mean(ic_pointwise), ncol(ylp)),
    pointwise = ic_pointwise
  )

  estimates <- rbind(
    elpd = c(base_elpd$estimate, base_elpd$se),
    ic = c(ic$estimate, ic$se)
  )
  pointwise <- cbind(
    elpd = base_elpd$pointwise,
    ic = ic$pointwise
  )

  if (include_p_eff) {
    lpd_input <- if (is.null(ylp_insample)) ylp else ylp_insample
    checkmate::assert_matrix(lpd_input, ncols = ncol(ylp))
    lpd_pointwise <- matrixStats::colLogSumExps(lpd_input) -
      log(nrow(lpd_input))
    p_eff_pointwise <- lpd_pointwise - base_elpd$pointwise
    p_eff <- list(
      estimate = sum(p_eff_pointwise),
      se = ncol(ylp) *
        .se_helper(p_eff_pointwise, mean(p_eff_pointwise), ncol(ylp)),
      pointwise = p_eff_pointwise
    )
    estimates <- rbind(
      elpd = c(base_elpd$estimate, base_elpd$se),
      p_eff = c(p_eff$estimate, p_eff$se),
      ic = c(ic$estimate, ic$se)
    )
    pointwise <- cbind(
      elpd = base_elpd$pointwise,
      p_eff = p_eff$pointwise,
      ic = ic$pointwise
    )
  }

  colnames(estimates) <- c("Estimate", "SE")
  list(estimates = estimates, pointwise = pointwise)
}

.merge_estimate_row <- function(estimates, row_name, estimate, se) {
  if (is.null(estimates)) {
    return(
      matrix(
        c(estimate, se),
        nrow = 1,
        dimnames = list(row_name, c("Estimate", "SE"))
      )
    )
  }

  # TODO: Should always update? error? warn? message?
  if (row_name %in% rownames(estimates)) {
    estimates[row_name, ] <- c(estimate, se)
    return(estimates)
  }

  rbind(
    estimates,
    matrix(
      c(estimate, se),
      nrow = 1,
      dimnames = list(row_name, colnames(estimates))
    )
  )
}

.merge_pointwise_col <- function(pointwise, col_name, values) {
  if (is.null(pointwise)) {
    return(matrix(values, ncol = 1, dimnames = list(NULL, col_name)))
  }
  if (col_name %in% colnames(pointwise)) {
    pointwise[, col_name] <- values
    return(pointwise)
  }
  cbind(pointwise, matrix(values, ncol = 1, dimnames = list(NULL, col_name)))
}

# internal function to coerce non-LOO matrix inputs to a matrix
# if the input is a vector, it is coerced to a matrix with one row
.coerce_nonloo_matrix_input <- function(x) {
  if (is.null(x) || !is.null(dim(x)) || !is.atomic(x)) {
    return(x)
  }
  matrix(x, nrow = 1)
}

# internal function to assert that the predperf object
# is compatible with the source
.assert_predperf_compatibility <- function(source, predperf) {
  if (is.null(predperf)) return(invisible(NULL))

  expected_classes <- list(
    loo = c("loo", "loo_pred_measure"),
    insample = "insample_pred_measure",
    kfold = "kfold_pred_measure",
    test = "test_pred_measure"
  )

  checkmate::assert_multi_class(predperf, expected_classes[[source]])
  invisible(NULL)
}
