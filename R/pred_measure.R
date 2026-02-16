#' In-Sample Predictive Measure
#'
#' @param y Placeholder
#' @param ypred Placeholder
#' @param mupred Placeholder
#' @param ylp Placeholder
#' @param predperf Placeholder
#' @param measure Placeholder
#' @param group_ids Placeholder
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
  measure = .pred_measure_choices(),
  group_ids = NULL
) {
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
    psis_object = NULL,
    save_psis = FALSE
  )
}

#' K-fold Predictive Measure
#'
#' @param y Placeholder
#' @param ypred Placeholder
#' @param mupred Placeholder
#' @param ylp Placeholder
#' @param ylp_insample Placeholder
#' @param fold_id Placeholder
#' @param predperf Placeholder
#' @param measure Placeholder
#' @param group_ids Placeholder
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
  measure = .pred_measure_choices(),
  group_ids = NULL
) {
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
    psis_object = NULL,
    save_psis = FALSE
  )
}

#' Test-Set Predictive Measure
#'
#' @param y Placeholder
#' @param ypred Placeholder
#' @param mupred Placeholder
#' @param ylp Placeholder
#' @param ylp_insample Placeholder
#' @param predperf Placeholder
#' @param measure Placeholder
#' @param group_ids Placeholder
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
  measure = .pred_measure_choices(),
  group_ids = NULL
) {
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
    psis_object = NULL,
    save_psis = FALSE
  )
}

#' @export
print.pred_measure <- function(x, digits = 1, ...) {
  dims <- attr(x, "dims")
  if (is.null(dims) && !is.null(x$log_weights)) {
    dims <- dim(x$log_weights)
  }
  source <- .pred_measure_source_label(x)

  cat("\n")
  if (!is.null(dims) && length(dims) == 2) {
    cat(
      sprintf(
        "Computed from %s posterior draws and %s observations.\n",
        dims[1],
        dims[2]
      )
    )
  }
  cat(sprintf("Data source: %s\n\n", source))
  print(loo:::.fr(as.data.frame(x$estimates), digits), quote = FALSE)
  invisible(x)
}

#' @export
print.insample_pred_measure <- function(x, digits = 1, ...) {
  print.pred_measure(x, digits = digits, ...)
}

#' @export
print.kfold_pred_measure <- function(x, digits = 1, ...) {
  print.pred_measure(x, digits = digits, ...)
  if (!is.null(x$metadata$fold_id)) {
    cat(sprintf("Folds: %s\n", length(unique(x$metadata$fold_id))))
  }
  invisible(x)
}

#' @export
print.test_pred_measure <- function(x, digits = 1, ...) {
  print.pred_measure(x, digits = digits, ...)
}

# ----------------------------- Engine -----------------------------

.pred_measure_choices <- function() {
  c(
    "logscore",
    "elpd",
    "r2",
    "rps",
    "crps",
    "srps",
    "scrps",
    "mae",
    "rmse",
    "mse",
    "acc",
    # "energy",
    "balanced_acc"
  )
}

.pred_measure_source_label <- function(x) {
  cls <- class(x)
  if ("loo_pred_measure" %in% cls) {
    return("loo")
  }
  if ("insample_pred_measure" %in% cls) {
    return("in-sample")
  }
  if ("kfold_pred_measure" %in% cls) {
    return("k-fold")
  }
  if ("test_pred_measure" %in% cls) {
    return("test")
  }
  if (!is.null(x$metadata$source)) {
    return(as.character(x$metadata$source))
  }
  "unknown"
}

.infer_measure_dims <- function(ypred, ylp, mupred, fallback_log_weights = NULL) {
  if (!is.null(ypred)) {
    return(c(nrow(ypred), ncol(ypred)))
  }
  if (!is.null(ylp)) {
    return(c(nrow(ylp), ncol(ylp)))
  }
  if (!is.null(mupred)) {
    return(c(nrow(mupred), ncol(mupred)))
  }
  if (!is.null(fallback_log_weights)) {
    return(c(nrow(fallback_log_weights), ncol(fallback_log_weights)))
  }
  stop("Could not infer dimensions from inputs. Supply at least one matrix input.")
}

.validate_measure_args <- function(
  measure,
  y,
  ypred,
  mupred,
  ylp,
  S,
  n
) {
  if (
    measure %in%
      c(
        "mae",
        "mse",
        "rmse",
        "r2",
        "acc",
        "balanced_acc"
      )
  ) {
    checkmate::assert_vector(y, len = n, any.missing = FALSE)
    checkmate::assert_matrix(mupred, nrows = S, ncols = n)
    return(list(y = y, mupred = mupred))
  }

  if (
    measure %in%
      c(
        "rps",
        "crps",
        "srps",
        "scrps"
      )
  ) {
    checkmate::assert_vector(y, len = n, any.missing = FALSE)
    checkmate::assert_matrix(ypred, nrows = S, ncols = n)
    return(list(y = y, ypred = ypred))
  }

  if (
    measure %in%
      c(
        "elpd",
        "logscore"
      )
  ) {
    checkmate::assert_matrix(ylp, nrows = S, ncols = n)
    return(list(ylp = ylp))
  }

  stop("Unsupported measure: ", measure)
}

.resolve_loo_log_weights <- function(loo, psis_object, ylp, S, n) {
  log_weights <- loo$log_weights
  psis_loo <- loo$psis_object
  has_psis_arg <- !missing(psis_object) && !is.null(psis_object)
  has_loo_psis <- !is.null(psis_loo)
  psis_used <- NULL

  if (has_psis_arg && has_loo_psis) {
    warning(
      "Passing both PSIS and loo objects is not advised--defaulting to getting log-weights from loo object."
    )
    psis_used <- psis_loo
  } else if (has_loo_psis) {
    if (!is.psis(psis_loo)) {
      stop(
        "No valid `psis` object found in provided `loo` object. Make sure you rerun `loo()` with the argument `save_psis = TRUE`.\nIf you want to use equal weighting, call `pred_measure()` instead."
      )
    }
    psis_used <- psis_loo
  } else if (has_psis_arg) {
    if (!is.psis(psis_object)) {
      stop("Provided `psis_object` is not a valid `psis` object.")
    }
    psis_used <- psis_object
  } else if (is.null(log_weights)) {
    if (is.null(ylp)) {
      stop(
        "No possible way to obtain log-weights. Pass psis object, loo object with the argument `save_psis = TRUE`, or ylp."
      )
    }
    psis_used <- psis(ylp)
  }

  if (is.null(log_weights)) {
    log_weights <- weights(psis_used)
  }
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  list(log_weights = log_weights, psis_used = psis_used)
}

.equal_log_weights <- function(S, n) {
  matrix(-log(S), nrow = S, ncol = n)
}

.build_base_measure <- function(
  ylp,
  log_weights_std,
  include_p_eff = TRUE,
  ylp_insample = NULL
) {
  base_elpd <- .elpd_summary(ylp, log_weights_std)
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
    lpd_pointwise <- matrixStats::colLogSumExps(lpd_input) - log(nrow(lpd_input))
    p_eff_pointwise <- lpd_pointwise - base_elpd$pointwise
    p_eff <- list(
      estimate = sum(p_eff_pointwise),
      se = ncol(ylp) * .se_helper(p_eff_pointwise, mean(p_eff_pointwise), ncol(ylp)),
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

  if (row_name %in% rownames(estimates)) {
    estimates[row_name, ] <- c(estimate, se)
    return(estimates)
  }

  new_row <- matrix(
    c(estimate, se),
    nrow = 1,
    dimnames = list(row_name, colnames(estimates))
  )
  rbind(estimates, new_row)
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
  psis_object = NULL,
  save_psis = FALSE
) {
  source <- match.arg(source)
  measure <- match.arg(measure, .pred_measure_choices())

  if (source == "kfold" && is.null(fold_id)) {
    stop("`fold_id` is required for `kfold_pred_measure()`.")
  }

  if (!is.null(predperf) && !is.list(predperf)) {
    stop("`predperf` must be a prediction measure object.")
  }

  fallback_log_weights <- if (!is.null(predperf)) predperf$log_weights else NULL
  dims <- .infer_measure_dims(ypred, ylp, mupred, fallback_log_weights)
  S <- dims[1]
  n <- dims[2]

  if (source == "kfold") {
    checkmate::assert_vector(fold_id, len = n, any.missing = FALSE)
  }

  if (source == "loo") {
    lw_info <- .resolve_loo_log_weights(
      loo = predperf,
      psis_object = psis_object,
      ylp = ylp,
      S = S,
      n = n
    )
    log_weights <- lw_info$log_weights
    psis_used <- lw_info$psis_used
  } else if (!is.null(predperf) && !is.null(predperf$log_weights)) {
    log_weights <- predperf$log_weights
    checkmate::assert_matrix(log_weights, nrows = S, ncols = n)
    psis_used <- NULL
  } else {
    log_weights <- .equal_log_weights(S, n)
    psis_used <- NULL
  }

  log_weights_std <- .standardize_log_weights(log_weights)

  pointwise_col <- .match_pointwise_column(measure)
  results_col <- .match_results_column(measure)
  pred_fun <- .loo_predictive_measure_fun(measure)

  existing_estimates <- if (!is.null(predperf)) predperf$estimates else NULL
  existing_pointwise <- if (!is.null(predperf)) predperf$pointwise else NULL

  use_existing_elpd <- source == "loo" &&
    measure == "elpd" &&
    !is.null(existing_estimates) &&
    !is.null(existing_pointwise) &&
    results_col %in% rownames(existing_estimates) &&
    pointwise_col %in% colnames(existing_pointwise)

  if (use_existing_elpd) {
    measure_values <- list(
      estimate = existing_estimates[results_col, "Estimate"],
      se = existing_estimates[results_col, "SE"],
      pointwise = existing_pointwise[, pointwise_col]
    )
  } else {
    args <- .validate_measure_args(
      measure = measure,
      y = y,
      ypred = ypred,
      mupred = mupred,
      ylp = ylp,
      S = S,
      n = n
    )
    args$log_weights <- log_weights_std
    if (!is.null(existing_pointwise) && pointwise_col %in% colnames(existing_pointwise)) {
      args$pointwise <- existing_pointwise[, pointwise_col]
    }
    measure_values <- do.call(pred_fun, args)
  }

  include_p_eff <- source == "loo" || (source %in% c("kfold", "test") && !is.null(ylp_insample))
  base_measure <- NULL
  if (!is.null(ylp) && (is.null(existing_estimates) || is.null(existing_pointwise))) {
    base_measure <- .build_base_measure(
      ylp = ylp,
      log_weights_std = log_weights_std,
      include_p_eff = include_p_eff,
      ylp_insample = if (source %in% c("kfold", "test")) ylp_insample else NULL
    )
  }

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
    row_name = results_col,
    estimate = measure_values$estimate,
    se = measure_values$se
  )
  pointwise <- .merge_pointwise_col(
    pointwise = pointwise,
    col_name = pointwise_col,
    values = measure_values$pointwise
  )

  diagnostics <- if (!is.null(predperf) && !is.null(predperf$diagnostics)) {
    predperf$diagnostics
  } else if (!is.null(psis_used)) {
    psis_used$diagnostics
  } else {
    NULL
  }

  metadata <- if (!is.null(predperf) && !is.null(predperf$metadata)) {
    predperf$metadata
  } else {
    list()
  }
  metadata$source <- source
  if (!is.null(group_ids)) {
    metadata$group_ids <- group_ids
  }
  if (source == "kfold") {
    metadata$fold_id <- fold_id
  }

  class_name <- switch(
    source,
    loo = "loo_pred_measure",
    insample = "insample_pred_measure",
    kfold = "kfold_pred_measure",
    test = "test_pred_measure"
  )

  structure(
    list(
      estimates = estimates,
      pointwise = pointwise,
      diagnostics = diagnostics,
      metadata = metadata,
      log_weights = log_weights,
      psis_object = if (source == "loo" && save_psis) psis_used else NULL
    ),
    class = c(class_name, "pred_measure"),
    dims = dim(log_weights)
  )
}
