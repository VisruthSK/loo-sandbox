# pred_measure_helpers.R: shared helper utilities for weight handling,
# summaries, and measure name lookup.


#' Weighted Mean
#'
#' A wrapper around `stats::weighted.mean` which treats `NULL` weights
#' as missing.
#'
#' @noRd
.loo_weighted_mean <- function(x, log_weights) {
  sum(x * exp(log_weights - matrixStats::logSumExp(log_weights)))
}

.se_helper <- function(x, xbar, n) {
  sqrt(sum((x - xbar)^2) / (n * (n - 1)))
}

#' Simple Summary
#'
#' Function to get common estimate and associated SE.
#'
#' @noRd
.simple_pointwise_summary <- function(pointwise) {
  est <- mean(pointwise)
  list(
    estimate = est,
    se = .se_helper(pointwise, est, length(pointwise)),
    pointwise = pointwise
  )
}

# internal function to normalizes each column of a log-weight matrix so
# the weights sum to 1 on the probability scale.
.normalize_log_weights <- function(log_weights) {
  sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights),
    FUN = "-",
    check.margin = FALSE
  )
}

# internal helper to validate one or multiple measure names
.validate_measure_names <- function(measure, .pred_measure_choices) {
  invalid <- setdiff(measure, .pred_measure_choices)

  if (length(invalid) > 0) {
    stop(sprintf(
      "Invalid measure(s): %s\nMust be one of: %s",
      paste(shQuote(invalid), collapse = ", "),
      paste(shQuote(.pred_measure_choices), collapse = ", ")
    ))
  }

  invisible(measure)
}

# internal helper used by tests and wrapper parity checks
.match_results_column <- function(measure) {
  spec <- .get_measure_spec(measure)
  if (is.null(spec)) {
    stop("Unsupported measure: ", measure)
  }
  spec$results_col
}

# internal helper used by tests and wrapper parity checks
.match_pointwise_column <- function(measure) {
  spec <- .get_measure_spec(measure)
  if (is.null(spec)) {
    stop("Unsupported measure: ", measure)
  }
  spec$pointwise_col
}

# internal helper to ensure required measure arguments are provided
.validate_measure_arg_not_null <- function(measure, arg) {
  arg_name <- if (is.list(arg)) {
    names(arg)
  } else {
    deparse(substitute(arg))
  }

  if (is.list(arg)) {
    missing_args <- vapply(arg, is.null, logical(1))
  } else {
    missing_args <- is.null(arg)
  }

  if (!any(missing_args)) {
    return(invisible(NULL))
  }

  missing_arg_names <- arg_name[missing_args]
  if (length(missing_arg_names) == 1L) {
    stop(
      sprintf(
        "For %s, the `%s` argument must be specified but is currently NULL.",
        measure,
        missing_arg_names
      )
    )
  }

  stop(
    sprintf(
      "For %s, the following arguments must be specified but are currently NULL: %s.",
      measure,
      paste(sprintf("`%s`", missing_arg_names), collapse = ", ")
    )
  )
}