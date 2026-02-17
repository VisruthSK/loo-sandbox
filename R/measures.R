#' Identify measure function function by `measure`
#'
#' @noRd
#' @param measure The measure used.
#' @return The function used to compute predictive error or accuracy specified
#' by the argument `measure`.
.loo_predictive_measure_fun <- function(measure) {
  switch(
    measure,
    "elpd" = .elpd_summary,
    "logscore" = .elpd_summary,
    "mlpd" = .mlpd_summary,
    "mae" = .mae_summary,
    "r2" = .r2_summary,
    "rmse" = .rmse_summary,
    "mse" = .mse_summary,
    "acc" = .accuracy_summary,
    "balanced_acc" = .balanced_accuracy_summary,
    "rps" = .rps_summary,
    "crps" = .rps_summary,
    "srps" = .srps_summary,
    "scrps" = .srps_summary,
    # , "energy" = .energy
  )
}

#' Identify pointwise function by `measure`
#'
#' @noRd
#' @param measure The measure used.
#' @return The function used to compute pointwise values by the argument `measure`.
.match_pointwise_function <- function(measure) {
  switch(
    measure,
    "elpd" = .pointwise_elpd,
    "logscore" = .pointwise_elpd,
    "mlpd" = .pointwise_elpd,
    "mae" = .pointwise_squared_error,
    "r2" = .pointwise_squared_error,
    "rmse" = .pointwise_squared_error,
    "mse" = .pointwise_squared_error,
    "acc" = .pointwise_accuracy,
    "balanced_acc" = .pointwise_accuracy,
    "rps" = .pointwise_rps,
    "crps" = .pointwise_rps,
    "srps" = .pointwise_rps,
    "scrps" = .pointwise_rps
  )
}

#' Match measure to pointwise column name
#'
#' @noRd
#' @param measure The measure used.
#' @return The column name in the loo pointwise matrix for the given measure.
.match_pointwise_column <- function(measure) {
  switch(
    measure,
    "elpd" = "elpd",
    "logscore" = "elpd",
    "mlpd" = "elpd",
    "mae" = "squared_error",
    "r2" = "squared_error",
    "rmse" = "squared_error",
    "mse" = "squared_error",
    "acc" = "accuracy",
    "balanced_acc" = "accuracy",
    # TODO: check if these pointwise values are actually disparate
    "rps" = "rps",
    "crps" = "crps",
    "srps" = "srps",
    "scrps" = "scrps"
  )
}

#' Match measure to results column name
#'
#' @noRd
#' @param measure The measure used.
#' @return The column name in the loo results matrix for the given measure.
.match_results_column <- function(measure) {
  switch(
    measure,
    "elpd" = "elpd",
    "logscore" = "elpd",
    "mlpd" = "mlpd",
    "mae" = "mae",
    "r2" = "r2",
    "rmse" = "rmse",
    "mse" = "mse",
    "acc" = "acc",
    "balanced_acc" = "bal_acc",
    "rps" = "rps",
    "crps" = "crps",
    "srps" = "srps",
    "scrps" = "scrps"
  )
}

#' Shared parameters for pointwise functions
#'
#' @param y scalar, leave one out value
#' @param ypred vector (S) of posterior predictive draws
#' @param ylp vector (S) of pointwise log predictive densities
#' @param mupred vector (S) of point predictions
#' @param log_weights vector of standardized loo weights (S) on the log scale
#'
#' @section Assumptions:
#' `log_weights` are on the log scale and standardized. `y` is a scalar and any relevant amongst `mupred`, `ypred`, `ylp`, and `log_weights` are vectors of length `S`.
#'
#' @keywords internal
#' @name pointwise_measure_params
NULL

#' Shared parameters for summary functions
#'
#' @param y vector of observed values (n)
#' @param ypred matrix of posterior draws (S x n) of posterior predictive draws
#' @param ylp matrix of posterior draws (S x n) of pointwise log predictive densities
#' @param mupred matrix of posterior draws (S x n) of point predictions
#' @param log_weights matrix of standardized loo weights (S x n) on the log scale
#'
#' @param pointwise optional precomputed pointwise squared errors (n)
#'
#' @section Assumptions:
#' `log_weights` are on the log scale and standardized. `y` is a vector of length `n` and any relevant amongst `mupred`, `ypred`, `ylp`, and `log_weights` are matrices of size `S x n`.
#'
#' @keywords internal
#' @name summary_measure_params
NULL

# ----------------------------- Scores -------------------------------

# https://github.com/stan-dev/loo/blob/master/R/E_loo.R#L260
# https://github.com/stan-dev/loo/blob/master/R/helpers.R#L36
# https://github.com/stan-dev/loo/blob/master/R/loo.R#L427

# TODO: .energy multivariate CRPS, add later

#' Pointwise log score
#'
#' @noRd
#' @inheritParams pointwise_measure_params
#' @inheritSection pointwise_measure_params Assumptions
.pointwise_elpd <- function(ylp, log_weights) {
  if (is.null(log_weights)) {
    matrixStats::logSumExp(ylp) - log(length(ylp))
  } else {
    matrixStats::logSumExp(ylp + log_weights)
  }
}

#' Expected log-predictive density
#'
#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.elpd_summary <- function(ylp, log_weights, pointwise = NULL) {
  n <- ncol(ylp)
  pointwise <- if (is.null(pointwise)) {
    vapply(
      seq_len(n),
      function(i) {
        .pointwise_elpd(ylp[, i], log_weights[, i])
      },
      numeric(1)
    )
  } else {
    pointwise
  }
  list(
    estimate = sum(pointwise),
    se = n * .se_helper(pointwise, mean(pointwise), n),
    pointwise = pointwise
  )
}

#' Mean log-predictive density
#'
#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.mlpd_summary <- function(ylp, log_weights, pointwise = NULL) {
  n <- ncol(ylp)
  l <- .elpd_summary(ylp, log_weights, pointwise)
  # tranformation of elpd estimates, same pointwise values
  list(
    estimate = l$estimate / n,
    se = l$se / n,
    pointwise = l$pointwise
  )
}

#' Pointwise CRPS
#'
#' @noRd
#' @param scaled logical. If true, computes SRPS/SCRPS
#' @inheritParams pointwise_measure_params
#' @inheritSection pointwise_measure_params Assumptions
.pointwise_rps <- function(y, ypred, log_weights, scaled) {
  if (is.null(log_weights)) {
    EXy <- mean(abs(y - ypred))
    y <- sort(y)
    n <- length(y)
    EXX <- -2 * mean(y - 2 * y * (0:(n - 1)) / (n - 1))
  } else {
    EXy <- sum(log_weights * abs(y - ypred))
    ord <- order(y)
    y <- y[ord]
    log_weights <- log_weights[ord]
    cw <- (cumsum(log_weights) - log_weights) / (1 - log_weights)
    EXX <- -2 * sum(log_weights * (y - 2 * y * cw))
  }

  if (!scaled) {
    # Gneiting & Raftery (2007)
    rps <- -EXy + 0.5 * EXX
  } else {
    # Scaled version by Bolin & Wallin (2023)
    rps <- -EXy / EXX - 0.5 * log(EXX)
  }

  rps
}

#' @noRd
#' @param scaled logical. If true, computes SRPS/SCRPS
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.rps_summary <- function(
  y,
  ypred,
  log_weights,
  scaled = FALSE,
  pointwise = NULL
) {
  n <- length(y)

  pointwise <- if (is.null(pointwise)) {
    vapply(
      seq_len(n),
      function(i) {
        .pointwise_rps(
          y[i],
          ypred[, i],
          log_weights[, i],
          scaled
        )
      },
      numeric(1)
    )
  } else {
    pointwise
  }

  .simple_pointwise_summary(pointwise)
}

#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.srps_summary <- function(y, ypred, log_weights, pointwise = NULL) {
  .rps_summary(
    y = y,
    ypred = ypred,
    log_weights = log_weights,
    scaled = TRUE,
    pointwise = pointwise
  )
}

# ----------------------------- Metrics -----------------------------

#' Pointwise squared error
#'
#' @noRd
#' @inheritParams pointwise_measure_params
#' @inheritSection pointwise_measure_params Assumptions
.pointwise_squared_error <- function(y, mupred, log_weights) {
  (y - .loo_weighted_mean(mupred, log_weights))^2
}

#' Mean absolute error
#'
#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.mae_summary <- function(y, mupred, log_weights, pointwise = NULL) {
  pointwise <- if (is.null(pointwise)) {
    vapply(
      seq_len(length(y)),
      function(i) {
        .pointwise_squared_error(
          y[i],
          mupred[, i],
          log_weights[, i]
        )
      },
      numeric(1)
    )
  } else {
    pointwise
  }
  .simple_pointwise_summary(sqrt(pointwise)) |>
    modifyList(list(pointwise = pointwise)) # store squared error
}

#' Mean squared error
#'
#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.mse_summary <- function(
  y,
  mupred,
  log_weights,
  pointwise = NULL
) {
  pointwise <- if (is.null(pointwise)) {
    vapply(
      seq_len(length(y)),
      function(i) {
        .pointwise_squared_error(
          y[i],
          mupred[, i],
          log_weights[, i]
        )
      },
      numeric(1)
    )
  } else {
    pointwise
  }
  .simple_pointwise_summary(pointwise)
}

#' Root mean squared error
#'
#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.rmse_summary <- function(
  y,
  mupred,
  log_weights,
  pointwise = NULL
) {
  pointwise <- if (is.null(pointwise)) {
    vapply(
      seq_len(length(y)),
      function(i) {
        .pointwise_squared_error(
          y[i],
          mupred[, i],
          log_weights[, i]
        )
      },
      numeric(1)
    )
  } else {
    pointwise
  }
  mean_mse <- mean(pointwise)
  list(
    estimate = sqrt(mean_mse),
    se = .se_helper(pointwise, mean_mse, length(y)) / sqrt(mean_mse) / 2, # Comes from the first order Taylor approx.
    pointwise = pointwise
  )
}

#' R^2
#'
#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.r2_summary <- function(
  y,
  mupred,
  log_weights,
  pointwise = NULL
) {
  n <- length(y)
  pointwise <- if (is.null(pointwise)) {
    vapply(
      seq_len(n),
      function(i) {
        .pointwise_squared_error(
          y[i],
          mupred[, i],
          log_weights[, i]
        )
      },
      numeric(1)
    )
  } else {
    pointwise
  }

  mse_loo <- mean(pointwise)
  se_mse_loo <- .se_helper(pointwise, mse_loo, n)

  squared_error_y_pointwise <- (y - mean(y))^2
  mse_y <- mean(squared_error_y_pointwise)

  se_r2 <- sqrt(
    se_mse_loo^2 -
      2 *
        (mse_loo / mse_y) *
        cov(pointwise, squared_error_y_pointwise) /
        n +
      (mse_loo^2 / mse_y^2) * var(squared_error_y_pointwise) / n
  ) /
    mse_y

  list(
    estimate = 1 - mse_loo / mse_y,
    se = se_r2,
    pointwise = pointwise
  )
}

#' Classification accuracy
#'
#' Assumes values in `mupred` only take on 0 or 1
#'
#' @noRd
#' @inheritParams pointwise_measure_params
#' @inheritSection pointwise_measure_params Assumptions
.pointwise_accuracy <- function(y, mupred, log_weights) {
  .loo_weighted_mean(mupred == y, log_weights)
}

#' Classification accuracy
#'
#' Assumes values in `mupred` only take on 0 or 1
#'
#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.accuracy_summary <- function(y, mupred, log_weights, pointwise = NULL) {
  checkmate::assert_subset(mupred, choices = c(0, 1))
  pointwise <- if (is.null(pointwise)) {
    vapply(
      seq_len(length(y)),
      function(i) {
        .pointwise_accuracy(
          y[i],
          mupred[, i],
          log_weights[, i]
        )
      },
      numeric(1)
    )
  } else {
    pointwise
  }

  .simple_pointwise_summary(pointwise)
}

#' Balanced classification accuracy
#'
#' @noRd
#' @inheritParams summary_measure_params
#' @inheritSection summary_measure_params Assumptions
.balanced_accuracy_summary <- function(
  y,
  mupred,
  log_weights,
  pointwise = NULL
) {
  # TODO: check implementation again
  n <- length(y)
  cls_counts <- table(y)

  pointwise <- if (is.null(pointwise)) {
    vapply(
      seq_len(n),
      function(i) {
        .pointwise_accuracy(
          y[i],
          mupred[, i],
          log_weights[, i]
        )
      },
      numeric(1)
    )
  } else {
    pointwise
  }

  .simple_pointwise_summary(
    pointwise *
      n /
      (length(cls_counts) * as.numeric(cls_counts[match(y, names(cls_counts))]))
  )
}

# ----------------------------- Helpers -----------------------------

#' Weighted Mean
#'
#' A wrapper around `stats::weighted.mean` which treats `NULL` weights as missing.
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

.standardize_log_weights <- function(log_weights) {
  sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights),
    FUN = "-",
    check.margin = FALSE
  )
}
