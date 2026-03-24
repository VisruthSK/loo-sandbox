# pred_measure_measures.R: internal pointwise and summary implementations
# for all supported predictive scores and metrics.


# internal constant
.pred_measure_choices <- c("logscore", "mlpd", "elpd", "r2", "rps", "crps",
"srps", "scrps", "mae", "rmse", "mse", "acc", "balanced_acc") # "energy",

#' Shared parameters for pointwise functions
#'
#' @param y scalar, leave one out value
#' @param ypred vector (S) of posterior predictive draws
#' @param ylp vector (S) of pointwise log predictive densities
#' @param mupred vector (S) of point predictions
#' @param log_weights vector of standardized loo weights (S) on the log scale
#'
#' @section Assumptions:
#' `log_weights` are on the log scale and standardized.
#' `y` is a scalar and any relevant amongst `mupred`, `ypred`, `ylp`, and
#' `log_weights` are vectors of length `S`.
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
    ypred <- sort(ypred)
    n <- length(ypred)
    EXX <- -2 * mean(ypred - 2 * ypred * (0:(n - 1)) / (n - 1))
  } else {
    weights <- exp(log_weights - matrixStats::logSumExp(log_weights))
    EXy <- sum(weights * abs(y - ypred))
    ord <- order(ypred)
    ypred <- ypred[ord]
    weights <- weights[ord]
    cw <- (cumsum(weights) - weights) / (1 - weights)
    EXX <- -2 * sum(weights * (ypred - 2 * ypred * cw))
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
      seq_along(y),
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
      seq_along(y),
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
      seq_along(y),
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
      seq_along(y),
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

  (pointwise *
    n /
    (length(cls_counts) *
      as.numeric(cls_counts[match(y, names(cls_counts))]))) |>
    .simple_pointwise_summary() |>
    modifyList(list(pointwise = pointwise))
}

# internal function to get the measure specification
# @noRd
# @param measure The measure used.
# @return The measure specification.
.measure_spec <- list(
  elpd = list(
    fun = .pointwise_elpd,
    summary_fun = .elpd_summary,
    pointwise_col = "elpd",
    results_col = "elpd"
  ),
  logscore = list(
    fun = .pointwise_elpd,
    summary_fun = .elpd_summary,
    pointwise_col = "elpd",
    results_col = "elpd"
  ),
  mlpd = list(
    fun = .pointwise_elpd,
    summary_fun = .mlpd_summary,
    pointwise_col = "elpd",
    results_col = "mlpd"
  ),
  mae = list(
    fun = .pointwise_squared_error,
    summary_fun = .mae_summary,
    pointwise_col = "squared_error",
    results_col = "mae"
  ),
  r2 = list(
    fun = .pointwise_squared_error,
    summary_fun = .r2_summary,
    pointwise_col = "squared_error",
    results_col = "r2"
  ),
  rmse = list(
    fun = .pointwise_squared_error,
    summary_fun = .rmse_summary,
    pointwise_col = "squared_error",
    results_col = "rmse"
  ),
  mse = list(
    fun = .pointwise_squared_error,
    summary_fun = .mse_summary,
    pointwise_col = "squared_error",
    results_col = "mse"
  ),
  acc = list(
    fun = .pointwise_accuracy,
    summary_fun = .accuracy_summary,
    pointwise_col = "accuracy",
    results_col = "acc"
  ),
  balanced_acc = list(
    fun = .pointwise_accuracy,
    summary_fun = .balanced_accuracy_summary,
    pointwise_col = "accuracy",
    results_col = "bal_acc"
  ),
  rps = list(
    fun = .pointwise_rps,
    summary_fun = .rps_summary,
    pointwise_col = "rps",
    results_col = "rps"
  ),
  crps = list(
    fun = .pointwise_rps,
    summary_fun = .rps_summary,
    pointwise_col = "crps",
    results_col = "crps"
  ),
  srps = list(
    fun = .pointwise_rps,
    summary_fun = .srps_summary,
    pointwise_col = "srps",
    results_col = "srps"
  ),
  scrps = list(
    fun = .pointwise_rps,
    summary_fun = .srps_summary,
    pointwise_col = "scrps",
    results_col = "scrps"
  )
)
.get_measure_spec <- function(measure) .measure_spec[[measure]]
