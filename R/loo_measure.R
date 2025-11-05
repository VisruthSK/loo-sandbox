#' @export
loo_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  loo = NULL,
  log_weights = NULL,
  measure = c(
    "logscore",
    "elpd",
    "r2",
    "rps",
    "crps",
    # TODO: is srps option needed
    "srps",
    "scrps",
    "mae",
    "rmse",
    "mse",
    "acc",
    # "energy",
    "balanced_acc"
  ),
  group_ids = NULL,
  psis_object = NULL,
  save_psis = FALSE
) {
  # TODO: document function

  n <- length(y)
  # TODO: refactor this--ypred etc could be functions
  S <- if (!is.null(ypred)) {
    nrow(ypred)
  } else if (!is.null(ylp)) {
    nrow(ylp)
  } else if (!is.null(mupred)) {
    nrow(mupred)
  }
  measure <- match.arg(measure)
  pred_fun <- .loo_predictive_measure_fun(measure)

  # TODO: check warning msgs
  # check appropriate arguments for measure
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
    checkmate::assert_matrix(mupred, nrows = S, ncols = n)
    args <- list(y, mupred)
  } else if (
    measure %in%
      c(
        "rps",
        "crps",
        # "energy",
        "srps",
        "scrps"
      )
  ) {
    checkmate::assert_matrix(ypred, nrows = S, ncols = n)
    args <- list(y, ypred)
  } else if (
    measure %in%
      c(
        "elpd",
        "logscore"
      )
  ) {
    checkmate::assert_matrix(ylp, nrows = S, ncols = n)
    args <- list(ylp)
  }

  # provide equal weights or normalize given weights
  if (is.null(log_weights) || missing(log_weights)) {
    log_weights <- matrix(-log(S), nrow = S, ncol = n)
  } else {
    checkmate::assert_matrix(log_weights, nrows = S, ncols = n)
    log_weights <- sweep(
      log_weights,
      2,
      matrixStats::colLogSumExps(log_weights)
    )
  }

  do.call(pred_fun, append(args, list(log_weights)))
}

# ----------------------------- Measures -----------------------------

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
    "logscore" = .logscore_summary,
    "mlpd" = .logscore_summary,
    "r2" = .r2_summary,
    "mae" = .mae_summary,
    "rmse" = .rmse_summary,
    "mse" = .mse_summary,
    "acc" = .accuracy_summary,
    "balanced_acc" = .balanced_accuracy_summary,
    "rps" = .rps_summary,
    "crps" = .rps_summary,
    # TODO: is srps option needed
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
    "r2" = .pointwise_squared_error,
    "mae" = .pointwise_absolute_error,
    "rmse" = .pointwise_squared_error,
    "mse" = .pointwise_squared_error,
    "acc" = .pointwise_accuracy,
    "balanced_acc" = .pointwise_accuracy,
    "rps" = .pointwise_rps,
    "crps" = .pointwise_rps,
    "scrps" = .pointwise_rps
  )
}

#' Shared parameters for pointwise functions
#'
#' @param y scalar, leave one out value
#' @param ypred vector (S) of posterior predictive draws
#' @param ylp vector (S) of pointwise LOO log predictive densities
#' @param mupred vector (S) of point predictions
#' @param log_weights vector of loo weights (S) on the log scale
#'
#' @keywords internal
#' @name pointwise_measure_params
NULL

# Note that `pointwise` below is only for mse, rmse, and r2

#' Shared parameters for summary functions
#'
#' @param y vector of observed values (n)
#' @param ypred matrix of posterior draws (S x n) of posterior predictive draws
#' @param ylp matrix of posterior draws (S x n) of pointwise LOO log predictive densities
#' @param mupred matrix of posterior draws (S x n) of point predictions
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @param pointwise optional precomputed pointwise squared errors (n)
#'
#' @keywords internal
#' @name summary_measure_params
NULL

# ----------------------------- Metrics -----------------------------

#' Pointwise absolute error
#'
#' @noRd
#' @inheritParams pointwise_measure_params
.pointwise_absolute_error <- function(y, mupred, log_weights) {
  abs(y - .loo_weighted_mean(mupred, log_weights))
}

#' Pointwise squared error
#'
#' @noRd
#' @inheritParams pointwise_measure_params
.pointwise_squared_error <- function(y, mupred, log_weights) {
  (y - .loo_weighted_mean(mupred, log_weights))^2
}

#' Mean absolute error
#'
#' @noRd
#' @inheritParams summary_measure_params
.mae_summary <- function(y, mupred, log_weights) {
  .simple_pointwise_summary(vapply(
    seq_len(length(y)),
    function(i) {
      .pointwise_absolute_error(
        y[i],
        mupred[, i],
        log_weights[, i]
      )
    },
    numeric(1)
  ))
}

# TODO: take pointwise values from loo object--wrapper function shouldn't expose `pointwise`

#' Mean squared error
#'
#' @noRd
#' @inheritParams summary_measure_params
.mse_summary <- function(
  y,
  mupred,
  log_weights,
  pointwise = vapply(
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
) {
  .simple_pointwise_summary(pointwise)
}

#' Root mean squared error
#'
#' @noRd
#' @inheritParams summary_measure_params
.rmse_summary <- function(
  y,
  mupred,
  log_weights,
  pointwise = vapply(
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
) {
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
.r2_summary <- function(
  y,
  mupred,
  log_weights,
  pointwise = vapply(
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
) {
  n <- length(y)

  mse_loo_pointwise <- pointwise
  mse_loo <- mean(mse_loo_pointwise)
  se_mse_loo <- .se_helper(mse_loo_pointwise, mse_loo, n)

  squared_error_y_pointwise <- (y - mean(y))^2
  mse_y <- mean(squared_error_y_pointwise)

  se_r2 <- sqrt(
    se_mse_loo^2 -
      2 *
        (mse_loo / mse_y) *
        cov(mse_loo_pointwise, squared_error_y_pointwise) /
        n +
      (mse_loo^2 / mse_y^2) * var(squared_error_y_pointwise) / n
  ) /
    mse_y

  list(
    estimate = 1 - mse_loo / mse_y,
    se = se_r2,
    pointwise = mse_loo_pointwise
  )
}

#' Classification accuracy
#'
#' Assumes values in `mupred` only take on 0 or 1
#'
#' @noRd
#' @inheritParams pointwise_measure_params
.pointwise_accuracy <- function(y, mupred, log_weights) {
  .loo_weighted_mean(mupred == y, log_weights)
}

#' Classification accuracy
#'
#' Assumes values in `mupred` only take on 0 or 1
#'
#' @noRd
#' @inheritParams summary_measure_params
.accuracy_summary <- function(y, mupred, log_weights) {
  checkmate::assert_subset(mupred, choices = c(0, 1))
  .simple_pointwise_summary(vapply(
    seq_len(length(y)),
    function(i) {
      .pointwise_accuracy(
        y[i],
        mupred[, i],
        log_weights[, i]
      )
    },
    numeric(1)
  ))
}

#' Balanced classification accuracy
#'
#' @noRd
#' @inheritParams summary_measure_params
.balanced_accuracy_summary <- function(y, mupred, log_weights) {
  # TODO: check implementation again
  n <- length(y)
  cls_counts <- table(y)

  .simple_pointwise_summary(
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
    ) *
      n /
      (length(cls_counts) * as.numeric(cls_counts[match(y, names(cls_counts))]))
  )
}

# ----------------------------- Scores -------------------------------

# https://github.com/stan-dev/loo/blob/master/R/E_loo.R#L260
# https://github.com/stan-dev/loo/blob/master/R/helpers.R#L36
# https://github.com/stan-dev/loo/blob/master/R/loo.R#L427

# TODO: .energy multivariate CRPS, add later

#' Pointwise log score
#'
#' @noRd
#' @inheritParams pointwise_measure_params
.pointwise_elpd <- function(ylp, log_weights) {
  if (is.null(log_weights)) {
    matrixStats::logSumExp(ylp) - log(length(ylp))
  } else {
    matrixStats::logSumExp(ylp + log_weights)
  }
}

#' Log score
#'
#' @noRd
#' @inheritParams summary_measure_params
.logscore_summary <- function(ylp, log_weights) {
  # TODO: rewrite to remove hard dep on R > 4.1.0
  .elpd_summary(ylp, log_weights) |>
    (\(l) {
      n <- ncol(ylp)
      modifyList(
        l,
        list(estimate = l$estimate / n, se = l$se / n)
      )
    })()
}

#' Expected log-predictive density
#'
#' @noRd
#' @inheritParams summary_measure_params
.elpd_summary <- function(ylp, log_weights) {
  n <- ncol(ylp)
  pointwise <- vapply(
    seq_len(n),
    function(i) {
      .pointwise_elpd(ylp[, i], log_weights[, i])
    },
    numeric(1)
  )
  list(
    estimate = sum(pointwise),
    se = n * .se_helper(pointwise, mean(pointwise), n),
    pointwise = pointwise
  )
}

#' Pointwise CRPS
#'
#' @noRd
#' @param scaled logical. If true, computes SRPS/SCRPS
#' @inheritParams pointwise_measure_params
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

#' (Continuous) Ranked Probability Score
#'
#' Given draws ypred from the predictive distribution and the observations y,
#' computes
#'   - rank probability score (Epstein, 1969) for discrete ypred and y and
#'     continuous rank probability score (Matheson and Winkler, 1976;
#'     Gneiting & Raftery, 2007) for continuous ypred and y, using the
#'     probability weighted moment form (Taillardat et al., 2016;
#'     Zamo & Naveau, 2017)
#'   - scaled versions of these, if `scaled=TRUE` (Bolin & Wallin, 2023).
#'
#' Utility version of the score is returned, that is, bigger is
#' better, to match the utility version of log score / elpd (the original
#' rank probability score by Epstein (1969) was also in this direction).
#'
#' The same sample based $L$-moment estimator is used for continuous and
#' discrete variables. It is commonly stated that the probability
#' weighted moment form assumes F(x) is continuous. However, Hosking
#' (1990) states that $L$-moments can be used with discrete
#' distributions ``provided that the quantile function is `normalized'
#' in the sense of Widder (1941).''  Hosking (1996) states the same
#' condition more simply as ``A discrete random variable can be
#' approximated arbitrarily closely by a continuous random variable,
#' so the result is also valid for discrete random variables''.
#'
#' @references
#' \itemize{
#'   \item Bolin, D. and Wallin, J. (2023). Local scale invariance and
#'   robustness of proper scoring rules. \emph{Statistical Science},
#'   38(1):140-159.
#'
#'   \item Epstein, E.S. (1969). A scoring system for probability
#'   forecasts of ranked categories. \emph{Journal of Applied Meteorology},
#'   8(6):985-987.
#'
#'   \item Gneiting, T. and Raftery, A.E. (2007). Strictly Proper
#'   Scoring Rules, Prediction, and Estimation. \emph{Journal of the
#'   American Statistical Association}, 102(477):359-378.
#'
#'   \item Hosking, J.R.M. (1990). $L$-moments: analysis and estimation
#'   of distributions using linear combinations of order statistics.
#'   \emph{Journal of the Royal Statistical Society Series B: Statistical
#'   Methodology}, 52(1):105-124.
#'
#'   \item Hosking, J.R.M. (1996). Some theoretical results concerning
#'   $L$-moments. Research report RC 14492. IBM Thomas J. Watson Research
#'   Division.
#'
#'   \item Matheson, J.E., and Winkler, R.L. (1976). Scoring Rules for
#'   Continuous Probability Distributions. \emph{Management Science},
#'   22(10), 1087-1096.
#'
#'   \item Taillardat, M., Mestre, O., Zamo, M., and Naveau, P. (2016).
#'   Calibrated Ensemble Forecasts Using Quantile Regression Forests and
#'   Ensemble Model Output Statistics. \emph{Monthly Weather Review},
#'   144(6), 2375-2393.
#'
#'   \item Widder, D.V. (1941). \emph{The Laplace Transform}. Princeton:
#'   Princeton University Press.
#'
#'   \item Zamo, M., and Naveau, P. (2018). Estimation of the Continuous
#'   Ranked Probability Score with Limited Information and Applications
#'   to Ensemble Weather Forecasts. \emph{Mathematical Geosciences},
#'   50, 209–234.
#' }
#'
#' @noRd
#' @param scaled logical. If true, computes SRPS/SCRPS
#' @inheritParams summary_measure_params
.rps_summary <- function(y, ypred, log_weights, scaled = FALSE) {
  n <- length(y)

  pointwise <- vapply(
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

  .simple_pointwise_summary(pointwise)
}

#' @noRd
#' @inheritParams summary_measure_params
.srps_summary <- function(y, ypred, log_weights) {
  .rps_summary(y, ypred, log_weights, scaled = TRUE)
}

# ----------------------------- Helpers -----------------------------

# TODO: weights always provided
#' A wrapper around `stats::weighted.mean` which treats `NULL` weights as missing.
#'
#' @noRd
#' @param x vector to take mean of
#' @param log_weights log weights
#' @return weighted mean of `x`
.loo_weighted_mean <- function(x, log_weights) {
  sum(x * exp(log_weights - matrixStats::logSumExp(log_weights)))
}

.se_helper <- function(x, x_mean, n) {
  sqrt(sum((x - x_mean)^2) / (n * (n - 1)))
}

.simple_pointwise_summary <- function(pointwise) {
  est <- mean(pointwise)
  list(
    estimate = est,
    se = .se_helper(pointwise, est, length(pointwise)),
    pointwise = pointwise
  )
}
