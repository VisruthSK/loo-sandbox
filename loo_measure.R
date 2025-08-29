library(checkmate)

# TODO: write documentation
loo_pred_measure <- function(x, ...) {
  UseMethod("loo_pred_measure")
}

loo_pred_measure.matrix <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  loo = NULL,
  weights = NULL,
  measure = c(
    "logscore",
    "elpd",
    "r2",
    "rps",
    "crps",
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
  stopifnot(
    is.numeric(y),
    # TODO: check that for scores there are y and ypred and for metrics there are y and mupred
    (is.numeric(y) || is.function(ypred) && is.null(mupred)) ||
      (is.numeric(mupred) || is.function(mupred) && is.null(ypred))
  ) # TODO: flesh out checks
  measure <- match.arg(measure)
  pred_fun <- .loo_predictive_measure_fun(measure)

  if (
    measure %in%
      c(
        "mae",
        "mse",
        "rmse",
        "r2",
        "acc",
        "balanced_acc",
        "rps",
        "crps",
        "scrps"
        # ,"energy"
      )
  ) {
    assert_matrix(ypred, ncols = length(y))
    args <- list(y, ypred)
  } else if (measure %in% c("elpd", "logscore")) {
    assert_matrix(ylp, ncols = length(y))
    args <- list(y, ylp)
  }

  # TODO: check if weights is correct length or null
  assert_numeric(weights, len = length(y))

  do.call(pred_fun, append(args, weights))
}

# ----------------------------- Metrics -----------------------------

#' Select predictive measure function based on user's `measure` argument
#'
#' @noRd
#' @param measure The measure used.
#' @return The function used to compute predictive error or accuracy specified
#' by the argument `measure`.
.loo_predictive_measure_fun <- function(measure) {
  switch(
    measure,
    "elpd" = .elpd,
    "logscore" = .logscore,
    "mlpd" = .logscore,
    "r2" = .r2,
    "mae" = .mae,
    "rmse" = .rmse,
    "mse" = .mse,
    "acc" = .accuracy,
    "balanced_acc" = .balanced_accuracy,
    "rps" = .rps,
    "crps" = .rps,
    "scrps" = function(y, yhat) .rps(y, yhat, scaled = TRUE)
    # , "energy" = .energy
  )
}

#' @param y A scalar, leave one out value
#' @param yhat A vector of posterior draws of length S
#' @param loo_weights optional loo weights for calculation of metric. Set to NULL if unweighted
#'
#' @keywords internal
#' @name .metric_common_params
NULL

#' @param y A vector of observed values
#' @param yhat A matrix of posterior draws (S x n)
#' @param loo_weights optional weights for calculation of metric. Set to NULL if unweighted
#'
#' @keywords internal
#' @name .summary_metric_common_params
NULL

#' Mean absolute error
#'
#' @noRd
#' @inheritParams .metric_common_params
.mae <- function(y, yhat, loo_weights) {
  abs(y - .loo_weighted_mean(yhat, loo_weights))
}

#' Mean absolute error
#'
#' @noRd
#' @inheritParams .summary_metric_common_params
.mae_summary <- function(y, yhat, loo_weights) {
  .simple_pointwise_summary(vapply(
    seq_len(length(y)),
    function(i) .mae(y[i], yhat[, i], loo_weights),
    numeric(1)
  ))
}

#' Mean squared error
#'
#' @noRd
#' @inheritParams .metric_common_params
.mse <- function(y, yhat, loo_weights) {
  (y - .loo_weighted_mean(yhat, loo_weights))^2
}

#' Mean squared error
#'
#' @noRd
#' @inheritParams .summary_metric_common_params
.mse_summary <- function(y, yhat, loo_weights) {
  .simple_pointwise_summary(vapply(
    seq_len(length(y)),
    function(i) .mse(y[i], yhat[, i], loo_weights),
    numeric(1)
  ))
}

# TODO: maybe memoize?
# TODO: maybe this should be in a closure to share env with .rmse and .r2 so `.mse(y, yhat, weights)` is only calculated once on those?
# do.call() and pass around env?

#' Root mean squared error
#'
#' @noRd
#' @inheritParams .metric_common_params
.rmse <- function(y, yhat, loo_weights) {
  .mse(y, yhat, loo_weights)
}

#' Root mean squared error
#'
#' @noRd
#' @inheritParams .summary_metric_common_params
.rmse_summary <- function(y, yhat, loo_weights) {
  n <- length(y)
  sq <- vapply(
    seq_len(n),
    function(i) .rmse(y[i], yhat[, i], loo_weights),
    numeric(1)
  )^2
  mean_mse <- mean(sq)
  list(
    estimate = sqrt(mean_mse),
    se = sqrt(.se_helper(sq, mean_mse, n)^2 / mean_mse / 4), # Comes from the first order Taylor approx.
    pointwise = sq
  )
}

#' R^2
#'
#' @noRd
#' @inheritParams .metric_common_params
.r2 <- function(y, yhat, loo_weights) {
  .mse(y, yhat, loo_weights)
}

#' R^2
#'
#' @noRd
#' @inheritParams .summary_metric_common_params
.r2_summary <- function(y, yhat, loo_weights) {
  n <- length(y)

  mse_loo_pointwise <- vapply(
    seq_len(n),
    function(i) .r2(y[i], yhat[, i], loo_weights),
    numeric(1)
  )
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
#' @noRd
#' @inheritParams .metric_common_params
.accuracy <- function(y, yhat, loo_weights) {
  .loo_weighted_mean(yhat == y, loo_weights)
}

#' Classification accuracy
#'
#' Assuming values in `yhat` only take on 0 or 1
#'
#' @noRd
#' @inheritParams .summary_metric_common_params
.accuracy_summary <- function(y, yhat, loo_weights) {
  assert_subset(yhat, choices = c(0, 1))
  .simple_pointwise_summary(vapply(
    seq_len(length(y)),
    function(i) .accuracy(y[i], yhat[, i], loo_weights),
    numeric(1)
  ))
}

# TODO: write summaries

#' Balanced classification accuracy
#'
#' @noRd
#' @inheritParams .metric_common_params
.balanced_accuracy <- function(y, yhat, loo_weights) {
  .accuracy(y, yhat, loo_weights)
}

#' Balanced classification accuracy
#'
#' @noRd
#' @inheritParams .summary_metric_common_params
.balanced_accuracy_summary <- function(y, yhat, loo_weights) {
  n <- length(y)
  cls_counts <- table(y)

  .simple_pointwise_summary(
    vapply(
      seq_len(n),
      function(i) .balanced_accuracy(y[i], yhat[, i], loo_weights),
      numeric(1)
    ) *
      n /
      (length(cls_counts) * as.numeric(cls_counts[match(y, names(cls_counts))]))
  )
}

# ----------------------------- Scores -------------------------------

# TODO: .energy multivariate CRPS, add later

# TODO: ylp is a matrix
# TODO: Do pointwise and summary
#' Log score
#'
#' @noRd
#' @param ylp Numeric vector of pointwise LOO log predictive densities.
.logscore <- function(y, ylp, weights) {
  w <- if (is.null(weights)) 1 else weights
  .elpd(y, ylp * w, NULL) |>
    (\(l) {
      d <- if (is.null(weights)) length(y) else sum(weights)
      modifyList(
        l,
        list(estimate = l$estimate / d, se = l$se / d)
      )
    })()
}

#' Expected log-predictive density
#'
#' @noRd
#' @param ylp Numeric vector of pointwise LOO log predictive densities.
.elpd <- function(y, ylp, weights) {
  n <- length(y)
  list(
    estimate = sum(ylp),
    se = sqrt(n / (n - 1) * sum((ylp - .loo_weighted_mean(ylp, weights))^2)),
    pointwise = ylp
  )
}

# TODO: change measures to compute pointwise; summarize elsewhere
# https://github.com/stan-dev/loo/blob/master/R/E_loo.R#L260
# https://github.com/stan-dev/loo/blob/master/R/helpers.R#L36
# https://github.com/stan-dev/loo/blob/master/R/loo.R#L427

#' CRPS
#'
#' @inheritParams .metric_common_params
.rps <- function(y, yhat, loo_weights, scaled) {
  if (is.null(loo_weights)) {
    EXy <- mean(abs(y - yhat))
    y <- sort(y)
    n <- length(y)
    EXX <- -2 * mean(y - 2 * y * (0:(n - 1)) / (n - 1))
  } else {
    EXy <- sum(loo_weights * abs(y - yhat))
    ord <- order(y)
    y <- y[ord]
    loo_weights <- loo_weights[ord]
    cw <- (cumsum(loo_weights) - loo_weights) / (1 - loo_weights)
    EXX <- -2 * sum(loo_weights * (y - 2 * y * cw))
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
#' Given draws x from the predictive distribution and the observations y,
#' computes
#'   - rank probability score (Epstein, 1969) for discrete x and y and
#'     continuous rank probability score (Matheson and Winkler, 1976;
#'     Gneiting & Raftery, 2007) for continuous x and y, using the
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
#' @noRd
#' @param y A vector of observed values
#' @param ypred Predictive draws matrix
#' @param loo_weights Optional nonnegative loo_weightss for draws
#' @param scaled logical. If true, computes SRPS/SCRPS
.rps <- function(y, yhat, loo_weights, scaled) {
  n <- length(y)

  pointwise <- vapply(
    seq_len(n),
    FUN = function(i) .rps(y[i], yhat[, i], loo_weights, scaled),
    FUN.VALUE = numeric(1)
  )

  list(
    estimate = mean(pointwise),
    se = .se_helper(pointwise, mean(pointwise), n),
    pointwise = pointwise
  )
}

# ----------------------------- Helpers -----------------------------

#' A wrapper around `stats::weighted.mean` which treats `NULL` weights as missing.
#'
#' @noRd
#' @param x vector to take mean of
#' @param weights optional weights. Set to `NULL` if unweighted mean is desired.
#' @return The mean of `x`, optionally weighted by `weights` if specified.
.loo_weighted_mean <- function(x, weights) {
  if (missing(weights) || is.null(weights)) {
    mean(x)
  }
  weighted.mean(x, weights)
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
