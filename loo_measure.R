# TODO: write documentation
loo_pred_measure <- function(x, ...) {
  UseMethod("loo_pred_measure")
}

# TODO: use `{checkmate}`
loo_pred_measure.matrix <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  loo = NULL,
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
    "energy",
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
    stopifnot(is.matrix(ypred))
    args <- list(y, ypred)
  } else if (measure %in% c("elpd", "logscore")) {
    stopifnot(is.numeric(ylp))
    args <- list(y, ylp)
  }

  do.call(pred_fun, args)
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

#' @param y A vector of observed values
#' @param yhat A matrix of posterior draws (S x n)
#' @param weights optional weights for calculation of metric. Set to NULL if unweighted
#'
#' @keywords internal
#' @name .metric_common_params
NULL

# TODO: change name
#` Helper for .mae and .mse
.l12norm <- function(y, yhat, weights, penaltyFunc) {
  n <- length(y)
  stopifnot(
    is.matrix(yhat),
    is.function(penaltyFunc),
    ncol(yhat) == n
  )
  mu <- matrixStats::colWeightedMeans(yhat, w = weights)
  pointwise <- penaltyFunc(y - mu)
  est <- mean(pointwise)
  se <- sqrt(sum((pointwise - est)^2) / (n * (n - 1)))
  list(estimate = est, se = se, pointwise = pointwise)
}

#' Mean absolute error
#'
#' @noRd
#' @inheritParams .metric_common_params
.mae <- function(y, yhat, weights) {
  .l12norm(y, yhat, weights, abs)
}

#' Mean squared error
#'
#' @noRd
#' @inheritParams .metric_common_params
.mse <- function(y, yhat, weights) {
  .l12norm(y, yhat, weights, function(x) x^2)
}
# TODO: maybe memoize?

#' Root mean squared error
#'
#' @noRd
#' @inheritParams .metric_common_params
.rmse <- function(y, yhat, weights) {
  est <- .mse(y, yhat, weights)
  mean_mse <- est$estimate
  var_mse <- est$se^2
  var_rmse <- var_mse / mean_mse / 4 # Comes from the first order Taylor approx.
  list(
    estimate = sqrt(mean_mse),
    se = sqrt(var_rmse),
    pointwise = est$pointwise
  )
}

# TODO: double check stopifnots to see if ncol(yhat) is being compared to length of y instead of length of yhat

#' R^2
#'
#' @noRd
#' @inheritParams .metric_common_params
.r2 <- function(y, yhat, weights) {
  n <- length(y)
  stopifnot(is.matrix(yhat), ncol(yhat) == n)

  mse_loo_res <- .mse(y, yhat, weights)
  mse_loo <- mse_loo_res$estimate
  se_mse_loo <- mse_loo_res$se
  mse_loo_pointwise <- mse_loo_res$pointwise

  ybar <- mean(y)
  mse_y_res <- .mse(y, rep(ybar, n))
  mse_y <- mse_y_res$estimate
  mse_y_pointwise <- mse_y_res$pointwise

  # TODO: Could you please check this Aki
  se_r2 <- sqrt(
    se_mse_loo^2 -
      2 *
        (mse_loo / mse_y) *
        cov(mse_loo_pointwise, mse_y_pointwise) /
        n +
      (mse_loo^2 / mse_y^2) * var(mse_y_pointwise) / n
  ) /
    mse_y

  list(
    estimate = 1 - mse_loo / mse_y,
    se = se_r2,
    pointwise = mse_loo_pointwise
  )
}


# TODO: Ask Aki in slack
#' Classification accuracy
#'
#' Assuming values in `yhat` only take on 0 or 1
#'
#' @noRd
#' @inheritParams .metric_common_params
.accuracy <- function(y, yhat, weights) {
  # TODO: check yhat vals are 0, 1
  n <- length(y)
  stopifnot(is.matrix(yhat), ncol(yhat) == n)
  pointwise <- vapply(
    seq_len(n),
    function(j) weighted.mean(yhat[, j] == y[j], weights), # TODO: check if w = NULL is same as missing
    numeric(1)
  ) # TODO: make this weighted--maybe just weightedMean?
  est <- mean(pointwise)
  se <- sqrt(sum((pointwise - est)^2) / (n * (n - 1)))
  list(estimate = est, se = se, pointwise = pointwise)
}

# TODO: Check with Aki?
#' Balanced classification accuracy
#'
#' @noRd
#' @inheritParams .metric_common_params
.balanced_accuracy <- function(y, yhat, weights) {
  # TODO: add weights
  n <- length(y)
  stopifnot(is.matrix(yhat), ncol(yhat) == n)
  r <- vapply(seq_len(n), function(j) mean(yhat[, j] == y[j]), numeric(1))
  classes <- unique(y)
  recalls <- vapply(classes, function(cl) mean(r[y == cl]), numeric(1))
  names(recalls) <- as.character(classes)
  C <- length(recalls)
  est <- mean(recalls)
  se <- if (C > 1) sqrt(sum((recalls - est)^2) / (C * (C - 1))) else 0
  list(estimate = est, se = se, pointwise = recalls)
}

# ----------------------------- Scores -------------------------------

# TODO: .energy multivariate CRPS, add later

#' Log score
#'
#' @noRd
#' @param ylp Numeric vector of pointwise LOO log predictive densities.
.logscore <- function(y, ylp, weights) {
  n <- length(y)
  stopifnot(is.numeric(ylp), length(ylp) == n)
  est <- weighted.mean(ylp, weights) # TODO: check w
  se <- sqrt(n / (n - 1) * sum((ylp - est)^2)) # TODO: double check
  list(estimate = est, se = se, pointwise = ylp)
}

#' Expected log-predictive density
#'
#' @noRd
#' @param ylp Numeric vector of pointwise LOO log predictive densities.
.elpd <- function(y, ylp, weights) {
  n <- length(y)
  stopifnot(is.numeric(ylp), length(ylp) == n)
  est <- sum(ylp)
  se <- sqrt(n / (n - 1) * sum((ylp - mean(ylp))^2))
  list(estimate = est, se = se, pointwise = ylp)
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
#' @param w Optional nonnegative weights for draws
#' @param scaled logical. If true, computes SRPS/SCRPS
.rps <- function(y, yhat, w = NULL, scaled = FALSE) {
  # TODO: move this up to main function
  n <- length(y)
  stopifnot(ncol(yhat) == n)
  # TODO: change w to weight
  # TODO: no defaults

  pointwise <- vapply(
    seq_len(n),
    FUN = function(i) {
      x <- yhat[, i]
      obs <- y[i]

      if (is.null(w)) {
        EXy <- mean(abs(x - obs))
        x <- sort(x)
        n_draws <- length(x)
        EXX <- -2 * mean(x - 2 * x * (0:(n_draws - 1)) / (n_draws - 1))
      } else {
        EXy <- sum(w * abs(x - obs))
        ord <- order(x)
        x <- x[ord]
        w <- w[ord]
        cw <- cumsum(w)
        cw <- (cw - w) / (1 - w)
        EXX <- -2 * sum(w * (x - 2 * x * cw))
      }

      if (!scaled) {
        # Gneiting & Raftery (2007)
        rps <- -EXy + 0.5 * EXX
      } else {
        # Scaled version by Bolin & Wallin (2023)
        rps <- -EXy / EXX - 0.5 * log(EXX)
      }
      rps
    },
    FUN.VALUE = numeric(1)
  )

  list(
    estimate = mean(pointwise),
    se = sqrt(n / (n - 1) * sum((pointwise - mean(pointwise))^2)),
    pointwise = pointwise
  )
}
