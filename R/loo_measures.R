#' Mean squared error
#'
#' @inheritParams summary_measure_params
#'
#' @export
mse <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .mse_summary(
    y = y,
    mupred = mupred,
    log_weights = log_weights,
    pointwise = NULL
  )
}

#' Root mean squared error
#'
#' @inheritParams summary_measure_params
#'
#' @export
rmse <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .rmse_summary(
    y = y,
    mupred = mupred,
    log_weights = log_weights,
    pointwise = NULL
  )
}

#' Mean absolute error
#'
#' @inheritParams summary_measure_params
#'
#' @export
mae <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .mae_summary(y = y, mupred = mupred, log_weights = log_weights)
}

#' R^2
#'
#' @inheritParams summary_measure_params
#'
#' @export
r2 <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .r2_summary(
    y = y,
    mupred = mupred,
    log_weights = log_weights,
    pointwise = NULL
  )
}

#' Classification accuracy
#'
#' @inheritParams summary_measure_params
#'
#' @export
acc <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .accuracy_summary(
    y = y,
    mupred = mupred,
    log_weights = log_weights,
    pointwise = NULL
  )
}

#' Balanced classification accuracy
#'
#' @inheritParams summary_measure_params
#'
#' @export
balanced_acc <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .balanced_accuracy_summary(
    y = y,
    mupred = mupred,
    log_weights = log_weights,
    pointwise = NULL
  )
}

#' Expected log-predictive density
#'
#' @inheritParams summary_measure_params
#'
#' @export
elpd <- function(ylp, log_weights) {
  n <- ncol(ylp)
  S <- nrow(ylp)
  checkmate::assert_matrix(ylp, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .elpd_summary(ylp = ylp, log_weights = log_weights, pointwise = NULL)
}

#' Log score
#'
#' @inheritParams summary_measure_params
#'
#' @export
logscore <- function(ylp, log_weights) {
  n <- ncol(ylp)
  S <- nrow(ylp)
  checkmate::assert_matrix(ylp, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .logscore_summary(ylp = ylp, log_weights = log_weights, pointwise = NULL)
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
#' @inheritParams summary_measure_params
#'
#' @export
rps <- function(y, ypred, log_weights) {
  n <- length(y)
  S <- nrow(ypred)
  checkmate::assert_matrix(ypred, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .rps_summary(
    y = y,
    ypred = ypred,
    log_weights = log_weights,
    scaled = FALSE,
    pointwise = NULL
  )
}

#' Continuous Ranked Probability Score
#'
#' @rdname rps
#' @export
crps <- function(y, ypred, log_weights) {
  rps(y, ypred, log_weights)
}

#' Scaled Ranked Probability Score
#'
#' @rdname rps
#' @export
srps <- function(y, ypred, log_weights) {
  n <- length(y)
  S <- nrow(ypred)
  checkmate::assert_matrix(ypred, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .srps_summary(
    y = y,
    ypred = ypred,
    log_weights = log_weights,
    pointwise = NULL
  )
}

#' Scaled Continuous Ranked Probability Score
#'
#' @rdname rps
#' @export
scrps <- function(y, ypred, log_weights) {
  srps(y, ypred, log_weights)
}
