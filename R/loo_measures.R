#' Mean squared error
#'
#' @param y vector of observed values (n)
#' @param mupred matrix of posterior draws (S x n) of point predictions
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
mse <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .mse_summary(y, mupred, log_weights, pointwise = NULL)
}

#' Root mean squared error
#'
#' @param y vector of observed values (n)
#' @param mupred matrix of posterior draws (S x n) of point predictions
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
rmse <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .rmse_summary(y, mupred, log_weights, pointwise = NULL)
}

#' Mean absolute error
#'
#' @param y vector of observed values (n)
#' @param mupred matrix of posterior draws (S x n) of point predictions
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
mae <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .mae_summary(y, mupred, log_weights)
}

#' R^2
#'
#' @param y vector of observed values (n)
#' @param mupred matrix of posterior draws (S x n) of point predictions
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
r2 <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .r2_summary(y, mupred, log_weights, pointwise = NULL)
}

#' Classification accuracy
#'
#' @param y vector of observed values (n)
#' @param mupred matrix of posterior draws (S x n) of point predictions
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
acc <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .accuracy_summary(y, mupred, log_weights, pointwise = NULL)
}

#' Balanced classification accuracy
#'
#' @param y vector of observed values (n)
#' @param mupred matrix of posterior draws (S x n) of point predictions
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
balanced_acc <- function(y, mupred, log_weights) {
  n <- length(y)
  S <- nrow(mupred)
  checkmate::assert_matrix(mupred, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .balanced_accuracy_summary(y, mupred, log_weights, pointwise = NULL)
}

#' Expected log-predictive density
#'
#' @param ylp matrix of posterior draws (S x n) of pointwise LOO log predictive densities
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
elpd <- function(ylp, log_weights) {
  n <- ncol(ylp)
  S <- nrow(ylp)
  checkmate::assert_matrix(ylp, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .elpd_summary(ylp, log_weights, pointwise = NULL)
}

#' Log score
#'
#' @param ylp matrix of posterior draws (S x n) of pointwise LOO log predictive densities
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
logscore <- function(ylp, log_weights) {
  n <- ncol(ylp)
  S <- nrow(ylp)
  checkmate::assert_matrix(ylp, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .logscore_summary(ylp, log_weights, pointwise = NULL)
}

#' (Continuous) Ranked Probability Score
#'
#' @param y vector of observed values (n)
#' @param ypred matrix of posterior draws (S x n) of posterior predictive draws
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
rps <- function(y, ypred, log_weights) {
  n <- length(y)
  S <- nrow(ypred)
  checkmate::assert_matrix(ypred, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .rps_summary(y, ypred, log_weights, scaled = FALSE, pointwise = NULL)
}

#' (Continuous) Ranked Probability Score
#'
#' @param y vector of observed values (n)
#' @param ypred matrix of posterior draws (S x n) of posterior predictive draws
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
crps <- function(y, ypred, log_weights) {
  rps(y, ypred, log_weights)
}

#' Scaled (Continuous) Ranked Probability Score
#'
#' @param y vector of observed values (n)
#' @param ypred matrix of posterior draws (S x n) of posterior predictive draws
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
srps <- function(y, ypred, log_weights) {
  n <- length(y)
  S <- nrow(ypred)
  checkmate::assert_matrix(ypred, nrows = S, ncols = n)
  checkmate::assert_matrix(log_weights, nrows = S, ncols = n)

  log_weights <- sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights)
  )

  .srps_summary(y, ypred, log_weights, pointwise = NULL)
}

#' Scaled (Continuous) Ranked Probability Score
#'
#' @param y vector of observed values (n)
#' @param ypred matrix of posterior draws (S x n) of posterior predictive draws
#' @param log_weights matrix of loo weights (S x n) on the log scale
#'
#' @export
scrps <- function(y, ypred, log_weights) {
  srps(y, ypred, log_weights)
}
