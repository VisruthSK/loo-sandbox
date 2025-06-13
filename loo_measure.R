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
  # TODO: write function
  stopifnot(
    is.numeric(y),
    # TODO: check that for scores there are y and ypred and for metrics there are y and mupred
    (is.numeric(ypred) || is.function(ypred) && is.null(mupred)) ||
      (is.numeric(mupred) || is.function(mupred) && is.null(ypred))
  ) # TODO: flesh out checks
  measure <- match.arg(measure)
  # is_metric <- measure %in% c("mae", "rmse", "mse", "acc", "balanced_acc", "r2")

  predictive_measure_func <- .loo_predictive_measure_fun(measure)
}


# ----------------------------- Internals -----------------------------

#' Select predictive measure function based on user's `measure` argument
#'
#' @noRd
#' @param measure The measure used.
#' @return The function used to compute predictive error or accuracy specified
#' by the argument `measure`.
.loo_predictive_measure_fun <- function(measure) {
  switch(
    measure,
    "r2" = .r2,
    "mae" = .mae,
    "rmse" = .rmse,
    "mse" = .mse,
    "acc" = .accuracy,
    'balanced_acc' = .balanced_accuracy
  )
}

#' Mean absolute error
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.mae <- function(y, yhat) {
  stopifnot(length(y) == length(yhat))
  n <- length(y)
  e <- abs(y - yhat)
  list(estimate = mean(e), se = sd(e) / sqrt(n))
}

#' Mean squared error
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.mse <- function(y, yhat) {
  stopifnot(length(y) == length(yhat))
  n <- length(y)
  e <- (y - yhat)^2
  list(estimate = mean(e), se = sd(e) / sqrt(n))
}

#' Root mean squared error
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.rmse <- function(y, yhat) {
  est <- .mse(y, yhat)
  mean_mse <- est$estimate
  var_mse <- est$se^2
  var_rmse <- var_mse / mean_mse / 4 # Comes from the first order Taylor approx.
  return(list(estimate = sqrt(mean_mse), se = sqrt(var_rmse)))
}

#' Classification accuracy
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.accuracy <- function(y, yhat) {
  stopifnot(
    length(y) == length(yhat),
    all(y <= 1 & y >= 0),
    all(yhat <= 1 & yhat >= 0)
  )
  n <- length(y)
  yhat <- as.integer(yhat > 0.5)
  acc <- as.integer(yhat == y)
  est <- mean(acc)
  list(estimate = est, se = sqrt(est * (1 - est) / n))
}

#' Balanced classification accuracy
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.balanced_accuracy <- function(y, yhat) {
  stopifnot(
    length(y) == length(yhat),
    all(y <= 1 & y >= 0),
    all(yhat <= 1 & yhat >= 0)
  )
  n <- length(y)
  yhat <- as.integer(yhat > 0.5)
  mask <- y == 0

  tn <- mean(yhat[mask] == y[mask]) # True negatives
  tp <- mean(yhat[!mask] == y[!mask]) # True positives

  bls_acc <- (tp + tn) / 2
  # This approximation has quite large bias for small samples
  bls_acc_var <- (tp * (1 - tp) + tn * (1 - tn)) / 4
  list(estimate = bls_acc, se = sqrt(bls_acc_var / n))
}

#' R^2
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.r2 <- function(y, yhat) {
  stopifnot(length(y) == length(yhat))
  1 - (sum((y - yhat)^2) / sum((y - mean(y))^2))
}
