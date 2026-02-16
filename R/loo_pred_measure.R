# TODO: document function--look at loo #281
#' LOO Predictive Measure
#'
#' @param y Placeholder
#' @param ypred Placeholder
#' @param mupred Placeholder
#' @param ylp Placeholder
#' @param loo Placeholder
#' @param measure Placeholder
#' @param group_ids Placeholder
#' @param model_name Placeholder
#' @param psis_object Placeholder
#' @param save_psis Placeholder
#'
#' @return Placeholder
#'
#' @export
loo_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measure = .pred_measure_choices(),
  group_ids = NULL,
  model_name = NULL,
  loo = NULL, # This `loo` object needs to be fit with `save_psis = TRUE`
  psis_object = NULL,
  save_psis = FALSE
) {
  .pred_measure_engine(
    source = "loo",
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_insample = NULL,
    measure = measure,
    predperf = loo,
    fold_id = NULL,
    group_ids = group_ids,
    model_name = model_name,
    psis_object = psis_object,
    save_psis = save_psis
  )
}

#' @export
print.loo_pred_measure <- function(x, digits = 1, plot_k = FALSE, ...) {
  print.pred_measure(x, digits = digits, ...)
  cat("------\n")
  pareto_k <- x$diagnostics$pareto_k
  if (is.null(pareto_k)) {
    cat("No Pareto-k diagnostics available.\n")
    return(invisible(x))
  }
  n <- length(pareto_k)
  bins <- cut(
    pareto_k,
    breaks = c(-Inf, 0.7, 1, Inf),
    labels = c("good", "bad", "very bad"),
    include.lowest = TRUE
  )
  counts <- table(bins)
  count_good <- if ("good" %in% names(counts)) {
    as.integer(counts["good"])
  } else {
    0L
  }
  count_bad <- if ("bad" %in% names(counts)) as.integer(counts["bad"]) else 0L
  count_vbad <- if ("very bad" %in% names(counts)) {
    as.integer(counts["very bad"])
  } else {
    0L
  }
  cat("Pareto k diagnostic values:\n")
  cat(sprintf(
    "  good (k <= 0.7): %d (%.1f%%)\n",
    count_good,
    100 * count_good / n
  ))
  cat(sprintf(
    "  bad (0.7 < k <= 1): %d (%.1f%%)\n",
    count_bad,
    100 * count_bad / n
  ))
  cat(sprintf(
    "  very bad (k > 1): %d (%.1f%%)\n",
    count_vbad,
    100 * count_vbad / n
  ))
  if (plot_k) {
    graphics::plot(
      pareto_k,
      ylab = "Pareto-k",
      xlab = "Observation",
      main = "Pareto-k diagnostics",
      pch = 16
    )
  }
  invisible(x)
}
