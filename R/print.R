# print.R: S3 print methods for predictive measure objects and source
# labeling helpers.


#' @export
print.pred_measure <- function(x, digits = 1, ...) {
  # TODO: ppl should be able to choose what is printed
  dims <- attr(x, "dims")
  if (is.null(dims) && !is.null(x$log_weights)) {
    dims <- dim(x$log_weights)
  }
  source <- .pred_measure_source_label(x)

  cat("\n")
  if (!is.null(dims) && length(dims) == 2) {
    cat(
      sprintf(
        "Computed from %s posterior draws and %s observations.\n",
        dims[1],
        dims[2]
      )
    )
  }
  cat(sprintf("Data source: %s\n\n", source))
  print(
    format(round(as.data.frame(x$estimates), digits), nsmall = digits),
    quote = FALSE
  )
  invisible(x)
}

#' @export
print.kfold_pred_measure <- function(x, digits = 1, ...) {
  print.pred_measure(x, digits = digits, ...)
  if (!is.null(x$metadata$fold_id)) {
    cat(sprintf("Folds: %s\n", length(unique(x$metadata$fold_id))))
  }
  invisible(x)
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
  labels <- c("good", "bad", "very bad")
  ranges <- c("k <= 0.7", "0.7 < k <= 1", "k > 1")
  bins <- cut(
    pareto_k,
    breaks = c(-Inf, 0.7, 1, Inf),
    labels = labels,
    include.lowest = TRUE
  )
  counts <- tabulate(as.integer(bins[!is.na(bins)]), nbins = length(labels))

  cat("Pareto k diagnostic values:\n")
  for (i in seq_along(labels)) {
    cat(sprintf(
      "  %s (%s): %d (%.1f%%)\n",
      labels[i],
      ranges[i],
      counts[i],
      100 * counts[i] / n
    ))
  }

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

.pred_measure_source_label <- function(x) {
  cls <- class(x)
  if ("loo_pred_measure" %in% cls) {
    return("loo")
  }
  if ("insample_pred_measure" %in% cls) {
    return("in-sample")
  }
  if ("kfold_pred_measure" %in% cls) {
    return("k-fold")
  }
  if ("test_pred_measure" %in% cls) {
    return("test")
  }
  if (!is.null(x$metadata$source)) {
    return(as.character(x$metadata$source))
  }
  "unknown"
}
