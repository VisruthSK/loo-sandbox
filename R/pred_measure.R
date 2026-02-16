#' In-Sample Predictive Measure
#'
#' @param y Placeholder
#' @param ypred Placeholder
#' @param mupred Placeholder
#' @param ylp Placeholder
#' @param predperf Placeholder
#' @param measure Placeholder
#' @param group_ids Placeholder
#'
#' @return Placeholder
#'
#' @export
insample_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  predperf = NULL,
  measure = .pred_measure_choices(),
  group_ids = NULL,
  model_name = NULL
) {
  .pred_measure_engine(
    source = "insample",
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_insample = NULL,
    measure = measure,
    predperf = predperf,
    fold_id = NULL,
    group_ids = group_ids,
    model_name = model_name,
    psis_object = NULL,
    save_psis = FALSE
  )
}

#' K-fold Predictive Measure
#'
#' @param y Placeholder
#' @param ypred Placeholder
#' @param mupred Placeholder
#' @param ylp Placeholder
#' @param ylp_insample Placeholder
#' @param fold_id Placeholder
#' @param predperf Placeholder
#' @param measure Placeholder
#' @param group_ids Placeholder
#'
#' @return Placeholder
#'
#' @export
kfold_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_insample = NULL,
  fold_id = NULL,
  predperf = NULL,
  measure = .pred_measure_choices(),
  group_ids = NULL,
  model_name = NULL
) {
  .pred_measure_engine(
    source = "kfold",
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_insample = ylp_insample,
    measure = measure,
    predperf = predperf,
    fold_id = fold_id,
    group_ids = group_ids,
    model_name = model_name,
    psis_object = NULL,
    save_psis = FALSE
  )
}

#' Test-Set Predictive Measure
#'
#' @param y Placeholder
#' @param ypred Placeholder
#' @param mupred Placeholder
#' @param ylp Placeholder
#' @param ylp_insample Placeholder
#' @param predperf Placeholder
#' @param measure Placeholder
#' @param group_ids Placeholder
#'
#' @return Placeholder
#'
#' @export
test_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_insample = NULL,
  predperf = NULL,
  measure = .pred_measure_choices(),
  group_ids = NULL,
  model_name = NULL
) {
  .pred_measure_engine(
    source = "test",
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_insample = ylp_insample,
    measure = measure,
    predperf = predperf,
    fold_id = NULL,
    group_ids = group_ids,
    model_name = model_name,
    psis_object = NULL,
    save_psis = FALSE
  )
}

#' @export
print.pred_measure <- function(x, digits = 1, ...) {
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
print.insample_pred_measure <- function(x, digits = 1, ...) {
  print.pred_measure(x, digits = digits, ...)
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
print.test_pred_measure <- function(x, digits = 1, ...) {
  print.pred_measure(x, digits = digits, ...)
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
