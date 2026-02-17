test_that("loo_pred_measure adds base loo summaries without loo object", {
  inputs <- roaches_models$model1
  result <- loo_pred_measure(
    y = inputs$y,
    mupred = inputs$mupred,
    ylp = inputs$ylp,
    psis_object = inputs$loo$psis_object,
    measure = "r2"
  )

  expect_true(
    all(c("elpd", "p_eff", "ic", "r2") %in% rownames(result$estimates))
  )
  expect_true(
    all(
      c("elpd", "p_eff", "ic", "squared_error") %in%
        colnames(result$pointwise)
    )
  )
})

test_that("loo_pred_measure aligns with loo summaries", {
  inputs <- roaches_models$model1
  loo_obj <- inputs$loo
  result <- loo_pred_measure(
    ylp = inputs$ylp,
    loo = loo_obj,
    measure = "elpd"
  )

  expect_equal(
    result$estimates[rownames(loo_obj$estimates), , drop = FALSE],
    loo_obj$estimates
  )
  expect_equal(
    result$pointwise[, colnames(loo_obj$pointwise), drop = FALSE],
    loo_obj$pointwise
  )
})

test_that("loo_pred_measure logscore alias uses loo elpd summaries", {
  inputs <- roaches_models$model1
  result <- loo_pred_measure(
    loo = inputs$loo,
    measure = "logscore"
  )
  expect_false("logscore" %in% rownames(result$estimates))
  expect_equal(
    result$estimates["elpd", , drop = FALSE],
    inputs$loo$estimates["elpd", , drop = FALSE]
  )
  expect_equal(result$pointwise[, "elpd"], inputs$loo$pointwise[, "elpd"])
})

test_that("loo_pred_measure rejects non-loo pred_measure objects", {
  inputs <- roaches_models$model1
  nonloo <- insample_pred_measure(
    y = inputs$y,
    ylp = inputs$ylp,
    measure = "elpd"
  )

  expect_error(
    loo_pred_measure(
      ylp = inputs$ylp,
      loo = nonloo,
      measure = "elpd"
    ),
    "`loo` must be a `loo` object"
  )
})

test_that("loo_pred_measure computes finite SRPS and keeps score aliases", {
  inputs <- roaches_models$model1
  result_rps <- loo_pred_measure(
    y = inputs$y,
    ypred = inputs$ypred,
    loo = inputs$loo,
    measure = "rps"
  )
  result_crps <- loo_pred_measure(
    y = inputs$y,
    ypred = inputs$ypred,
    loo = inputs$loo,
    measure = "crps"
  )
  result_srps <- loo_pred_measure(
    y = inputs$y,
    ypred = inputs$ypred,
    loo = inputs$loo,
    measure = "srps"
  )
  result_scrps <- loo_pred_measure(
    y = inputs$y,
    ypred = inputs$ypred,
    loo = inputs$loo,
    measure = "scrps"
  )

  pointwise_rps_ref <- vapply(
    seq_along(inputs$y),
    function(i) {
      log_weights <- result_rps$log_weights[, i]
      w <- exp(log_weights - matrixStats::logSumExp(log_weights))
      ypred <- inputs$ypred[, i]
      exy <- sum(w * abs(inputs$y[i] - ypred))
      ord <- order(ypred)
      ypred <- ypred[ord]
      w <- w[ord]
      cw <- (cumsum(w) - w) / (1 - w)
      exx <- -2 * sum(w * (ypred - 2 * ypred * cw))
      -exy + 0.5 * exx
    },
    numeric(1)
  )
  pointwise_srps_ref <- vapply(
    seq_along(inputs$y),
    function(i) {
      log_weights <- result_rps$log_weights[, i]
      w <- exp(log_weights - matrixStats::logSumExp(log_weights))
      ypred <- inputs$ypred[, i]
      exy <- sum(w * abs(inputs$y[i] - ypred))
      ord <- order(ypred)
      ypred <- ypred[ord]
      w <- w[ord]
      cw <- (cumsum(w) - w) / (1 - w)
      exx <- -2 * sum(w * (ypred - 2 * ypred * cw))
      -exy / exx - 0.5 * log(exx)
    },
    numeric(1)
  )

  expect_true(all(is.finite(result_srps$pointwise[, "srps"])))
  expect_true(is.finite(result_srps$estimates["srps", "Estimate"]))
  expect_true(is.finite(result_srps$estimates["srps", "SE"]))
  expect_equal(result_rps$pointwise[, "rps"], pointwise_rps_ref)
  expect_equal(result_srps$pointwise[, "srps"], pointwise_srps_ref)
  expect_equal(result_rps$pointwise[, "rps"], result_crps$pointwise[, "crps"])
  expect_equal(result_rps$estimates["rps", ], result_crps$estimates["crps", ])
  expect_equal(
    result_srps$pointwise[, "srps"],
    result_scrps$pointwise[, "scrps"]
  )
  expect_equal(
    result_srps$estimates["srps", ],
    result_scrps$estimates["scrps", ]
  )
})

test_that("loo_pred_measure keeps mlpd distinct and aliases logscore to elpd", {
  inputs <- roaches_models$model1
  result_elpd <- loo_pred_measure(
    ylp = inputs$ylp,
    loo = inputs$loo,
    measure = "elpd"
  )
  result_mlpd <- loo_pred_measure(
    ylp = inputs$ylp,
    loo = inputs$loo,
    measure = "mlpd"
  )
  result_logscore <- loo_pred_measure(
    ylp = inputs$ylp,
    loo = inputs$loo,
    measure = "logscore"
  )

  mlpd_row <- .match_results_column("mlpd")
  mlpd_col <- .match_pointwise_column("mlpd")
  logscore_row <- .match_results_column("logscore")
  logscore_col <- .match_pointwise_column("logscore")
  expect_true(mlpd_row %in% rownames(result_mlpd$estimates))
  expect_true(mlpd_col %in% colnames(result_mlpd$pointwise))
  expect_false("logscore" %in% rownames(result_logscore$estimates))
  expect_true(logscore_row %in% rownames(result_logscore$estimates))
  expect_true(logscore_col %in% colnames(result_logscore$pointwise))
  expect_equal(
    result_logscore$estimates[logscore_row, "Estimate"],
    result_elpd$estimates["elpd", "Estimate"]
  )
  expect_equal(
    result_logscore$estimates[logscore_row, "SE"],
    result_elpd$estimates["elpd", "SE"]
  )
  expect_equal(
    result_logscore$pointwise[, logscore_col],
    result_elpd$pointwise[, "elpd"]
  )
  expect_equal(
    result_mlpd$estimates[mlpd_row, "Estimate"],
    result_elpd$estimates["elpd", "Estimate"] /
      length(result_mlpd$pointwise[, mlpd_col])
  )
  expect_equal(
    result_mlpd$estimates[mlpd_row, "SE"],
    result_elpd$estimates["elpd", "SE"] /
      length(result_mlpd$pointwise[, mlpd_col])
  )
})

test_that("loo_pred_measure stores model_name metadata", {
  inputs <- roaches_models$model1
  result <- loo_pred_measure(
    ylp = inputs$ylp,
    loo = inputs$loo,
    measure = "elpd",
    model_name = "fit_lin"
  )
  expect_equal(result$metadata$model_name, "fit_lin")
})

test_that("loo_pred_measure fallback PSIS path matches loo::loo", {
  inputs <- roaches_models$model1
  ll <- inputs$ylp

  ref <- suppressWarnings(loo::loo(ll, r_eff = 1, save_psis = TRUE))
  ref_log_weights <- stats::weights(
    ref$psis_object,
    normalize = TRUE,
    log = TRUE
  )
  result <- suppressWarnings(
    loo_pred_measure(
      ylp = ll,
      measure = "elpd"
    )
  )

  expect_equal(result$log_weights, ref_log_weights)
  expect_equal(result$pointwise[, "elpd"], ref$pointwise[, "elpd_loo"])
  expect_equal(
    result$estimates["elpd", "Estimate"],
    ref$estimates["elpd_loo", "Estimate"]
  )
  expect_equal(
    result$estimates["elpd", "SE"],
    ref$estimates["elpd_loo", "SE"]
  )
})
