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

test_that("loo_pred_measure keeps mlpd alias behavior", {
  inputs <- roaches_models$model1
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
  expect_true(mlpd_row %in% rownames(result_mlpd$estimates))
  expect_true(mlpd_col %in% colnames(result_mlpd$pointwise))
  expect_equal(
    result_mlpd$estimates[mlpd_row, "Estimate"],
    result_logscore$estimates["logscore", "Estimate"]
  )
  expect_equal(
    result_mlpd$estimates[mlpd_row, "SE"],
    result_logscore$estimates["logscore", "SE"]
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
