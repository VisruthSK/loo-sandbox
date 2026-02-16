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
