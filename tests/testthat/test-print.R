test_that("loo_pred_measure print snapshots", {
  temp <- roaches_models$model1
  loo1 <- loo_pred_measure(
    temp$y,
    mupred = temp$mupred,
    measure = "r2",
    loo = temp$loo
  )
  loo2 <- loo_pred_measure(
    temp$y,
    mupred = temp$mupred,
    measure = "rmse",
    loo = loo1
  )

  expect_snapshot_output(print(loo1))
  expect_snapshot_output(print(loo2))
})
