measure_specs <- function(temp) {
  list(
    list(measure = "r2", args = list(y = temp$y, mupred = temp$mupred)),
    list(measure = "rmse", args = list(y = temp$y, mupred = temp$mupred)),
    list(measure = "mse", args = list(y = temp$y, mupred = temp$mupred)),
    list(measure = "mae", args = list(y = temp$y, mupred = temp$mupred)),
    list(
      measure = "acc",
      args = list(y = temp$y_binary, mupred = temp$binary_mupred)
    ),
    list(
      measure = "balanced_acc",
      args = list(y = temp$y_binary, mupred = temp$binary_mupred)
    ),
    list(measure = "rps", args = list(y = temp$y, ypred = temp$ypred)),
    list(measure = "crps", args = list(y = temp$y, ypred = temp$ypred)),
    list(measure = "srps", args = list(y = temp$y, ypred = temp$ypred)),
    list(measure = "scrps", args = list(y = temp$y, ypred = temp$ypred)),
    list(measure = "elpd", args = list(ylp = temp$ylp)),
    list(measure = "mlpd", args = list(ylp = temp$ylp)),
    list(measure = "logscore", args = list(ylp = temp$ylp))
  )
}

run_measure_snapshots <- function(loo_start, measures) {
  loo_iter <- loo_start
  for (measure in measures) {
    loo_iter <- do.call(
      loo_pred_measure,
      c(
        measure$args,
        list(measure = measure$measure, loo = loo_iter)
      )
    )
    expect_snapshot_output(print(loo_iter))
  }
  loo_iter
}

test_that("loo_pred_measure print snapshots", {
  temp <- roaches_models$model1
  loo_ordered <- run_measure_snapshots(temp$loo, measure_specs(temp))
  loo_shuffled <- run_measure_snapshots(
    temp$loo,
    with(set.seed(0), sample(measure_specs(temp)))
  )
  expect_setequal(
    rownames(loo_ordered$estimates),
    rownames(loo_shuffled$estimates)
  )
  expect_equal(
    loo_ordered$estimates[rownames(loo_shuffled$estimates), , drop = FALSE],
    loo_shuffled$estimates
  )
})
