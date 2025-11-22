# Helper to run the comparison loop
run_wrapper_test <- function(
  measure,
  wrapper_fn,
  type = c("default", "binary", "ylp", "ypred")
) {
  type <- match.arg(type)

  test_that(paste(measure, "wrapper matches loo_pred_measure"), {
    for (model_name in names(roaches_models)) {
      inputs <- roaches_models[[model_name]]
      log_weights <- weights(inputs$loo$psis_object)

      if (type == "default") {
        # mse, rmse, mae, r2
        expected <- loo_pred_measure(
          y = inputs$y,
          mupred = inputs$mupred,
          loo = inputs$loo,
          measure = measure
        )
        result <- wrapper_fn(
          y = inputs$y,
          mupred = inputs$mupred,
          log_weights = log_weights
        )
      } else if (type == "binary") {
        # acc, balanced_acc
        expected <- loo_pred_measure(
          y = inputs$y_binary,
          mupred = inputs$binary_mupred,
          loo = inputs$loo,
          measure = measure
        )
        result <- wrapper_fn(
          y = inputs$y_binary,
          mupred = inputs$binary_mupred,
          log_weights = log_weights
        )
      } else if (type == "ylp") {
        # elpd, logscore
        expected <- loo_pred_measure(
          y = inputs$y,
          ylp = inputs$ylp,
          loo = inputs$loo,
          measure = measure
        )
        result <- wrapper_fn(
          ylp = inputs$ylp,
          log_weights = log_weights
        )
      } else if (type == "ypred") {
        # rps, crps, srps, scrps
        expected <- loo_pred_measure(
          y = inputs$y,
          ypred = inputs$ypred,
          loo = inputs$loo,
          measure = measure
        )
        result <- wrapper_fn(
          y = inputs$y,
          ypred = inputs$ypred,
          log_weights = log_weights
        )
      }

      expect_equal(result, expected, info = model_name)
    }
  })
}

run_wrapper_test("mse", mse, "default")
run_wrapper_test("rmse", rmse, "default")
run_wrapper_test("mae", mae, "default")
run_wrapper_test("r2", r2, "default")

run_wrapper_test("acc", acc, "binary")
run_wrapper_test("balanced_acc", balanced_acc, "binary")

run_wrapper_test("elpd", elpd, "ylp")
run_wrapper_test("logscore", logscore, "ylp")

run_wrapper_test("rps", rps, "ypred")
run_wrapper_test("crps", crps, "ypred")
run_wrapper_test("srps", srps, "ypred")
run_wrapper_test("scrps", scrps, "ypred")
