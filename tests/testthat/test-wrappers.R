# Helper to run the comparison loop
run_wrapper_test <- function(
  measure,
  wrapper_fn
) {
  wrapper_name <- deparse(substitute(wrapper_fn))

  test_that(paste(measure, "wrapper matches loo_pred_measure"), {
    for (model_name in names(roaches_models)) {
      inputs <- roaches_models[[model_name]]
      log_weights <- weights(inputs$loo$psis_object)

      switch(
        wrapper_name,
        mse = ,
        rmse = ,
        mae = ,
        r2 = {
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
        },
        acc = ,
        balanced_acc = {
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
        },
        elpd = ,
        logscore = {
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
        },
        rps = ,
        crps = ,
        srps = ,
        scrps = {
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
        },
        stop("Unknown wrapper: ", wrapper_name)
      )

      expect_equal(result, expected, info = model_name)
    }
  })
}

run_wrapper_test("mse", mse)
run_wrapper_test("rmse", rmse)
run_wrapper_test("mae", mae)
run_wrapper_test("r2", r2)

run_wrapper_test("acc", acc)
run_wrapper_test("balanced_acc", balanced_acc)

run_wrapper_test("elpd", elpd)
run_wrapper_test("logscore", logscore)

run_wrapper_test("rps", rps)
run_wrapper_test("crps", crps)
run_wrapper_test("srps", srps)
run_wrapper_test("scrps", scrps)
