test_that("balanced accuracy matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.balanced_accuracy_summary(
      inputs$y_binary,
      inputs$binary_mupred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y_binary,
      mupred = inputs$binary_mupred,
      measure = "balanced_acc"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("accuracy matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.accuracy_summary(
      inputs$y_binary,
      inputs$binary_mupred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y_binary,
      mupred = inputs$binary_mupred,
      measure = "acc"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("mae matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.mae_summary(
      inputs$y,
      inputs$mupred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      mupred = inputs$mupred,
      measure = "mae"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("mse matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.mse_summary(
      inputs$y,
      inputs$mupred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      mupred = inputs$mupred,
      measure = "mse"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("rmse matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.rmse_summary(
      inputs$y,
      inputs$mupred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      mupred = inputs$mupred,
      measure = "rmse"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("r2 matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.r2_summary(
      inputs$y,
      inputs$mupred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      mupred = inputs$mupred,
      measure = "r2"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("logscore matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.logscore_summary(
      inputs$ylp,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      ylp = inputs$ylp,
      measure = "logscore"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("rps matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.rps_summary(
      inputs$y,
      inputs$ypred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      ypred = inputs$ypred,
      measure = "rps"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("crps matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.rps_summary(
      inputs$y,
      inputs$ypred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      ypred = inputs$ypred,
      measure = "crps"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("srps matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.srps_summary(
      inputs$y,
      inputs$ypred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      ypred = inputs$ypred,
      measure = "srps"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("scrps matches internal summary for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    manual <- loosandbox:::.srps_summary(
      inputs$y,
      inputs$ypred,
      inputs$log_weights
    )
    result <- loo_pred_measure(
      y = inputs$y,
      ypred = inputs$ypred,
      measure = "scrps"
    )
    expect_equal(
      unname(result$pointwise),
      unname(manual$pointwise),
      info = model_name
    )
    expect_equal(result$estimate, manual$estimate, info = model_name)
    expect_equal(result$se, manual$se, info = model_name)
  }
})

test_that("elpd matches loo package values for roaches models", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]
    loo_log_weights <- weights(
      inputs$loo$psis_object,
      normalize = TRUE
    )
    result <- loo_pred_measure(
      y = inputs$y,
      ylp = inputs$ylp,
      log_weights = loo_log_weights,
      measure = "elpd"
    )
    expect_equal(
      unname(result$pointwise),
      unname(inputs$loo$pointwise[, "elpd_loo"]),
      info = model_name
    )
    expect_equal(
      c(result$estimate, result$se),
      unname(inputs$loo$estimate["elpd_loo", c("Estimate", "SE")]),
      info = model_name
    )
  }
})

test_that("logscore with loo object reuses pointwise values and updates loo object", {
  for (model_name in names(roaches_models)) {
    inputs <- roaches_models[[model_name]]

    # Get PSIS weights from the loo object (same as used in elpd test)
    loo_log_weights <- weights(
      inputs$loo$psis_object,
      normalize = TRUE
    )

    # Create a copy of loo object to test that it gets updated
    loo_copy <- inputs$loo

    # Add an elpd_loo column to the loo object (initially with dummy values)
    # Note: logscore uses "elpd_loo" column per .match_pointwise_column()
    loo_copy$pointwise <- cbind(
      loo_copy$pointwise,
      elpd_loo = rep(NA_real_, nrow(loo_copy$pointwise))
    )

    # Call loo_pred_measure with the loo object and PSIS weights
    result <- loo_pred_measure(
      y = inputs$y,
      ylp = inputs$ylp,
      loo = loo_copy,
      log_weights = loo_log_weights,
      measure = "logscore"
    )

    # Expected pointwise values: should match the loo object's elpd_loo values
    # (since logscore and elpd share the same pointwise values)
    expected_pointwise <- inputs$loo$pointwise[, "elpd_loo"]

    expect_equal(
      unname(result$pointwise),
      unname(expected_pointwise),
      info = paste(model_name, "- pointwise values")
    )

    # Check that the loo object's pointwise column was updated
    expect_equal(
      unname(loo_copy$pointwise[, "elpd_loo"]),
      unname(expected_pointwise),
      info = paste(model_name, "- loo object updated")
    )

    # Check that estimate and SE are correct
    # logscore is elpd divided by n
    n <- length(expected_pointwise)
    expected_estimate <- sum(expected_pointwise) / n
    expected_se <- n *
      sqrt(
        sum((expected_pointwise - mean(expected_pointwise))^2) / (n * (n - 1))
      ) /
      n

    expect_equal(
      result$estimate,
      expected_estimate,
      info = paste(model_name, "- estimate")
    )

    expect_equal(
      result$se,
      expected_se,
      info = paste(model_name, "- se")
    )
  }
})
