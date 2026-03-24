test_that("insample_pred_measure returns class and metadata", {
  temp <- roaches_models$model1
  result <- insample_pred_measure(
    y = temp$y,
    mupred = temp$mupred,
    ylp = temp$ylp,
    measure = "r2"
  )

  expect_s3_class(result, "insample_pred_measure")
  expect_s3_class(result, "pred_measure")
  expect_equal(result$metadata$source, "insample")
  expect_true(all(c("elpd", "ic", "r2") %in% rownames(result$estimates)))
  expect_false("p_eff" %in% rownames(result$estimates))
})

test_that("kfold_pred_measure stores fold ids and optional p_eff", {
  temp <- roaches_models$model1
  fold_id <- rep(1:5, length.out = length(temp$y))

  result_no_peff <- kfold_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    measure = "elpd",
    fold_id = fold_id
  )
  expect_equal(result_no_peff$metadata$fold_id, fold_id)
  expect_false("p_eff" %in% rownames(result_no_peff$estimates))

  result_peff <- kfold_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    ylp_insample = temp$ylp,
    measure = "elpd",
    fold_id = fold_id
  )
  expect_true("p_eff" %in% rownames(result_peff$estimates))
})

test_that("test_pred_measure supports accumulation via predperf", {
  temp <- roaches_models$model1
  first <- test_pred_measure(
    y = temp$y,
    mupred = temp$mupred,
    ylp = temp$ylp,
    ylp_insample = temp$ylp,
    measure = "r2"
  )
  second <- test_pred_measure(
    y = temp$y,
    ypred = temp$ypred,
    measure = "rps",
    predperf = first
  )

  expect_s3_class(second, "test_pred_measure")
  expect_true(all(c("r2", "rps") %in% rownames(second$estimates)))
  expect_true(all(c("squared_error", "rps") %in% colnames(second$pointwise)))
})

test_that("kfold_pred_measure requires fold_id", {
  temp <- roaches_models$model1
  expect_error(
    kfold_pred_measure(
      y = temp$y,
      ylp = temp$ylp,
      measure = "elpd"
    ),
    "`fold_id` is required"
  )
})

test_that("test_pred_measure includes p_eff only when ylp_insample is provided", {
  temp <- roaches_models$model1

  without_peff <- test_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    measure = "elpd"
  )
  expect_false("p_eff" %in% rownames(without_peff$estimates))

  with_peff <- test_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    ylp_insample = temp$ylp,
    measure = "elpd"
  )
  expect_true("p_eff" %in% rownames(with_peff$estimates))
})

test_that("insample_pred_measure supports accumulation via predperf", {
  temp <- roaches_models$model1
  first <- insample_pred_measure(
    y = temp$y,
    mupred = temp$mupred,
    ylp = temp$ylp,
    measure = "mse"
  )
  second <- insample_pred_measure(
    y = temp$y,
    ypred = temp$ypred,
    measure = "rps",
    predperf = first
  )

  expect_s3_class(second, "insample_pred_measure")
  expect_true(all(c("mse", "rps") %in% rownames(second$estimates)))
  expect_true(all(c("squared_error", "rps") %in% colnames(second$pointwise)))
})

test_that("non-LOO wrappers accept multiple measures in one call", {
  temp <- roaches_models$model1

  multi <- insample_pred_measure(
    y = temp$y,
    mupred = temp$mupred,
    ylp = temp$ylp,
    measure = c("mse", "r2", "elpd", "logscore")
  )

  sequential <- insample_pred_measure(
    y = temp$y,
    mupred = temp$mupred,
    ylp = temp$ylp,
    measure = "mse"
  )
  sequential <- insample_pred_measure(
    y = temp$y,
    mupred = temp$mupred,
    measure = "r2",
    predperf = sequential
  )
  sequential <- insample_pred_measure(
    y = temp$y,
    mupred = temp$mupred,
    measure = "r2",
    predperf = sequential
  )
  sequential <- insample_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    measure = "elpd",
    predperf = sequential
  )

  expect_equal(
    multi$estimates[rownames(sequential$estimates), , drop = FALSE],
    sequential$estimates
  )
  expect_equal(
    multi$pointwise[, colnames(sequential$pointwise), drop = FALSE],
    sequential$pointwise
  )
  expect_true("logscore" %in% rownames(multi$estimates))
  expect_true("logscore" %in% colnames(multi$pointwise))
  expect_equal(multi$estimates["elpd", "Estimate"], sequential$estimates["elpd", "Estimate"])
})

test_that("pred_measure print includes source labels and kfold fold count", {
  temp <- roaches_models$model1
  fold_id <- rep(1:4, length.out = length(temp$y))

  ins <- insample_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    measure = "elpd"
  )
  kf <- kfold_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    measure = "elpd",
    fold_id = fold_id
  )
  tst <- test_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    measure = "elpd"
  )

  ins_out <- paste(capture.output(print(ins)), collapse = "\n")
  kf_out <- paste(capture.output(print(kf)), collapse = "\n")
  tst_out <- paste(capture.output(print(tst)), collapse = "\n")

  expect_match(ins_out, "Data source: in-sample", fixed = TRUE)
  expect_match(kf_out, "Data source: k-fold", fixed = TRUE)
  expect_match(kf_out, "Folds: 4", fixed = TRUE)
  expect_match(tst_out, "Data source: test", fixed = TRUE)
})

test_that("non-LOO functions reject incompatible predperf objects", {
  temp <- roaches_models$model1
  fold_id <- rep(1:4, length.out = length(temp$y))

  ins <- insample_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    measure = "elpd"
  )

  expect_error(
    test_pred_measure(
      y = temp$y,
      ylp = temp$ylp,
      measure = "elpd",
      predperf = ins
    ),
    "test_pred_measure"
  )
  expect_error(
    insample_pred_measure(
      y = temp$y,
      ylp = temp$ylp,
      measure = "elpd",
      predperf = temp$loo
    ),
    "Must inherit from class 'insample_pred_measure'"
  )
  expect_error(
    kfold_pred_measure(
      y = temp$y,
      ylp = temp$ylp,
      measure = "elpd",
      fold_id = fold_id,
      predperf = ins
    ),
    "kfold_pred_measure"
  )
})

test_that("non-LOO paths accept vector mupred and ylp", {
  temp <- roaches_models$model1
  fold_id <- rep(1:4, length.out = length(temp$y))

  mupred_vec <- temp$mupred[1, ]
  mupred_mat <- matrix(mupred_vec, nrow = 1)
  ylp_vec <- temp$ylp[1, ]
  ylp_mat <- matrix(ylp_vec, nrow = 1)

  ins_vec_metric <- insample_pred_measure(
    y = temp$y,
    mupred = mupred_vec,
    measure = "mse"
  )
  ins_mat_metric <- insample_pred_measure(
    y = temp$y,
    mupred = mupred_mat,
    measure = "mse"
  )
  expect_equal(ins_vec_metric$estimates, ins_mat_metric$estimates)
  expect_equal(ins_vec_metric$pointwise, ins_mat_metric$pointwise)

  test_vec_metric <- test_pred_measure(
    y = temp$y,
    mupred = mupred_vec,
    measure = "rmse"
  )
  test_mat_metric <- test_pred_measure(
    y = temp$y,
    mupred = mupred_mat,
    measure = "rmse"
  )
  expect_equal(test_vec_metric$estimates, test_mat_metric$estimates)
  expect_equal(test_vec_metric$pointwise, test_mat_metric$pointwise)

  kfold_vec_metric <- kfold_pred_measure(
    y = temp$y,
    mupred = mupred_vec,
    measure = "mae",
    fold_id = fold_id
  )
  kfold_mat_metric <- kfold_pred_measure(
    y = temp$y,
    mupred = mupred_mat,
    measure = "mae",
    fold_id = fold_id
  )
  expect_equal(kfold_vec_metric$estimates, kfold_mat_metric$estimates)
  expect_equal(kfold_vec_metric$pointwise, kfold_mat_metric$pointwise)

  ins_vec_score <- insample_pred_measure(
    ylp = ylp_vec,
    measure = "elpd"
  )
  ins_mat_score <- insample_pred_measure(
    ylp = ylp_mat,
    measure = "elpd"
  )
  expect_equal(ins_vec_score$estimates, ins_mat_score$estimates)
  expect_equal(ins_vec_score$pointwise, ins_mat_score$pointwise)

  test_vec_score <- test_pred_measure(
    ylp = ylp_vec,
    measure = "logscore"
  )
  test_mat_score <- test_pred_measure(
    ylp = ylp_mat,
    measure = "logscore"
  )
  expect_equal(test_vec_score$estimates, test_mat_score$estimates)
  expect_equal(test_vec_score$pointwise, test_mat_score$pointwise)

  test_vec_mlpd <- test_pred_measure(
    ylp = ylp_vec,
    measure = "mlpd"
  )
  test_mat_mlpd <- test_pred_measure(
    ylp = ylp_mat,
    measure = "mlpd"
  )
  expect_equal(test_vec_mlpd$estimates, test_mat_mlpd$estimates)
  expect_equal(test_vec_mlpd$pointwise, test_mat_mlpd$pointwise)

  kfold_vec_score <- kfold_pred_measure(
    ylp = ylp_vec,
    measure = "elpd",
    fold_id = fold_id
  )
  kfold_mat_score <- kfold_pred_measure(
    ylp = ylp_mat,
    measure = "elpd",
    fold_id = fold_id
  )
  expect_equal(kfold_vec_score$estimates, kfold_mat_score$estimates)
  expect_equal(kfold_vec_score$pointwise, kfold_mat_score$pointwise)
})

test_that("model_name is stored and preserved across non-LOO accumulation", {
  temp <- roaches_models$model1

  first <- insample_pred_measure(
    y = temp$y,
    ylp = temp$ylp,
    measure = "elpd",
    model_name = "m1"
  )
  expect_equal(first$metadata$model_name, "m1")

  second <- insample_pred_measure(
    y = temp$y,
    ypred = temp$ypred,
    measure = "rps",
    predperf = first
  )
  expect_equal(second$metadata$model_name, "m1")

  third <- insample_pred_measure(
    y = temp$y,
    ypred = temp$ypred,
    measure = "crps",
    predperf = second,
    model_name = "m2"
  )
  expect_equal(third$metadata$model_name, "m2")
})


test_that(".assert_predperf_compatibility returns invisibly for NULL predperf", {
  expect_invisible(.assert_predperf_compatibility("insample", NULL))
  expect_null(.assert_predperf_compatibility("insample", NULL))
})

test_that(".assert_predperf_compatibility accepts supported classes per source", {
  expect_invisible(
    .assert_predperf_compatibility(
      source = "loo",
      predperf = structure(list(), class = "loo")
    )
  )
  expect_invisible(
    .assert_predperf_compatibility(
      source = "loo",
      predperf = structure(list(), class = c("loo_pred_measure", "pred_measure"))
    )
  )
  expect_invisible(
    .assert_predperf_compatibility(
      source = "insample",
      predperf = structure(list(), class = c("insample_pred_measure", "pred_measure"))
    )
  )
  expect_invisible(
    .assert_predperf_compatibility(
      source = "kfold",
      predperf = structure(list(), class = c("kfold_pred_measure", "pred_measure"))
    )
  )
  expect_invisible(
    .assert_predperf_compatibility(
      source = "test",
      predperf = structure(list(), class = c("test_pred_measure", "pred_measure"))
    )
  )
})

test_that(".assert_predperf_compatibility rejects incompatible predperf classes", {
  expect_error(
    .assert_predperf_compatibility(
      source = "loo",
      predperf = structure(list(), class = c("insample_pred_measure", "pred_measure"))
    ),
    regexp = "inherit|class|loo"
  )
  expect_error(
    .assert_predperf_compatibility(
      source = "insample",
      predperf = structure(list(), class = "test_pred_measure")
    ),
    regexp = "inherit|class|insample_pred_measure"
  )
  expect_error(
    .assert_predperf_compatibility(
      source = "kfold",
      predperf = structure(list(), class = "insample_pred_measure")
    ),
    regexp = "inherit|class|kfold_pred_measure"
  )
  expect_error(
    .assert_predperf_compatibility(
      source = "test",
      predperf = structure(list(), class = "kfold_pred_measure")
    ),
    regexp = "inherit|class|test_pred_measure"
  )
})

test_that(".assert_finite_log_weight_columns rejects non-normalizable columns", {
  bad_log_weights <- matrix(
    c(-Inf, -Inf, -2, -3),
    nrow = 2
  )
  expect_error(
    .assert_finite_log_weight_columns(bad_log_weights),
    "has no finite values"
  )
})

test_that(".validate_measure_arg_not_null validates required args", {
  ylp <- NULL
  expect_invisible(.validate_measure_arg_not_null("elpd", 1))
  expect_error(
    .validate_measure_arg_not_null("elpd", ylp),
    "For elpd, the `ylp` argument must be specified but is currently NULL\\."
  )
  expect_error(
    .validate_measure_arg_not_null("elpd", NULL),
    "For elpd, the `NULL` argument must be specified but is currently NULL\\."
  )
})

test_that(".validate_measure_arg_not_null validates multiple args", {
  expect_invisible(
    .validate_measure_arg_not_null(
      "elpd",
      list(ylp = 1, log_weights = "ok")
    )
  )
  expect_error(
    .validate_measure_arg_not_null(
      "elpd",
      list(ylp = NULL, log_weights = 1)
    ),
    "For elpd, the `ylp` argument must be specified but is currently NULL\\."
  )
  expect_error(
    .validate_measure_arg_not_null(
      "elpd",
      list(ylp = NULL, log_weights = NULL)
    ),
    "For elpd, the following arguments must be specified but are currently NULL: `ylp`, `log_weights`\\."
  )
})

test_that(".validate_measure_names supports vectors and deduplicates aliases", {
  choices <- c("logscore", "mlpd", "elpd", "r2", "mse")

  # Invalid values throw an error
  expect_error(
    .validate_measure_names(c("logscor", "lpd", "mse", "mse"), choices),
    regexp = "Invalid measure\\(s\\): 'logscor', 'lpd'"
  )

  expect_error(
    .validate_measure_names(c("logscore", "lpd", "mse", "mse"), choices),
    regexp = "Invalid measure\\(s\\): 'lpd'"
  )

  expect_no_error(
    .validate_measure_names(c("logscore", "elpd", "mse", "mse"), choices)
  )

  expect_no_error(
    .validate_measure_names("logscore", choices)
  )

  expect_error(
    .validate_measure_names("logscor", choices)
  )
})

