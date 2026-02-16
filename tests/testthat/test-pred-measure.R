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
