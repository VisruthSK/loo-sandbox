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
