roaches_manual_inputs <- local({
  load_model <- function(id) {
    fit <- readRDS(test_path("fits", sprintf("fit%d.RDS", id)))
    y <- as.numeric(rstanarm::get_y(fit))
    draws_pred <- rstanarm::posterior_predict(fit)
    list(
      y = y,
      y_binary = as.integer(y > 0),
      ypred = draws_pred,
      mupred = rstanarm::posterior_epred(fit),
      binary_mupred = ifelse(draws_pred > 0, 1L, 0L),
      ylp = rstanarm::log_lik(fit)
    )
  }

  list(
    model1 = load_model(1),
    model2 = load_model(2)
  )
})

manual_equal_weights_se <- function(values) {
  n <- length(values)
  dm <- posterior::as_draws_matrix(matrix(values, ncol = 1))
  sd_val <- posterior::summarise_draws(dm, "sd")$sd
  sd_val / sqrt(posterior::ndraws(dm))
}

manual_col_means <- function(draws) {
  dm <- posterior::as_draws_matrix(draws)
  posterior::summarise_draws(dm, "mean")$mean
}

manual_accuracy_pointwise <- function(y, mupred) {
  eq_mat <- sweep(mupred, 2, y, `==`)
  manual_col_means(eq_mat)
}

manual_balanced_accuracy_pointwise <- function(y, mupred) {
  acc <- manual_accuracy_pointwise(y, mupred)
  cls_counts <- table(y)
  cls_lookup <- as.numeric(cls_counts[match(y, names(cls_counts))])
  acc * length(y) / (length(cls_counts) * cls_lookup)
}

manual_mae_summary <- function(y, mupred) {
  pointwise <- abs(y - manual_col_means(mupred))
  estimate <- mean(pointwise)
  se <- manual_equal_weights_se(pointwise)
  list(pointwise = pointwise, estimate = estimate, se = se)
}

manual_mse_summary <- function(y, mupred) {
  pointwise <- (y - manual_col_means(mupred))^2
  estimate <- mean(pointwise)
  se <- manual_equal_weights_se(pointwise)
  list(pointwise = pointwise, estimate = estimate, se = se)
}

manual_rmse_summary <- function(y, mupred) {
  mse_pointwise <- (y - manual_col_means(mupred))^2
  mean_mse <- mean(mse_pointwise)
  estimate <- sqrt(mean_mse)
  se <- manual_equal_weights_se(mse_pointwise) / (2 * sqrt(mean_mse))
  list(pointwise = mse_pointwise, estimate = estimate, se = se)
}

manual_r2_summary <- function(y, mupred) {
  mse_pointwise <- (y - manual_col_means(mupred))^2
  n <- length(y)
  mse_loo <- mean(mse_pointwise)
  se_mse_loo <- manual_equal_weights_se(mse_pointwise)
  centered_y <- y - mean(y)
  squared_error_y <- centered_y^2
  mse_y <- mean(squared_error_y)
  cov_term <- stats::cov(mse_pointwise, squared_error_y)
  var_term <- stats::var(squared_error_y)
  se_r2 <- sqrt(
    se_mse_loo^2 -
      2 * (mse_loo / mse_y) * cov_term / n +
      (mse_loo^2 / mse_y^2) * var_term / n
  ) /
    mse_y
  list(
    pointwise = mse_pointwise,
    estimate = 1 - mse_loo / mse_y,
    se = se_r2
  )
}

manual_elpd_summary <- function(ylp) {
  S <- nrow(ylp)
  pointwise <- matrixStats::colLogSumExps(ylp) - log(S)
  n <- length(pointwise)
  estimate <- sum(pointwise)
  se <- n * manual_equal_weights_se(pointwise)
  list(pointwise = pointwise, estimate = estimate, se = se)
}

manual_logscore_summary <- function(ylp) {
  elpd <- manual_elpd_summary(ylp)
  n <- length(elpd$pointwise)
  list(
    pointwise = elpd$pointwise,
    estimate = elpd$estimate / n,
    se = elpd$se / n
  )
}

manual_rps_point <- function(obs, draws, scaled = FALSE) {
  S <- length(draws)
  log_weight <- -log(S)
  exy <- log_weight * sum(abs(obs - draws))
  exx <- -2 * log_weight * obs
  if (!scaled) {
    -exy + 0.5 * exx
  } else {
    if (exx <= 0) {
      return(NaN)
    }
    -exy / exx - 0.5 * log(exx)
  }
}

manual_rps_summary <- function(y, ypred, scaled = FALSE) {
  pointwise <- vapply(
    seq_along(y),
    function(i) {
      manual_rps_point(y[i], ypred[, i], scaled)
    },
    numeric(1)
  )
  estimate <- mean(pointwise)
  se <- manual_equal_weights_se(pointwise)
  list(pointwise = pointwise, estimate = estimate, se = se)
}

manual_specs <- list(
  balanced_acc = list(
    prepare = function(inputs) {
      pointwise <- manual_balanced_accuracy_pointwise(
        inputs$y_binary,
        inputs$binary_mupred
      )
      estimate <- mean(pointwise)
      se <- manual_equal_weights_se(pointwise)
      list(pointwise = pointwise, estimate = estimate, se = se)
    },
    call_args = function(inputs) {
      list(
        y = inputs$y_binary,
        mupred = inputs$binary_mupred,
        measure = "balanced_acc"
      )
    }
  ),
  acc = list(
    prepare = function(inputs) {
      pointwise <- manual_accuracy_pointwise(
        inputs$y_binary,
        inputs$binary_mupred
      )
      estimate <- mean(pointwise)
      se <- manual_equal_weights_se(pointwise)
      list(pointwise = pointwise, estimate = estimate, se = se)
    },
    call_args = function(inputs) {
      list(
        y = inputs$y_binary,
        mupred = inputs$binary_mupred,
        measure = "acc"
      )
    }
  ),
  mae = list(
    prepare = function(inputs) manual_mae_summary(inputs$y, inputs$mupred),
    call_args = function(inputs) {
      list(y = inputs$y, mupred = inputs$mupred, measure = "mae")
    }
  ),
  mse = list(
    prepare = function(inputs) manual_mse_summary(inputs$y, inputs$mupred),
    call_args = function(inputs) {
      list(y = inputs$y, mupred = inputs$mupred, measure = "mse")
    }
  ),
  rmse = list(
    prepare = function(inputs) manual_rmse_summary(inputs$y, inputs$mupred),
    call_args = function(inputs) {
      list(y = inputs$y, mupred = inputs$mupred, measure = "rmse")
    }
  ),
  r2 = list(
    prepare = function(inputs) manual_r2_summary(inputs$y, inputs$mupred),
    call_args = function(inputs) {
      list(y = inputs$y, mupred = inputs$mupred, measure = "r2")
    }
  ),
  elpd = list(
    prepare = function(inputs) manual_elpd_summary(inputs$ylp),
    call_args = function(inputs) {
      list(y = inputs$y, ylp = inputs$ylp, measure = "elpd")
    }
  ),
  logscore = list(
    prepare = function(inputs) manual_logscore_summary(inputs$ylp),
    call_args = function(inputs) {
      list(y = inputs$y, ylp = inputs$ylp, measure = "logscore")
    }
  ),
  rps = list(
    prepare = function(inputs) manual_rps_summary(inputs$y, inputs$ypred),
    call_args = function(inputs) {
      list(y = inputs$y, ypred = inputs$ypred, measure = "rps")
    }
  ),
  crps = list(
    prepare = function(inputs) manual_rps_summary(inputs$y, inputs$ypred),
    call_args = function(inputs) {
      list(y = inputs$y, ypred = inputs$ypred, measure = "crps")
    }
  ),
  srps = list(
    prepare = function(inputs) {
      manual_rps_summary(inputs$y, inputs$ypred, scaled = TRUE)
    },
    call_args = function(inputs) {
      list(y = inputs$y, ypred = inputs$ypred, measure = "srps")
    }
  ),
  scrps = list(
    prepare = function(inputs) {
      manual_rps_summary(inputs$y, inputs$ypred, scaled = TRUE)
    },
    call_args = function(inputs) {
      list(y = inputs$y, ypred = inputs$ypred, measure = "scrps")
    }
  )
)

for (measure in names(manual_specs)) {
  spec <- manual_specs[[measure]]
  test_that(sprintf("manual %s matches loo_pred_measure", measure), {
    for (model_name in names(roaches_manual_inputs)) {
      inputs <- roaches_manual_inputs[[model_name]]
      manual <- spec$prepare(inputs)
      args <- spec$call_args(inputs)
      result <- do.call(loo_pred_measure, args)
      expect_equal(
        unname(result$pointwise),
        unname(manual$pointwise),
        info = paste(measure, model_name),
        tolerance = 1e-8
      )
      expect_equal(
        result$estimate,
        manual$estimate,
        info = paste(measure, model_name),
        tolerance = 1e-8
      )
      expect_equal(
        result$se,
        manual$se,
        info = paste(measure, model_name),
        tolerance = 1e-8
      )
    }
  })
}
