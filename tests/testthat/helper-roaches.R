# Shared roaches model inputs used across multiple measure tests
roaches_models <- local({
  load_model <- function(id) {
    fit <- readRDS(testthat::test_path("fits", sprintf("fit%d.RDS", id)))
    loo_obj <- readRDS(testthat::test_path("fits", sprintf("loo%d.RDS", id)))
    y <- as.numeric(rstanarm::get_y(fit))
    mupred <- rstanarm::posterior_epred(fit)
    ypred <- with(set.seed(0), rstanarm::posterior_predict(fit))
    S <- nrow(mupred)
    list(
      y = y,
      y_binary = as.integer(y > 0),
      ypred = ypred,
      mupred = mupred,
      binary_mupred = ifelse(ypred > 0, 1L, 0L),
      loo = loo_obj,
      ylp = rstanarm::log_lik(fit)
    )
  }

  list(
    model1 = load_model(1),
    model2 = load_model(2)
  )
})
