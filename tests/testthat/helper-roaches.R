# Shared roaches model inputs used across multiple measure tests
roaches_models <- local({
  rename_loo_columns <- function(loo_obj) {
    old_to_new <- c(
      elpd_loo = "elpd",
      p_loo = "p_eff",
      looic = "ic",
      mcse_elpd_loo = "mcse_elpd",
      n_eff = "n_eff",
      pareto_k = "pareto_k"
    )
    rename_vec <- function(x) {
      if (is.null(x)) {
        return(x)
      }
      keep <- x %in% names(old_to_new)
      x[keep] <- unname(old_to_new[x[keep]])
      x
    }

    rownames(loo_obj$estimates) <- rename_vec(rownames(loo_obj$estimates))
    colnames(loo_obj$pointwise) <- rename_vec(colnames(loo_obj$pointwise))
    loo_obj
  }

  load_model <- function(id) {
    fit <- readRDS(testthat::test_path("fits", sprintf("fit%d.RDS", id)))
    loo_obj <- readRDS(testthat::test_path("fits", sprintf("loo%d.RDS", id)))
    loo_obj <- rename_loo_columns(loo_obj)
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
