# library(loo)
# library(rstanarm)
# options(mc.cores = 4)
# data(roaches)

# fit1 <- stan_glm(
#   formula = y ~ roach1 + treatment + senior,
#   offset = log(exposure2),
#   data = roaches,
#   family = poisson(link = "log"),
#   prior = normal(0, 2.5, autoscale = TRUE),
#   prior_intercept = normal(0, 5, autoscale = TRUE),
#   seed = 12345
# )

# loo1 <- loo(fit1, save_psis = TRUE)

# fit2 <- update(fit1, family = neg_binomial_2)
# loo2 <- loo(fit2, save_psis = TRUE)

# loo_compare(loo1, loo2)

# # code from https://mc-stan.org/loo/articles/loo2-example.html

# saveRDS(loo1, "loo1.RDS")
# saveRDS(loo2, "loo2.RDS")

# saveRDS(fit1, "fit1.RDS")
# saveRDS(fit2, "fit2.RDS")

test_that("first test using fit from roaches", {
  loo <- readRDS(test_path("fits", "loo1.RDS"))
  fit <- readRDS(test_path("fits", "fit1.RDS"))

  ylp <- rstanarm::log_lik(fit)
  loo_weights <- weights(loo$psis_object, normalize = TRUE)
  new <- .elpd_summary(ylp, loo_weights)
  expect_equal(new$pointwise, loo$pointwise[, "elpd_loo"])
  expect_equal(
    c(Estimate = new$estimate, SE = new$se),
    loo$estimate["elpd_loo", ]
  )
})

test_that("second test using fit from roaches", {
  loo <- readRDS(test_path("fits", "loo2.RDS"))
  fit <- readRDS(test_path("fits", "fit2.RDS"))

  ylp <- rstanarm::log_lik(fit)
  loo_weights <- weights(loo$psis_object, normalize = TRUE)
  new <- .elpd_summary(ylp, loo_weights)
  expect_equal(new$pointwise, loo$pointwise[, "elpd_loo"])
  expect_equal(
    c(Estimate = new$estimate, SE = new$se),
    loo$estimate["elpd_loo", ]
  )
})
