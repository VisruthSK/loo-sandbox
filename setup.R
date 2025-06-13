library(loo)
library(rstanarm)
options(mc.cores = 4)
data(roaches)

fit1 <- stan_glm(
  formula = y ~ roach1 + treatment + senior,
  offset = log(exposure2),
  data = roaches,
  family = poisson(link = "log"),
  prior = normal(0, 2.5, autoscale = TRUE),
  prior_intercept = normal(0, 5, autoscale = TRUE),
  seed = 12345
)

loo1 <- loo(fit1, save_psis = TRUE)

fit2 <- update(fit1, family = neg_binomial_2)
loo2 <- loo(fit2, save_psis = TRUE)

loo_compare(loo1, loo2)

saveRDS(loo1, "loo1.RDS")
saveRDS(loo2, "loo2.RDS")
