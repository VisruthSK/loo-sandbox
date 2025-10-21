library(loo)
source("loo_measure.R")

loo1 <- readRDS("loo1.RDS")
fit1 <- readRDS("fit1.RDS")

loo1$estimates
loo1$pointwise[, "elpd_loo"]

ylp <- rstanarm::log_lik(fit1)
loo_weights <- weights(loo1$psis_object, log = FALSE, normalize = TRUE)
.elpd_summary(ylp, loo_weights)

loo2 <- readRDS("loo2.RDS")
fit2 <- readRDS("fit2.RDS")

loo2$estimates
loo2$pointwise[, "elpd_loo"]

ylp <- rstanarm::log_lik(fit2)
loo_weights <- weights(loo2$psis_object, log = FALSE, normalize = TRUE)
.elpd_summary(ylp, loo_weights)
