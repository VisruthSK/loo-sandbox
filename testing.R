library(loo)
options(mc.cores = 4)

loo1 <- readRDS("loo1.RDS")
loo2 <- readRDS("loo2.RDS")

loo_compare(loo1, loo2)

str(loo1)

# keep only first 3 elements
new_loo1 <- loo1[1:3]

str(new_loo1)

add_measure <- function(my_loo, y, ypred, measure=c("RPS", "SRPS", "CRPS", "SCRPS", "mae", "rmse", "mse", "R2", "acc", "balanced_acc"){
  measure <- match.arg(measure)
  # dispatch on measure to calculate new results
  # measure funcs should return vector: estimate, SE
}