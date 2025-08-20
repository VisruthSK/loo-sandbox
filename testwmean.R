set.seed(seed)
y <- rnorm(100) |> sort()
yhat <- rnorm(100, 0.1) |> sort()
weights <- rexp(100) |> sort() |> prop.table()

tmp <- function(weights) {
  print(any(missing(weights)))
  list(mean((y - yhat)^2), weighted.mean((y - yhat)^2, weights))
}

tmp(NULL)
