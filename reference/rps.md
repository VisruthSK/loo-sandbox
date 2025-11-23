# (Continuous) Ranked Probability Score

Given draws ypred from the predictive distribution and the observations
y, computes

- rank probability score (Epstein, 1969) for discrete ypred and y and
  continuous rank probability score (Matheson and Winkler, 1976;
  Gneiting & Raftery, 2007) for continuous ypred and y, using the
  probability weighted moment form (Taillardat et al., 2016; Zamo &
  Naveau, 2017)

- scaled versions of these, if `scaled=TRUE` (Bolin & Wallin, 2023).

## Usage

``` r
rps(y, ypred, log_weights)

crps(y, ypred, log_weights)

srps(y, ypred, log_weights)

scrps(y, ypred, log_weights)
```

## Arguments

- y:

  vector of observed values (n)

- ypred:

  matrix of posterior draws (S x n) of posterior predictive draws

- log_weights:

  matrix of standardized loo weights (S x n) on the log scale

## Details

Utility version of the score is returned, that is, bigger is better, to
match the utility version of log score / elpd (the original rank
probability score by Epstein (1969) was also in this direction).

The same sample based \$L\$-moment estimator is used for continuous and
discrete variables. It is commonly stated that the probability weighted
moment form assumes F(x) is continuous. However, Hosking (1990) states
that \$L\$-moments can be used with discrete distributions
`` provided that the quantile function is `normalized' in the sense of Widder (1941).'' Hosking (1996) states the same condition more simply as  ``A
discrete random variable can be approximated arbitrarily closely by a
continuous random variable, so the result is also valid for discrete
random variables”.

## References

- Bolin, D. and Wallin, J. (2023). Local scale invariance and robustness
  of proper scoring rules. *Statistical Science*, 38(1):140-159.

- Epstein, E.S. (1969). A scoring system for probability forecasts of
  ranked categories. *Journal of Applied Meteorology*, 8(6):985-987.

- Gneiting, T. and Raftery, A.E. (2007). Strictly Proper Scoring Rules,
  Prediction, and Estimation. *Journal of the American Statistical
  Association*, 102(477):359-378.

- Hosking, J.R.M. (1990). \$L\$-moments: analysis and estimation of
  distributions using linear combinations of order statistics. *Journal
  of the Royal Statistical Society Series B: Statistical Methodology*,
  52(1):105-124.

- Hosking, J.R.M. (1996). Some theoretical results concerning
  \$L\$-moments. Research report RC 14492. IBM Thomas J. Watson Research
  Division.

- Matheson, J.E., and Winkler, R.L. (1976). Scoring Rules for Continuous
  Probability Distributions. *Management Science*, 22(10), 1087-1096.

- Taillardat, M., Mestre, O., Zamo, M., and Naveau, P. (2016).
  Calibrated Ensemble Forecasts Using Quantile Regression Forests and
  Ensemble Model Output Statistics. *Monthly Weather Review*, 144(6),
  2375-2393.

- Widder, D.V. (1941). *The Laplace Transform*. Princeton: Princeton
  University Press.

- Zamo, M., and Naveau, P. (2018). Estimation of the Continuous Ranked
  Probability Score with Limited Information and Applications to
  Ensemble Weather Forecasts. *Mathematical Geosciences*, 50, 209–234.
