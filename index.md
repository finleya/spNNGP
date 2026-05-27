# spNNGP

[![CRAN status](https://www.r-pkg.org/badges/version/spNNGP)](https://CRAN.R-project.org/package=spNNGP)
[![R-CMD-check](https://github.com/finleya/spNNGP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/finleya/spNNGP/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/finleya/spNNGP/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/finleya/spNNGP/actions/workflows/pkgdown.yaml)

`spNNGP` fits univariate Bayesian spatial regression models for large
datasets using nearest neighbor Gaussian processes (NNGPs). NNGP models
replace dense Gaussian process calculations with sparse neighbor-based
approximations, making fully Bayesian spatial regression practical for
datasets that are too large for conventional Gaussian process methods.

## Installation

Install the CRAN release with:

```r
install.packages("spNNGP")
```

Install the development version from GitHub with:

```r
remotes::install_github("finleya/spNNGP")
```

## Basic Use

```r
library(spNNGP)

set.seed(1)
n <- 100
coords <- cbind(runif(n), runif(n))
x <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(n)

fit <- spNNGP(
  y ~ x,
  coords = coords,
  method = "latent",
  n.neighbors = 10,
  starting = list(beta = c(0, 0), sigma.sq = 1, tau.sq = 1, phi = 3),
  tuning = list(phi = 0.1),
  priors = list(
    sigma.sq.IG = c(2, 1),
    tau.sq.IG = c(2, 1),
    phi.Unif = c(0.1, 30)
  ),
  n.samples = 100,
  n.omp.threads = 1,
  verbose = FALSE
)

summary(fit)
```

## Functionality

The package provides:

- Response and latent NNGP models for Gaussian outcomes.
- Latent NNGP models for binomial outcomes.
- Conjugate NNGP models for faster Gaussian analyses.
- Prediction and posterior summaries through familiar R methods.
- MCMC diagnostics and model summaries for fitted spatial models.

## Citations

If you use `spNNGP`, please cite:

Finley, A. O., Datta, A., and Banerjee, S. (2022). spNNGP R Package for
Nearest Neighbor Gaussian Process Models. *Journal of Statistical
Software*, 103(5). doi:10.18637/jss.v103.i05.

Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E.,
and Banerjee, S. (2019). Efficient algorithms for Bayesian nearest
neighbor Gaussian processes. *Journal of Computational and Graphical
Statistics*, 28(2), 401-414. doi:10.1080/10618600.2018.1537924.
Preprint: <https://arxiv.org/pdf/1702.00434>
