make_gaussian_data <- function(n = 40, beta = c(1, 0.5), sd = 1) {
  coords <- cbind(runif(n), runif(n))
  x <- rnorm(n)
  y <- beta[1] + beta[2] * x + rnorm(n, sd = sd)

  list(
    y = y,
    x = x,
    coords = coords,
    starting = list(beta = c(0, 0), sigma.sq = 1, tau.sq = 1, phi = 3),
    tuning = list(phi = 0.1),
    response.tuning = list(sigma.sq = 0.1, tau.sq = 0.1, phi = 0.1),
    priors = list(
      sigma.sq.IG = c(2, 1),
      tau.sq.IG = c(2, 1),
      phi.Unif = c(0.1, 30)
    )
  )
}

make_binomial_data <- function(n = 40) {
  coords <- cbind(runif(n), runif(n))
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(-0.2 + 0.8 * x))

  list(
    y = y,
    x = x,
    coords = coords,
    starting = list(beta = c(0, 0), sigma.sq = 1, phi = 3),
    tuning = list(phi = 0.1),
    priors = list(
      sigma.sq.IG = c(2, 1),
      phi.Unif = c(0.1, 30)
    )
  )
}

expect_finite_matrix <- function(x, nrow = NULL, ncol = NULL) {
  expect_true(is.matrix(x))
  expect_true(all(is.finite(x)))
  if (!is.null(nrow)) {
    expect_equal(nrow(x), nrow)
  }
  if (!is.null(ncol)) {
    expect_equal(ncol(x), ncol)
  }
}
