test_that("PGLogit supports fitted, residual, and summary methods", {
  set.seed(14)
  n <- 35
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(-0.2 + 0.8 * x))

  fit <- PGLogit(
    y ~ x,
    n.samples = 8,
    n.omp.threads = 1,
    fit.rep = TRUE,
    sub.sample = list(start = 4, end = 8, thin = 2),
    verbose = FALSE
  )

  expect_s3_class(fit, "PGLogit")
  expect_equal(dim(fit$p.beta.samples), c(8L, 2L))
  expect_equal(fit$s.indx, c(4L, 6L, 8L))
  expect_finite_matrix(fit$y.hat.samples, nrow = 35L, ncol = 3L)
  expect_finite_matrix(fit$y.hat.quants, nrow = 35L, ncol = 3L)
  expect_finite_matrix(fit$y.rep.samples, nrow = 35L, ncol = 3L)

  fit_out <- fitted(fit)
  expect_equal(dim(fit_out$y.hat.samples), c(35L, 3L))
  expect_equal(dim(fit_out$y.rep.samples), c(35L, 3L))

  res_out <- residuals(fit)
  expect_equal(dim(res_out$residuals.samples), c(35L, 3L))

  capture.output(
    sum_out <- summary(fit, sub.sample = list(start = 4, end = 8, thin = 2))
  )
  expect_equal(dim(sum_out), c(2L, 5L))
})

test_that("PGLogit validates binomial weights", {
  set.seed(15)
  y <- rbinom(10, 1, 0.5)
  x <- rnorm(10)

  expect_error(
    PGLogit(y ~ x, weights = 1:3, n.samples = 5, n.omp.threads = 1, verbose = FALSE),
    "weights must be of length n or 1"
  )
})
