test_that("conjugate NNGP returns summaries, fitted values, and predictions", {
  set.seed(12)
  dat <- make_gaussian_data(n = 35)

  fit <- spConjNNGP(
    dat$y ~ dat$x,
    coords = dat$coords,
    n.neighbors = 6,
    theta.alpha = c(phi = 5, alpha = 0.5),
    sigma.sq.IG = c(2, 1),
    cov.model = "exponential",
    fit.rep = TRUE,
    n.samples = 5,
    n.omp.threads = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "NNGP")
  expect_equal(fit$type, c("conjugate", "gaussian"))
  expect_equal(dim(fit$p.beta.theta.samples), c(5L, 4L))
  expect_finite_matrix(fit$y.hat.samples, nrow = 35L, ncol = 5L)
  expect_finite_matrix(fit$y.rep.samples, nrow = 35L, ncol = 5L)
  expect_true(all(c("n.neighbors", "nn.indx", "nn.indx.lu", "ord") %in% names(fit$neighbor.info)))

  fit_out <- fitted(fit)
  expect_equal(dim(fit_out$y.hat.samples), c(35L, 5L))

  pred <- predict(
    fit,
    X.0 = cbind(1, dat$x[1:3]),
    coords.0 = dat$coords[1:3, ],
    n.omp.threads = 1,
    verbose = FALSE
  )

  expect_s3_class(pred, "predict.NNGP")
  expect_finite_matrix(pred$y.0.hat, nrow = 3L, ncol = 1L)
  expect_finite_matrix(pred$y.0.hat.var, nrow = 3L, ncol = 1L)
})

test_that("conjugate NNGP validates theta alpha names", {
  set.seed(13)
  dat <- make_gaussian_data(n = 20)

  expect_error(
    spConjNNGP(
      dat$y ~ dat$x,
      coords = dat$coords,
      n.neighbors = 5,
      theta.alpha = c(phi = 5),
      sigma.sq.IG = c(2, 1),
      cov.model = "exponential",
      n.omp.threads = 1,
      verbose = FALSE
    ),
    "column name alpha must be specified"
  )
})
