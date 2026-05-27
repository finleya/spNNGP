test_that("response Gaussian sampler supports fitted, residual, and prediction methods", {
  set.seed(10)
  dat <- make_gaussian_data(n = 35)

  fit <- spNNGP(
    dat$y ~ dat$x,
    coords = dat$coords,
    method = "response",
    n.neighbors = 6,
    starting = dat$starting,
    tuning = dat$response.tuning,
    priors = dat$priors,
    n.samples = 8,
    fit.rep = TRUE,
    sub.sample = list(start = 4, end = 8, thin = 2),
    n.omp.threads = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "NNGP")
  expect_equal(fit$type, c("response", "gaussian"))
  expect_equal(dim(fit$p.beta.samples), c(8L, 2L))
  expect_equal(dim(fit$p.theta.samples), c(8L, 3L))
  expect_equal(fit$s.indx, c(4L, 6L, 8L))
  expect_finite_matrix(fit$y.hat.samples, nrow = 35L, ncol = 3L)
  expect_finite_matrix(fit$y.rep.samples, nrow = 35L, ncol = 3L)

  fit_out <- fitted(fit)
  expect_equal(dim(fit_out$y.hat.samples), c(35L, 3L))
  expect_equal(dim(fit_out$y.rep.samples), c(35L, 3L))

  res_out <- residuals(fit)
  expect_equal(dim(res_out$residuals.samples), c(35L, 3L))

  pred <- predict(
    fit,
    X.0 = cbind(1, dat$x[1:3]),
    coords.0 = dat$coords[1:3, ],
    sub.sample = list(start = 4, end = 8, thin = 2),
    n.omp.threads = 1,
    verbose = FALSE
  )

  expect_s3_class(pred, "predict.NNGP")
  expect_finite_matrix(pred$p.y.0, nrow = 3L, ncol = 3L)
})

test_that("neighbor information can be reused across response fits", {
  set.seed(11)
  dat <- make_gaussian_data(n = 30)

  fit <- spNNGP(
    dat$y ~ dat$x,
    coords = dat$coords,
    method = "response",
    n.neighbors = 5,
    starting = dat$starting,
    tuning = dat$response.tuning,
    priors = dat$priors,
    n.samples = 5,
    return.neighbor.info = TRUE,
    n.omp.threads = 1,
    verbose = FALSE
  )

  expect_true(all(c("n.neighbors", "nn.indx", "nn.indx.lu", "ord") %in% names(fit$neighbor.info)))

  refit <- suppressWarnings(spNNGP(
    dat$y ~ dat$x,
    coords = dat$coords,
    method = "response",
    starting = dat$starting,
    tuning = dat$response.tuning,
    priors = dat$priors,
    n.samples = 5,
    neighbor.info = fit$neighbor.info,
    n.omp.threads = 1,
    verbose = FALSE
  ))

  expect_s3_class(refit, "NNGP")
  expect_equal(refit$n.neighbors, fit$n.neighbors)
  expect_equal(dim(refit$p.theta.samples), c(5L, 3L))
})
