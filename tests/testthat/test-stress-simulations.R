test_that("latent exponential sampler survives repeated small simulated runs", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("SPNNGP_RUN_STRESS"), "true"),
    "Set SPNNGP_RUN_STRESS=true to run repeated sampler stress tests."
  )

  for (seed in 101:125) {
    set.seed(seed)
    dat <- make_gaussian_data(n = 50)

    fit <- spNNGP(
      dat$y ~ dat$x,
      coords = dat$coords,
      method = "latent",
      cov.model = "exponential",
      n.neighbors = 8,
      starting = dat$starting,
      tuning = dat$tuning,
      priors = dat$priors,
      n.samples = 20,
      n.omp.threads = 1,
      verbose = FALSE
    )

    expect_equal(dim(fit$p.w.samples), c(50L, 20L))
    expect_true(all(is.finite(fit$p.theta.samples)))
  }
})
