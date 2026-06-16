test_that("run_cca recovers canonical correlation structure", {
  set.seed(1)
  n <- 200
  z <- rnorm(n)                                   # shared latent signal
  X <- cbind(z + rnorm(n, sd = 0.3), rnorm(n), rnorm(n))
  Y <- cbind(z + rnorm(n, sd = 0.3), rnorm(n))
  
  out <- run_cca(X, Y)
  
  # correlations live in [0, 1] and are non-increasing
  expect_true(all(out$cor >= 0 & out$cor <= 1))
  expect_true(all(diff(out$cor) <= 1e-8))
  
  # defining property of CCA: the k-th canonical variates correlate at cor[k]
  for (k in seq_along(out$cor)) {
    emp <- abs(stats::cor(out$scores$xscores[, k], out$scores$yscores[, k]))
    expect_equal(emp, out$cor[k], tolerance = 1e-6)
  }
  
  # the shared latent should drive a strong leading canonical correlation
  expect_gt(out$cor[1], 0.7)
})