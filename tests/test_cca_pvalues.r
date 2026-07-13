test_that("cca_pvalues returns correct structure", {

  set.seed(1)
  n <- 150; k <- 3
  Z         <- matrix(rnorm(n * k), n, k)
  X         <- scale(Z %*% matrix(rnorm(k * 8), k, 8) + matrix(rnorm(n * 8, sd = 0.3), n, 8))
  Y         <- scale(Z %*% matrix(rnorm(k * 6), k, 6) + matrix(rnorm(n * 6, sd = 0.3), n, 6))

  res <- cca_pvalues(X, Y, K = 3, nperm = 99, seed = 42, verbose = FALSE)

  expect_named(res, c("pvalues", "raw_p", "obs_cancor", "counts", "method", "meta"))
  expect_length(res$pvalues, 3)
  expect_length(res$obs_cancor, 3)
  expect_true(all(res$pvalues >= 0 & res$pvalues <= 1))
  expect_true(all(res$obs_cancor >= 0 & res$obs_cancor <= 1, na.rm = TRUE))
  expect_true(!is.unsorted(res$obs_cancor[!is.na(res$obs_cancor)], strictly = FALSE) == FALSE ||
                res$obs_cancor[1] >= res$obs_cancor[2])  # deflated cors non-increasing
  expect_equal(names(res$pvalues), c("CC1", "CC2", "CC3"))
})


test_that("stepdown p-values are monotone non-decreasing", {

  set.seed(2)
  n <- 150; k <- 3
  Z <- matrix(rnorm(n * k), n, k)
  X <- scale(Z %*% matrix(rnorm(k * 8), k, 8) + matrix(rnorm(n * 8, sd = 0.4), n, 8))
  Y <- scale(Z %*% matrix(rnorm(k * 6), k, 6) + matrix(rnorm(n * 6, sd = 0.4), n, 6))

  res <- cca_pvalues(X, Y, K = 4, nperm = 99, method = "stepdown", seed = 7, verbose = FALSE)

  expect_equal(res$method, "stepdown")
  # cummax property: each p-value >= the previous
  expect_true(all(diff(res$pvalues) >= 0))
})


test_that("planted signal produces small p-values for first component", {

  set.seed(3)
  n <- 200
  # Strong shared signal in first component only
  z1        <- rnorm(n)
  X         <- scale(cbind(z1 + rnorm(n, sd = 0.1), matrix(rnorm(n * 7), n, 7)))
  Y         <- scale(cbind(z1 + rnorm(n, sd = 0.1), matrix(rnorm(n * 5), n, 5)))

  res <- cca_pvalues(X, Y, K = 2, nperm = 199, seed = 99, verbose = FALSE)

  expect_lt(res$pvalues[["CC1"]], 0.05)
})


test_that("null data produces large p-values", {

  set.seed(4)
  n <- 150
  X <- scale(matrix(rnorm(n * 6), n, 6))
  Y <- scale(matrix(rnorm(n * 5), n, 5))

  res <- cca_pvalues(X, Y, K = 2, nperm = 199, seed = 11, verbose = FALSE)

  # Under the null at least CC1 should not be significant at 0.05
  expect_gt(res$pvalues[["CC1"]], 0.05)
})


test_that("block permutation runs without error and respects block structure", {

  set.seed(5)
  n     <- 180
  block <- rep(1:3, each = 60)
  Z     <- matrix(rnorm(n * 3), n, 3)
  X     <- scale(Z %*% matrix(rnorm(3 * 7), 3, 7) + matrix(rnorm(n * 7, sd = 0.3), n, 7))
  Y     <- scale(Z %*% matrix(rnorm(3 * 5), 3, 5) + matrix(rnorm(n * 5, sd = 0.3), n, 5))

  expect_no_error(
    res <- cca_pvalues(X, Y, K = 2, nperm = 49, blocks = block, seed = 6, verbose = FALSE)
  )
  expect_true(res$meta$blocks)
  expect_length(res$pvalues, 2)
})
