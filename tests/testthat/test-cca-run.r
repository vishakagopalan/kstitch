test_that("run_cca recovers known canonical correlation structure", {

  set.seed(1)
  n <- 300; k <- 3

  # Shared latent factors with independent, small per-column noise
  Z <- matrix(rnorm(n * k), n, k)
  X <- scale(Z %*% matrix(rnorm(k * 8), k, 8) + matrix(rnorm(n * 8, sd = 0.2), n, 8))
  Y <- scale(Z %*% matrix(rnorm(k * 6), k, 6) + matrix(rnorm(n * 6, sd = 0.2), n, 6))
  rownames(X) <- rownames(Y) <- paste0("cell_", seq_len(n))
  colnames(X) <- paste0("ShapePC", seq_len(8))
  colnames(Y) <- paste0("Factor",  seq_len(6))

  res <- run_cca(X, Y, scale = FALSE)

  K <- min(ncol(X), ncol(Y))
  expect_length(res$cor, K)

  # All correlations must lie in the valid, clamped range.
  expect_true(all(res$cor >= 0 & res$cor <= 1))

  # Canonical correlations should be non-increasing (this is a property of
  # SVD singular values, sv$d, which are returned in descending order).
  expect_true(all(diff(res$cor) <= 1e-8))

  # With a 3-factor shared latent structure and low noise, the first 3
  # canonical correlations should be substantially larger than the rest.
  expect_true(all(res$cor[1:3] > 0.8))
  if (K > 3) {
    expect_true(all(res$cor[1:3] > res$cor[4:K] + 0.2))
  }

  # Structural checks on return shape.
  expect_named(res, c("cor", "xcoef", "ycoef", "scores"))
  expect_equal(dim(res$xcoef), c(ncol(X), K))
  expect_equal(dim(res$ycoef), c(ncol(Y), K))
  expect_equal(rownames(res$xcoef), colnames(X))
  expect_equal(rownames(res$ycoef), colnames(Y))
  expect_equal(dim(res$scores$xscores), c(n, K))
  expect_equal(dim(res$scores$yscores), c(n, K))
  expect_equal(rownames(res$scores$xscores), rownames(X))
})


test_that("run_cca produces zero-ish correlations under the null", {

  set.seed(2)
  n <- 300
  X <- scale(matrix(rnorm(n * 6), n, 6))
  Y <- scale(matrix(rnorm(n * 5), n, 5))
  rownames(X) <- rownames(Y) <- paste0("cell_", seq_len(n))

  res <- run_cca(X, Y, scale = FALSE)

  # Under the null, canonical correlations are inflated upward by
  # overfitting (small n relative to K), so this is a loose upper bound,
  # not a tight one -- just confirming no structure is spuriously large.
  expect_true(all(res$cor < 0.5))
})


test_that("run_cca errors on mismatched row counts", {

  X <- matrix(rnorm(100 * 5), 100, 5)
  Y <- matrix(rnorm(90  * 4),  90, 4)

  expect_error(run_cca(X, Y), "same number of rows")
})


test_that("run_cca does not error or produce NaN on near-singular input", {

  set.seed(3)
  n <- 100
  # Two nearly identical (highly collinear) columns in X to stress the
  # eigenvalue floor in inv_sqrt().
  base <- rnorm(n)
  X <- cbind(base, base + rnorm(n, sd = 1e-6), matrix(rnorm(n * 4), n, 4))
  Y <- matrix(rnorm(n * 4), n, 4)
  rownames(X) <- rownames(Y) <- paste0("cell_", seq_len(n))

  res <- run_cca(X, Y, scale = TRUE)

  expect_false(anyNA(res$cor))
  expect_false(any(is.nan(res$cor)))
  expect_true(all(is.finite(res$xcoef)))
  expect_true(all(is.finite(res$ycoef)))
})
