test_that("pooled mode returns correct structure and high first canonical correlation", {

  set.seed(42)
  n <- 200
  k <- 3

  # Planted shared latent structure
  Z         <- matrix(rnorm(n * k), n, k)
  shape_mat <- Z %*% matrix(rnorm(k * 8), k, 8) + matrix(rnorm(n * 8, sd = 0.3), n, 8)
  nmf_mat   <- Z %*% matrix(rnorm(k * 6), k, 6) + matrix(rnorm(n * 6, sd = 0.3), n, 6)

  cell_names        <- paste0("cell_", seq_len(n))
  rownames(shape_mat) <- cell_names
  rownames(nmf_mat)   <- cell_names
  colnames(shape_mat) <- paste0("ShapePC", seq_len(8))
  colnames(nmf_mat)   <- paste0("Factor",  seq_len(6))

  # Minimal Seurat-like stub: only @meta.data rownames needed; the function
  # uses rownames(obj@meta.data) rather than colnames(obj) to avoid
  # requiring the Seurat S3 method in tests.
  obj <- list(meta.data = data.frame(row.names = cell_names))

  res <- link_shape_and_factors(obj, nmf_mat, shape_mat, group.by = NULL)

  expect_named(res, "all")
  grp <- res[["all"]]

  # Expected fields
  expect_true(all(c("CC_Corr_Coefs", "CSP_Scores", "CEP_Scores",
                    "CSP_Vectors", "CEP_Vectors",
                    "CSP_Self_Correlations", "CEP_Self_Correlations",
                    "Misc_CCA") %in% names(grp)))

  # Dimensions
  k_out <- length(grp$CC_Corr_Coefs)
  expect_equal(nrow(grp$CSP_Scores), n)
  expect_equal(ncol(grp$CSP_Scores), k_out)
  expect_equal(nrow(grp$CEP_Scores), n)
  expect_equal(ncol(grp$CEP_Scores), k_out)

  # Canonical correlations valid
  expect_true(all(grp$CC_Corr_Coefs >= 0 & grp$CC_Corr_Coefs <= 1))
  expect_true(grp$CC_Corr_Coefs[1] > 0.7)

  # Column names use CSP/CEP terminology
  expect_true(all(grepl("^CSP", colnames(grp$CSP_Scores))))
  expect_true(all(grepl("^CEP", colnames(grp$CEP_Scores))))

  # Self-correlations are not NULL (bug that existed in original)
  expect_false(is.null(grp$CSP_Self_Correlations))
  expect_false(is.null(grp$CEP_Self_Correlations))
})


test_that("group.by splits cells and runs CCA per group", {

  set.seed(7)
  n_a <- 150; n_b <- 180
  n   <- n_a + n_b

  make_block <- function(n, k_shape = 8, k_nmf = 6) {
    Z         <- matrix(rnorm(n * 3), n, 3)
    shape_mat <- Z %*% matrix(rnorm(3 * k_shape), 3, k_shape) + matrix(rnorm(n * k_shape, sd = 0.3), n, k_shape)
    nmf_mat   <- Z %*% matrix(rnorm(3 * k_nmf),   3, k_nmf)   + matrix(rnorm(n * k_nmf,   sd = 0.3), n, k_nmf)
    list(shape = shape_mat, nmf = nmf_mat)
  }

  blk_a <- make_block(n_a); blk_b <- make_block(n_b)

  shape_mat <- rbind(blk_a$shape, blk_b$shape)
  nmf_mat   <- rbind(blk_a$nmf,   blk_b$nmf)

  cell_names          <- paste0("cell_", seq_len(n))
  rownames(shape_mat) <- cell_names
  rownames(nmf_mat)   <- cell_names
  colnames(shape_mat) <- paste0("ShapePC", seq_len(8))
  colnames(nmf_mat)   <- paste0("Factor",  seq_len(6))

  meta <- data.frame(
    celltype = c(rep("TypeA", n_a), rep("TypeB", n_b)),
    row.names = cell_names
  )

  obj <- list(meta.data = meta)

  res <- link_shape_and_factors(obj, nmf_mat, shape_mat, group.by = "celltype")

  expect_named(res, c("TypeA", "TypeB"), ignore.order = TRUE)

  # Each group has the right number of cells
  expect_equal(nrow(res[["TypeA"]]$CSP_Scores), n_a)
  expect_equal(nrow(res[["TypeB"]]$CSP_Scores), n_b)

  # Both groups show high first canonical correlation (planted structure)
  expect_true(res[["TypeA"]]$CC_Corr_Coefs[1] > 0.7)
  expect_true(res[["TypeB"]]$CC_Corr_Coefs[1] > 0.7)
})


test_that("groups below min_cells are skipped", {

  set.seed(3)
  n <- 50  # below default min_cells = 100

  shape_mat <- matrix(rnorm(n * 5), n, 5)
  nmf_mat   <- matrix(rnorm(n * 4), n, 4)
  cell_names          <- paste0("cell_", seq_len(n))
  rownames(shape_mat) <- cell_names
  rownames(nmf_mat)   <- cell_names

  obj <- list(meta.data = data.frame(row.names = cell_names))

  expect_warning(
    res <- link_shape_and_factors(obj, nmf_mat, shape_mat),
    "empty list"
  )
  expect_equal(length(res), 0)
})
