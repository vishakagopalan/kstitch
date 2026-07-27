test_that("returns correct structure and high first canonical correlation", {

  set.seed(42)
  n <- 200; k <- 3
  Z           <- matrix(rnorm(n * k), n, k)
  shape_mat   <- scale(Z %*% matrix(rnorm(k * 8), k, 8) + matrix(rnorm(n * 8, sd = 0.3), n, 8))
  nmf_mat     <- scale(Z %*% matrix(rnorm(k * 6), k, 6) + matrix(rnorm(n * 6, sd = 0.3), n, 6))
  cell_names  <- paste0("cell_", seq_len(n))
  rownames(shape_mat) <- cell_names
  rownames(nmf_mat)   <- cell_names
  colnames(shape_mat) <- paste0("ShapePC", seq_len(8))
  colnames(nmf_mat)   <- paste0("Factor",  seq_len(6))

  obj <- make_seurat_stub(cell_names)
  res <- link_shape_and_factors(obj, nmf_mat, shape_mat)

  expect_true(is.list(res))
  expect_true(all(c("CC_Corr_Coefs", "CSP_Scores", "CEP_Scores",
                    "CSP_Vectors",   "CEP_Vectors",
                    "Shape_Corr_With_CSP", "Exp_Corr_With_CEP",
                    "Misc_CCA") %in% names(res)))

  k_out <- length(res$CC_Corr_Coefs)
  expect_equal(nrow(res$CSP_Scores), n)
  expect_equal(ncol(res$CSP_Scores), k_out)
  expect_true(all(res$CC_Corr_Coefs >= 0 & res$CC_Corr_Coefs <= 1))
  expect_true(res$CC_Corr_Coefs[1] > 0.7)
  expect_true(all(grepl("^CSP", colnames(res$CSP_Scores))))
  expect_true(all(grepl("^CEP", colnames(res$CEP_Scores))))
  expect_false(is.null(res$Shape_Corr_With_CSP))
  expect_false(is.null(res$Exp_Corr_With_CEP))
  # No group / is_groupwise fields on flat result
  expect_false("group" %in% names(res))
  expect_false("is_groupwise" %in% names(res))
})

test_that("errors when fewer than min_cells are available", {

  set.seed(3)
  n         <- 50
  shape_mat <- matrix(rnorm(n * 5), n, 5,
                      dimnames = list(paste0("cell_", 1:n), paste0("ShapePC", 1:5)))
  nmf_mat   <- matrix(rnorm(n * 4), n, 4,
                      dimnames = list(paste0("cell_", 1:n), paste0("Factor",  1:4)))
  obj <- make_seurat_stub(paste0("cell_", 1:n))

  expect_error(
    link_shape_and_factors(obj, nmf_mat, shape_mat),
    "minimum is"
  )
})

test_that("errors when no cells remain after intersection", {

  set.seed(4)
  shape_mat <- matrix(rnorm(200 * 5), 200, 5,
                      dimnames = list(paste0("shape_cell_", 1:200), paste0("ShapePC", 1:5)))
  nmf_mat   <- matrix(rnorm(200 * 4), 200, 4,
                      dimnames = list(paste0("nmf_cell_",   1:200), paste0("Factor",  1:4)))
  obj <- make_seurat_stub(paste0("obj_cell_", 1:200))

  expect_error(
    link_shape_and_factors(obj, nmf_mat, shape_mat),
    "No cells remain"
  )
})

test_that("return_results = FALSE writes a single RDS and returns its path", {

  set.seed(3)
  n         <- 300
  shape_mat <- matrix(rnorm(n * 5), n, 5,
                      dimnames = list(paste0("cell_", 1:n), paste0("ShapePC", 1:5)))
  nmf_mat   <- matrix(rnorm(n * 4), n, 4,
                      dimnames = list(paste0("cell_", 1:n), paste0("Factor",  1:4)))
  obj     <- make_seurat_stub(paste0("cell_", 1:n))
  out_dir <- file.path(tempdir(), paste0("test_cca_", kstitch:::.random_id()))

  res <- link_shape_and_factors(obj, nmf_mat, shape_mat,
                                return_results = FALSE, output_dir = out_dir)

  expect_named(res, "output_path")
  expect_true(file.exists(res$output_path))
  expect_equal(basename(res$output_path), "all.rds")

  loaded <- load_kstitch_results(res$output_path, type = "cca")
  expect_true(is.list(loaded))
  expect_true("CSP_Scores" %in% names(loaded))

  unlink(out_dir, recursive = TRUE)
})

test_that("NA rows in shape_mat are dropped and a message is emitted", {

  set.seed(5)
  n         <- 250
  shape_mat <- matrix(rnorm(n * 5), n, 5,
                      dimnames = list(paste0("cell_", 1:n), paste0("ShapePC", 1:5)))
  nmf_mat   <- matrix(rnorm(n * 4), n, 4,
                      dimnames = list(paste0("cell_", 1:n), paste0("Factor",  1:4)))
  shape_mat[1:10, 1] <- NA
  obj <- make_seurat_stub(paste0("cell_", 1:n))

  expect_message(
    res <- link_shape_and_factors(obj, nmf_mat, shape_mat),
    "Dropping"
  )
  expect_equal(nrow(res$CSP_Scores), n - 10L)
})
