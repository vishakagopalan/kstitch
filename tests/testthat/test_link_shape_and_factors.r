test_that("pooled mode returns correct structure and high first canonical correlation", {

  set.seed(42)
  n <- 200; k <- 3
  Z         <- matrix(rnorm(n * k), n, k)
  shape_mat <- scale(Z %*% matrix(rnorm(k * 8), k, 8) + matrix(rnorm(n * 8, sd = 0.3), n, 8))
  nmf_mat   <- scale(Z %*% matrix(rnorm(k * 6), k, 6) + matrix(rnorm(n * 6, sd = 0.3), n, 6))
  cell_names        <- paste0("cell_", seq_len(n))
  rownames(shape_mat) <- cell_names
  rownames(nmf_mat)   <- cell_names
  colnames(shape_mat) <- paste0("ShapePC", seq_len(8))
  colnames(nmf_mat)   <- paste0("Factor",  seq_len(6))

  obj <- make_seurat_stub(cell_names)

  res <- link_shape_and_factors(obj, nmf_mat, shape_mat, group.by = NULL)

  expect_named(res, "all")
  grp <- res[["all"]]

  expect_true(all(c("CC_Corr_Coefs", "CSP_Scores", "CEP_Scores",
                    "CSP_Vectors", "CEP_Vectors",
                    "Shape_Corr_With_CSP", "Exp_Corr_With_CEP",
                    "Misc_CCA") %in% names(grp)))

  k_out <- length(grp$CC_Corr_Coefs)
  expect_equal(nrow(grp$CSP_Scores), n)
  expect_equal(ncol(grp$CSP_Scores), k_out)
  expect_true(all(grp$CC_Corr_Coefs >= 0 & grp$CC_Corr_Coefs <= 1))
  expect_true(grp$CC_Corr_Coefs[1] > 0.7)
  expect_true(all(grepl("^CSP", colnames(grp$CSP_Scores))))
  expect_true(all(grepl("^CEP", colnames(grp$CEP_Scores))))
  expect_false(is.null(grp$Shape_Corr_With_CSP))
  expect_false(is.null(grp$Exp_Corr_With_CEP))
})


test_that("group.by splits cells and runs CCA per group", {

  set.seed(7)
  n_a <- 150; n_b <- 180; n <- n_a + n_b

  make_block <- function(n, k_shape = 8, k_nmf = 6) {
    Z         <- matrix(rnorm(n * 3), n, 3)
    shape_mat <- Z %*% matrix(rnorm(3 * k_shape), 3, k_shape) + matrix(rnorm(n * k_shape, sd = 0.3), n, k_shape)
    nmf_mat   <- Z %*% matrix(rnorm(3 * k_nmf),   3, k_nmf)   + matrix(rnorm(n * k_nmf,   sd = 0.3), n, k_nmf)
    list(shape = shape_mat, nmf = nmf_mat)
  }

  blk_a <- make_block(n_a); blk_b <- make_block(n_b)
  shape_mat <- scale(rbind(blk_a$shape, blk_b$shape))
  nmf_mat   <- scale(rbind(blk_a$nmf,   blk_b$nmf))
  cell_names          <- paste0("cell_", seq_len(n))
  rownames(shape_mat) <- cell_names
  rownames(nmf_mat)   <- cell_names
  colnames(shape_mat) <- paste0("ShapePC", seq_len(8))
  colnames(nmf_mat)   <- paste0("Factor",  seq_len(6))

  obj <- make_seurat_stub(cell_names, meta_cols = list(
    celltype = c(rep("TypeA", n_a), rep("TypeB", n_b))
  ))

  res <- link_shape_and_factors(obj, nmf_mat, shape_mat, group.by = "celltype")

  expect_named(res, c("TypeA", "TypeB"), ignore.order = TRUE)
  expect_equal(nrow(res[["TypeA"]]$CSP_Scores), n_a)
  expect_equal(nrow(res[["TypeB"]]$CSP_Scores), n_b)
  expect_true(res[["TypeA"]]$CC_Corr_Coefs[1] > 0.7)
  expect_true(res[["TypeB"]]$CC_Corr_Coefs[1] > 0.7)
})


test_that("groups below min_cells are skipped", {

  set.seed(3)
  n <- 50
  shape_mat <- matrix(rnorm(n * 5), n, 5)
  nmf_mat   <- matrix(rnorm(n * 4), n, 4)
  cell_names          <- paste0("cell_", seq_len(n))
  rownames(shape_mat) <- cell_names
  rownames(nmf_mat)   <- cell_names

  obj <- make_seurat_stub(cell_names)

  expect_warning(
    res <- link_shape_and_factors(obj, nmf_mat, shape_mat),
    "empty list"
  )
  expect_equal(length(res), 0)
})

test_that("link_shape_and_factors attaches group and is_groupwise — single group", {
  set.seed(1)
  n <- 200
  shape_mat <- matrix(rnorm(n * 5), n, 5, dimnames = list(paste0("cell_", 1:n), paste0("ShapePC", 1:5)))
  nmf_mat   <- matrix(rnorm(n * 4), n, 4, dimnames = list(paste0("cell_", 1:n), paste0("Factor",  1:4)))
  obj <- make_seurat_stub(paste0("cell_", 1:n))

  res <- link_shape_and_factors(obj, nmf_mat, shape_mat)

  expect_equal(res$all$group, "all")
  expect_false(res$all$is_groupwise)
})

test_that("link_shape_and_factors attaches group and is_groupwise — grouped", {
  set.seed(2)
  n <- 300
  shape_mat <- matrix(rnorm(n * 5), n, 5, dimnames = list(paste0("cell_", 1:n), paste0("ShapePC", 1:5)))
  nmf_mat   <- matrix(rnorm(n * 4), n, 4, dimnames = list(paste0("cell_", 1:n), paste0("Factor",  1:4)))
  obj <- make_seurat_stub(paste0("cell_", 1:n),
                          meta_cols = list(celltype = rep(c("TypeA", "TypeB"), each = n / 2)))

  res <- link_shape_and_factors(obj, nmf_mat, shape_mat, group.by = "celltype")

  expect_equal(res$TypeA$group, "TypeA")
  expect_equal(res$TypeB$group, "TypeB")
  expect_true(res$TypeA$is_groupwise)
  expect_true(res$TypeB$is_groupwise)
})

test_that("link_shape_and_factors return_results = FALSE writes RDS files and returns paths", {
  set.seed(3)
  n <- 300
  shape_mat <- matrix(rnorm(n * 5), n, 5, dimnames = list(paste0("cell_", 1:n), paste0("ShapePC", 1:5)))
  nmf_mat   <- matrix(rnorm(n * 4), n, 4, dimnames = list(paste0("cell_", 1:n), paste0("Factor",  1:4)))
  obj <- make_seurat_stub(paste0("cell_", 1:n),
                          meta_cols = list(celltype = rep(c("TypeA", "TypeB"), each = n / 2)))
  out_dir <- file.path(tempdir(), paste0("test_cca_", kstitch:::.random_id()))

  res <- link_shape_and_factors(obj, nmf_mat, shape_mat, group.by = "celltype",
                                return_results = FALSE, output_dir = out_dir)

  expect_named(res, "output_paths")
  expect_named(res$output_paths, c("TypeA", "TypeB"), ignore.order = TRUE)
  expect_true(all(file.exists(res$output_paths)))

  loaded <- load_kstitch_results(res$output_paths[["TypeA"]], type = "cca")
  expect_true(is.list(loaded))
  expect_true("CSP_Scores" %in% names(loaded))

  unlink(out_dir, recursive = TRUE)
})
