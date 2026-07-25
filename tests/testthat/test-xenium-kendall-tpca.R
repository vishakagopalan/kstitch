# ---- helpers ----------------------------------------------------------------

.load_fixture_cell_ids <- function() {
  fixture_dir <- system.file(
    "extdata", "xenium_test", "tpca_fixtures", "Unit_Test_Keratinocytes_Nuclei",
    package = "kstitch"
  )
  if (!nzchar(fixture_dir)) return(NULL)
  result <- load_kstitch_results(fixture_dir, type = "tpca")
  rownames(result$all$TPCA_Embedding)
}

.split_cell_ids <- function(cell_ids) {
  mid <- floor(length(cell_ids) / 2)
  list(
    groupA = cell_ids[seq_len(mid)],
    groupB = cell_ids[seq(mid + 1L, length(cell_ids))]
  )
}

.boundary_path <- function() {
  system.file(
    "extdata", "xenium_test", "Test_Keratinocyte_Nuclear_Coordinates.parquet",
    package = "kstitch"
  )
}

.fixture_dir <- function() {
  system.file(
    "extdata", "xenium_test", "tpca_fixtures", "Unit_Test_Keratinocytes_Nuclei",
    package = "kstitch"
  )
}

.std_run_args <- list(
  use_parallel = FALSE,
  x_col        = "x",
  y_col        = "y",
  cell_id_col  = "cell"
)

test_that("run_tpca load_pre_shape = FALSE sets pre_shape_embedding to NULL", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  result <- do.call(run_tpca, c(
    list(boundary_parquet_path = .boundary_path(),
         output_dir            = NULL,
         load_pre_shape        = FALSE),
    .std_run_args
  ))

  expect_null(result$all$Info$pre_shape_embedding)
  expect_true(is.matrix(result$all$TPCA_Embedding))
  expect_true(all(c("variances", "v_matrix", "frechet_mean") %in%
                    names(result$all$Info)))
})

test_that("run_tpca groupwise load_pre_shape = FALSE sets pre_shape_embedding to NULL for all groups", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  cell_ids <- .load_fixture_cell_ids()
  skip_if(is.null(cell_ids) || length(cell_ids) < 4L, "not enough fixture cells")

  splits      <- .split_cell_ids(cell_ids)
  cell_groups <- c(
    stats::setNames(rep("groupA", length(splits$groupA)), splits$groupA),
    stats::setNames(rep("groupB", length(splits$groupB)), splits$groupB)
  )

  result <- do.call(run_tpca, c(
    list(boundary_parquet_path = .boundary_path(),
         output_dir            = withr::local_tempdir(),
         cell_groups           = cell_groups,
         load_pre_shape        = FALSE),
    .std_run_args
  ))

  for (grp in c("groupA", "groupB")) {
    expect_null(result[[grp]]$Info$pre_shape_embedding,
                label = paste(grp, "pre_shape_embedding is NULL"))
    expect_true(is.matrix(result[[grp]]$TPCA_Embedding),
                label = paste(grp, "TPCA_Embedding present"))
  }
})

test_that("load_kstitch_results load_pre_shape = FALSE sets pre_shape_embedding to NULL", {
  skip_on_ci()
  skip_if(!nzchar(.fixture_dir()), "TPCA fixture not found")

  result <- load_kstitch_results(.fixture_dir(), type = "tpca",
                                 load_pre_shape = FALSE)

  expect_null(result$all$Info$pre_shape_embedding)
  expect_true(is.matrix(result$all$TPCA_Embedding))
})

# ---- single-group: temp dir -------------------------------------------------

test_that("run_tpca defaults to a temp directory and cleans it up", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  result <- do.call(run_tpca, c(
    list(boundary_parquet_path = .boundary_path(), output_dir = NULL),
    .std_run_args
  ))

  expect_named(result, "all")
  expect_null(result$all$output_dir)
  expect_equal(result$all$group, "all")
  expect_false(result$all$is_groupwise)
  expect_true(is.matrix(result$all$TPCA_Embedding))
  expect_true(nrow(result$all$TPCA_Embedding) > 0)
  expect_true(all(c("variances", "v_matrix", "frechet_mean") %in% names(result$all$Info)))
})

# ---- single-group: explicit output_dir --------------------------------------

test_that("run_tpca with an explicit output_dir leaves files on disk", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  out_dir <- withr::local_tempdir()
  result  <- do.call(run_tpca, c(
    list(boundary_parquet_path = .boundary_path(), output_dir = out_dir),
    .std_run_args
  ))

  expect_named(result, "all")
  expect_equal(result$all$output_dir, out_dir)
  expect_equal(result$all$group, "all")
  expect_false(result$all$is_groupwise)
  expect_true(file.exists(file.path(out_dir, "TPCA_Info.h5")))
  expect_true(file.exists(file.path(out_dir, "Shape_Metadata.csv.gz")))
})

# ---- single-group: return_results = FALSE -----------------------------------

test_that("run_tpca return_results = FALSE returns paths and preserves files", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  out_dir <- withr::local_tempdir()
  result  <- do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = out_dir,
      return_results        = FALSE
    ),
    .std_run_args
  ))

  expect_named(result, c("output_paths", "fields"))
  expect_true(is.character(result$output_paths))
  expect_true(length(result$output_paths) == 1L)
  expect_true(file.exists(file.path(result$output_paths[[1]], "TPCA_Info.h5")))
  expect_true(file.exists(file.path(result$output_paths[[1]], "Shape_Metadata.csv.gz")))
  expect_null(result$TPCA_Embedding)
})

test_that("run_tpca return_results = FALSE with NULL output_dir preserves temp files", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  result <- do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = NULL,
      return_results        = FALSE
    ),
    .std_run_args
  ))

  expect_named(result, c("output_paths", "fields"))
  expect_true(file.exists(file.path(result$output_paths[[1]], "TPCA_Info.h5")))
})

test_that("load_kstitch_results type='tpca' loads deferred results correctly", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  out_dir <- withr::local_tempdir()
  deferred <- do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = out_dir,
      return_results        = FALSE
    ),
    .std_run_args
  ))

  loaded <- load_kstitch_results(deferred$output_paths[[1]], type = "tpca")

  expect_true(is.matrix(loaded$all$TPCA_Embedding))
  expect_true(nrow(loaded$all$TPCA_Embedding) > 0)
  expect_true(all(c("variances", "v_matrix", "frechet_mean") %in% names(loaded$all$Info)))
  expect_s3_class(loaded$all$Metadata, "data.frame")
})

# ---- single-group: known output ---------------------------------------------

test_that("run_tpca reproduces known output on fixture contours", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")
  skip_if(!nzchar(.fixture_dir()),   "expected tpca fixture not found")

  tmp_out  <- withr::local_tempdir()
  result   <- do.call(run_tpca, c(
    list(boundary_parquet_path = .boundary_path(), output_dir = tmp_out),
    .std_run_args
  ))
  expected <- load_kstitch_results(.fixture_dir(), type = "tpca")$all

  expect_setequal(rownames(result$all$TPCA_Embedding), rownames(expected$TPCA_Embedding))
  expect_equal(ncol(result$all$TPCA_Embedding), ncol(expected$TPCA_Embedding))

  common_ids  <- rownames(expected$TPCA_Embedding)
  res_aligned <- result$all$TPCA_Embedding[common_ids, , drop = FALSE]
  exp_aligned <- expected$TPCA_Embedding[common_ids, , drop = FALSE]

  expect_equal(result$all$Info$variances, expected$Info$variances, tolerance = 1e-4)

  vars       <- expected$Info$variances
  rel_gap    <- abs(diff(vars)) / (vars[-length(vars)] + 1e-12)
  degenerate <- c(FALSE, rel_gap < 1e-3) | c(rel_gap < 1e-3, FALSE)
  negligible <- vars / sum(vars) < 1e-3
  stable_pcs <- which(!degenerate & !negligible)
  skip_if(length(stable_pcs) == 0, "no numerically stable PCs to compare")

  for (j in stable_pcs) {
    diff_same <- mean(abs(res_aligned[, j] - exp_aligned[, j]))
    diff_flip <- mean(abs(res_aligned[, j] + exp_aligned[, j]))
    expect_true(
      min(diff_same, diff_flip) < 1e-4,
      label = paste0("PC", j, " sign-adjusted difference")
    )
  }

  proj_result   <- res_aligned %*% t(res_aligned)
  proj_expected <- exp_aligned %*% t(exp_aligned)
  expect_true(mean(abs(proj_result - proj_expected)) < 1e-3)
})

# ---- single-group: Frechet mean finiteness ----------------------------------

test_that("Frechet mean is finite on fixture contours", {
  skip_on_ci()
  skip_if(!nzchar(.fixture_dir()), "expected tpca fixture not found")

  result <- load_kstitch_results(.fixture_dir(), type = "tpca")
  expect_true(all(is.finite(result$all$Info$frechet_mean)))
})

# ---- single-group: parallel vs serial ---------------------------------------

test_that("frechet accepts a warm-start mu and converges faster than cold start", {
  skip_on_ci()
  skip_if(!nzchar(.fixture_dir()), "TPCA fixture not found")

  np   <- reticulate::import("numpy", convert = FALSE)
  ktpy <- reticulate::import_from_path(
    "kendall_tpca",
    path = system.file("python", package = "kstitch")
  )

  pre_shape_path <- file.path(.fixture_dir(), "Pre_Shape_Space_Embedding.h5")
  skip_if(!file.exists(pre_shape_path), "pre-shape fixture not found")

  raw <- rhdf5::h5read(pre_shape_path, "pre_shape_space_embedding")
  X_r <- aperm(raw, c(3L, 2L, 1L))
  X   <- reticulate::r_to_py(X_r)

  mu_true <- ktpy$frechet(X, eta = 1, tol = 1e-4, max_iter = 1000L)

  L         <- dim(X_r)[3]
  set.seed(42L)
  mu_rand_r <- matrix(rnorm(2L * L), nrow = 2L, ncol = L)
  mu_rand_r <- mu_rand_r - rowMeans(mu_rand_r)
  mu_rand_r <- mu_rand_r / norm(mu_rand_r, type = "F")
  mu_rand   <- reticulate::r_to_py(mu_rand_r)

  mu_warm_hist <- ktpy$frechet(
    X, eta = 1, tol = 1e-4, max_iter = 1000L,
    store_history = TRUE, mu = mu_true
  )
  iters_warm <- length(mu_warm_hist[[2L]])

  mu_cold_hist <- ktpy$frechet(
    X, eta = 1, tol = 1e-4, max_iter = 1000L,
    store_history = TRUE
  )
  iters_cold <- length(mu_cold_hist[[2L]])

  expect_true(iters_warm < iters_cold,
              label = "warm start from true mean converges faster than cold start")

  mu_rand_result <- ktpy$frechet(
    X, eta = 1, tol = 1e-4, max_iter = 1000L, mu = mu_rand
  )
  expect_true(all(is.finite(reticulate::py_to_r(mu_rand_result))),
              label = "frechet with random mu warm start returns finite result")
})

test_that("run_kendall_tpca accepts a warm-start mu and matches cold start output", {
  skip_on_ci()
  skip_if(!nzchar(.fixture_dir()), "TPCA fixture not found")

  ktpy <- reticulate::import_from_path(
    "kendall_tpca",
    path = system.file("python", package = "kstitch")
  )

  pre_shape_path <- file.path(.fixture_dir(), "Pre_Shape_Space_Embedding.h5")
  skip_if(!file.exists(pre_shape_path), "pre-shape fixture not found")

  raw <- rhdf5::h5read(pre_shape_path, "pre_shape_space_embedding")
  X   <- reticulate::r_to_py(aperm(raw, c(3L, 2L, 1L)))

  mu_true <- ktpy$frechet(X, eta = 1, tol = 1e-4, max_iter = 1000L)

  out_cold <- withr::local_tempdir()
  out_warm <- withr::local_tempdir()

  file.copy(
    file.path(.fixture_dir(), "Shape_Metadata.csv.gz"),
    c(file.path(out_cold, "Shape_Metadata.csv.gz"),
      file.path(out_warm, "Shape_Metadata.csv.gz"))
  )

  ktpy$run_kendall_tpca(
    pre_shape_input_dir   = .fixture_dir(),
    output_dir            = out_cold,
    eta                   = 1,
    frechet_mean_tol      = 1e-4,
    max_frechet_mean_iter = 1000L
  )

  ktpy$run_kendall_tpca(
    pre_shape_input_dir   = .fixture_dir(),
    output_dir            = out_warm,
    eta                   = 1,
    frechet_mean_tol      = 1e-4,
    max_frechet_mean_iter = 1000L,
    mu                    = mu_true
  )

  res_cold <- load_kstitch_results(out_cold, type = "tpca")$all
  res_warm <- load_kstitch_results(out_warm, type = "tpca")$all

  common_ids <- intersect(rownames(res_cold$TPCA_Embedding),
                          rownames(res_warm$TPCA_Embedding))
  expect_true(length(common_ids) > 0)

  for (j in seq_len(ncol(res_cold$TPCA_Embedding))) {
    a <- res_cold$TPCA_Embedding[common_ids, j]
    b <- res_warm$TPCA_Embedding[common_ids, j]
    expect_true(
      min(mean(abs(a - b)), mean(abs(a + b))) < 1e-3,
      label = paste0("PC", j, " cold vs warm start agree")
    )
  }

  expect_equal(res_cold$Info$variances, res_warm$Info$variances, tolerance = 1e-3)
})

test_that("top PC is robust to use_parallel = TRUE vs FALSE", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  result_serial <- run_tpca(
    boundary_parquet_path = .boundary_path(),
    output_dir            = NULL,
    use_parallel          = FALSE,
    x_col                 = "x",
    y_col                 = "y",
    cell_id_col           = "cell"
  )
  result_parallel <- run_tpca(
    boundary_parquet_path = .boundary_path(),
    output_dir            = NULL,
    use_parallel          = TRUE,
    num_threads           = 2L,
    x_col                 = "x",
    y_col                 = "y",
    cell_id_col           = "cell"
  )

  common_ids <- intersect(
    rownames(result_serial$all$TPCA_Embedding),
    rownames(result_parallel$all$TPCA_Embedding)
  )
  expect_true(length(common_ids) > 0)

  a <- result_serial$all$TPCA_Embedding[common_ids, 1L]
  b <- result_parallel$all$TPCA_Embedding[common_ids, 1L]
  expect_true(min(mean(abs(a - b)), mean(abs(a + b))) < 1e-3)
})

# ---- groupwise: structure ---------------------------------------------------

test_that("run_tpca groupwise returns named list with one result per group", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  cell_ids <- .load_fixture_cell_ids()
  skip_if(is.null(cell_ids) || length(cell_ids) < 4L, "not enough fixture cells")

  splits      <- .split_cell_ids(cell_ids)
  cell_groups <- c(
    stats::setNames(rep("groupA", length(splits$groupA)), splits$groupA),
    stats::setNames(rep("groupB", length(splits$groupB)), splits$groupB)
  )

  result <- do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = withr::local_tempdir(),
      cell_groups           = cell_groups
    ),
    .std_run_args
  ))

  expect_named(result, c("groupA", "groupB"), ignore.order = TRUE)
  for (grp in c("groupA", "groupB")) {
    expect_true(is.matrix(result[[grp]]$TPCA_Embedding))
    expect_true(nrow(result[[grp]]$TPCA_Embedding) > 0)
    expect_true(all(c("variances", "v_matrix", "frechet_mean") %in%
                      names(result[[grp]]$Info)))
    expect_equal(result[[grp]]$contour_type, "cell")
    expect_equal(result[[grp]]$group, grp)
    expect_true(result[[grp]]$is_groupwise)
  }
})

test_that("run_tpca groupwise cells do not bleed across groups", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  cell_ids <- .load_fixture_cell_ids()
  skip_if(is.null(cell_ids) || length(cell_ids) < 4L, "not enough fixture cells")

  splits      <- .split_cell_ids(cell_ids)
  cell_groups <- c(
    stats::setNames(rep("groupA", length(splits$groupA)), splits$groupA),
    stats::setNames(rep("groupB", length(splits$groupB)), splits$groupB)
  )

  result <- do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = withr::local_tempdir(),
      cell_groups           = cell_groups
    ),
    .std_run_args
  ))

  expect_true(length(intersect(
    rownames(result$groupA$TPCA_Embedding),
    rownames(result$groupB$TPCA_Embedding)
  )) == 0L)
  expect_true(all(rownames(result$groupA$TPCA_Embedding) %in% splits$groupA))
  expect_true(all(rownames(result$groupB$TPCA_Embedding) %in% splits$groupB))
})

test_that("run_tpca groupwise writes to correctly named subdirectories", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  cell_ids <- .load_fixture_cell_ids()
  skip_if(is.null(cell_ids) || length(cell_ids) < 4L, "not enough fixture cells")

  splits      <- .split_cell_ids(cell_ids)
  cell_groups <- c(
    stats::setNames(rep("groupA", length(splits$groupA)), splits$groupA),
    stats::setNames(rep("groupB", length(splits$groupB)), splits$groupB)
  )
  out_dir <- withr::local_tempdir()

  do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = out_dir,
      cell_groups           = cell_groups
    ),
    .std_run_args
  ))

  for (grp in c("groupA", "groupB")) {
    sub <- file.path(out_dir, paste0(grp, "_cell"))
    expect_true(dir.exists(sub),                              label = paste(grp, "subdir exists"))
    expect_true(file.exists(file.path(sub, "TPCA_Info.h5")), label = paste(grp, "TPCA_Info.h5"))
    expect_true(file.exists(file.path(sub, "Shape_Metadata.csv.gz")),
                label = paste(grp, "Shape_Metadata.csv.gz"))
  }
})

test_that("run_tpca groupwise Frechet means are finite for both groups", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  cell_ids <- .load_fixture_cell_ids()
  skip_if(is.null(cell_ids) || length(cell_ids) < 4L, "not enough fixture cells")

  splits      <- .split_cell_ids(cell_ids)
  cell_groups <- c(
    stats::setNames(rep("groupA", length(splits$groupA)), splits$groupA),
    stats::setNames(rep("groupB", length(splits$groupB)), splits$groupB)
  )

  result <- do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = withr::local_tempdir(),
      cell_groups           = cell_groups
    ),
    .std_run_args
  ))

  expect_true(all(is.finite(result$groupA$Info$frechet_mean)))
  expect_true(all(is.finite(result$groupB$Info$frechet_mean)))
})

# ---- groupwise: return_results = FALSE --------------------------------------

test_that("run_tpca groupwise return_results = FALSE returns paths for all groups", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  cell_ids <- .load_fixture_cell_ids()
  skip_if(is.null(cell_ids) || length(cell_ids) < 4L, "not enough fixture cells")

  splits      <- .split_cell_ids(cell_ids)
  cell_groups <- c(
    stats::setNames(rep("groupA", length(splits$groupA)), splits$groupA),
    stats::setNames(rep("groupB", length(splits$groupB)), splits$groupB)
  )
  out_dir <- withr::local_tempdir()

  result <- do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = out_dir,
      cell_groups           = cell_groups,
      return_results        = FALSE
    ),
    .std_run_args
  ))

  expect_named(result, c("output_paths", "fields"))
  expect_named(result$output_paths, c("groupA", "groupB"), ignore.order = TRUE)
  expect_null(result$TPCA_Embedding)

  for (grp in c("groupA", "groupB")) {
    expect_true(file.exists(
      file.path(result$output_paths[[grp]], "TPCA_Info.h5")
    ))
  }
})

test_that("load_kstitch_results loads groupwise deferred results correctly", {
  skip_on_ci()
  skip_if(!nzchar(.boundary_path()), "boundary fixture not found")

  cell_ids <- .load_fixture_cell_ids()
  skip_if(is.null(cell_ids) || length(cell_ids) < 4L, "not enough fixture cells")

  splits      <- .split_cell_ids(cell_ids)
  cell_groups <- c(
    stats::setNames(rep("groupA", length(splits$groupA)), splits$groupA),
    stats::setNames(rep("groupB", length(splits$groupB)), splits$groupB)
  )
  out_dir <- withr::local_tempdir()

  deferred <- do.call(run_tpca, c(
    list(
      boundary_parquet_path = .boundary_path(),
      output_dir            = out_dir,
      cell_groups           = cell_groups,
      return_results        = FALSE
    ),
    .std_run_args
  ))

  for (grp in c("groupA", "groupB")) {
    loaded <- load_kstitch_results(deferred$output_paths[[grp]], type = "tpca")
    expect_true(is.matrix(loaded$all$TPCA_Embedding))
    expect_true(nrow(loaded$all$TPCA_Embedding) > 0)
    expect_true(all(c("variances", "v_matrix", "frechet_mean") %in% names(loaded$all$Info)))
  }
})
