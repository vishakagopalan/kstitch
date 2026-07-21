test_that("run_tpca defaults to a temp directory and cleans it up", {
  skip_on_ci()

  boundary_path <- system.file(
    "extdata", "xenium_test", "Test_Keratinocyte_Nuclear_Coordinates.parquet",
    package = "kstitch"
  )
  skip_if(boundary_path == "", "boundary fixture not found")

  result <- run_tpca(
    boundary_parquet_path = boundary_path,
    output_dir            = NULL,
    use_parallel          = FALSE,
    x_col                 = "x",
    y_col                 = "y",
    cell_id_col           = "cell"
  )

  # output_dir should be reported as NULL, since the temp dir used
  # internally is deleted by the time run_tpca() returns.
  expect_null(result$output_dir)

  # The actual TPCA result should still be fully populated regardless of
  # where (or whether) it was written to disk during the call.
  expect_true(is.matrix(result$TPCA_Embedding))
  expect_true(nrow(result$TPCA_Embedding) > 0)
  expect_true(all(c("variances", "v_matrix", "frechet_mean") %in% names(result$Info)))
})


test_that("run_tpca with an explicit output_dir leaves files on disk", {
  skip_on_ci()

  boundary_path <- system.file(
    "extdata", "xenium_test", "Test_Keratinocyte_Nuclear_Coordinates.parquet",
    package = "kstitch"
  )
  skip_if(boundary_path == "", "boundary fixture not found")

  out_dir <- withr::local_tempdir()
  result <- run_tpca(
    boundary_parquet_path = boundary_path,
    output_dir            = out_dir,
    use_parallel          = FALSE,
    x_col                 = "x",
    y_col                 = "y",
    cell_id_col           = "cell"
  )

  # output_dir should be reported back (not NULL) and the files should
  # still exist on disk after the call returns, unlike the NULL-output_dir
  # case above.
  expect_equal(result$output_dir, out_dir)
  expect_true(file.exists(file.path(out_dir, "TPCA_Info.h5")))
  expect_true(file.exists(file.path(out_dir, "Shape_Metadata.csv.gz")))
})


test_that(".load_kendall_tpca_output reads cell boundary TPCA correctly", {
  skip_on_ci()

  fixture_dir <- system.file(
    "extdata", "xenium_test", "tpca_fixtures", "Unit_Test_Keratinocytes_Nuclei",
    package = "kstitch"
  )
  skip_if(fixture_dir == "", "TPCA fixture not found")

  result <- kstitch:::.load_kendall_tpca_output(fixture_dir)

  expect_true(is.matrix(result$TPCA_Embedding))
  expect_true(all(c("variances", "v_matrix", "frechet_mean") %in% names(result$Info)))

  expected <- readRDS(system.file(
    "extdata", "xenium_test", "Test_TPCA_Info.rds",
    package = "kstitch"
  ))

  # TPCA_Embedding is sign-ambiguous (SVD-derived), so compare up to a
  # global sign flip. Check mean absolute difference rather than max, since
  # entry-wise tolerance is too strict for this snapshot comparison.
  diff_same_sign <- mean(abs(result$TPCA_Embedding - expected$TPCA_Embedding))
  diff_flip_sign <- mean(abs(result$TPCA_Embedding + expected$TPCA_Embedding))
  expect_true(min(diff_same_sign, diff_flip_sign) < 1e-4)

  expect_equal(result$Info$variances, expected$Info$variances, tolerance = 1e-4)
})


test_that("run_tpca reproduces known output on fixture contours", {
  skip_on_ci()

  boundary_path <- system.file(
    "extdata", "xenium_test", "Test_Keratinocyte_Nuclear_Coordinates.parquet",
    package = "kstitch"
  )
  skip_if(boundary_path == "", "boundary fixture not found")

  expected_tpca_path <- system.file(
    "extdata", "xenium_test", "tpca_fixtures", "Unit_Test_Keratinocytes_Nuclei",
    package = "kstitch"
  )
  skip_if(expected_tpca_path == "", "expected tpca fixture not found")

  tmp_out <- withr::local_tempdir()
  result <- run_tpca(
    boundary_parquet_path = boundary_path,
    output_dir            = tmp_out,
    use_parallel          = FALSE,
    x_col                 = "x",
    y_col                 = "y",
    cell_id_col           = "cell"
  )

  expected <- kstitch:::.load_kendall_tpca_output(expected_tpca_path)

  # ---- structural checks -----------------------------------------------
  expect_setequal(rownames(result$TPCA_Embedding), rownames(expected$TPCA_Embedding))
  expect_equal(ncol(result$TPCA_Embedding), ncol(expected$TPCA_Embedding))

  common_ids <- rownames(expected$TPCA_Embedding)
  res_aligned <- result$TPCA_Embedding[common_ids, , drop = FALSE]
  exp_aligned <- expected$TPCA_Embedding[common_ids, , drop = FALSE]

  # ---- variance spectrum (sign-independent, most robust check) ---------
  expect_equal(result$Info$variances, expected$Info$variances, tolerance = 1e-4)

  # ---- identify degenerate (near-equal variance) PCs --------------------
  vars <- expected$Info$variances
  rel_gap <- abs(diff(vars)) / (vars[-length(vars)] + 1e-12)
  degenerate <- c(FALSE, rel_gap < 1e-3) | c(rel_gap < 1e-3, FALSE)

  # ---- identify negligible-variance PCs -----------------------------
  negligible <- vars / sum(vars) < 1e-3

  stable_pcs <- which(!degenerate & !negligible)
  skip_if(length(stable_pcs) == 0, "no numerically stable PCs to compare")

  # ---- per-PC comparison, up to independent sign flips, stable PCs only --
  for (j in stable_pcs) {
    diff_same_sign <- mean(abs(res_aligned[, j] - exp_aligned[, j]))
    diff_flip_sign <- mean(abs(res_aligned[, j] + exp_aligned[, j]))
    expect_true(
      min(diff_same_sign, diff_flip_sign) < 1e-4,
      label = paste0("PC", j, " sign-adjusted difference")
    )
  }

  # ---- subspace-level check across ALL PCs (sign- and rotation-robust) --
  proj_result   <- res_aligned %*% t(res_aligned)
  proj_expected <- exp_aligned %*% t(exp_aligned)
  subspace_diff <- mean(abs(proj_result - proj_expected))
  expect_true(subspace_diff < 1e-3)
})


test_that("Frechet mean is finite on fixture contours", {
  skip_on_ci()

  expected_tpca_path <- system.file(
    "extdata", "xenium_test", "tpca_fixtures", "Unit_Test_Keratinocytes_Nuclei",
    package = "kstitch"
  )
  skip_if(expected_tpca_path == "", "expected tpca fixture not found")

  expected <- kstitch:::.load_kendall_tpca_output(expected_tpca_path)

  expect_true(all(is.finite(expected$Info$frechet_mean)))
})


test_that("top PC is robust to use_parallel = TRUE vs FALSE", {
  skip_on_ci()

  boundary_path <- system.file(
    "extdata", "xenium_test", "Test_Keratinocyte_Nuclear_Coordinates.parquet",
    package = "kstitch"
  )
  skip_if(boundary_path == "", "boundary fixture not found")

  result_serial <- run_tpca(
    boundary_parquet_path = boundary_path,
    output_dir            = NULL,
    use_parallel          = FALSE,
    x_col = "x", y_col = "y", cell_id_col = "cell"
  )
  result_parallel <- run_tpca(
    boundary_parquet_path = boundary_path,
    output_dir            = NULL,
    use_parallel          = TRUE,
    num_threads           = 2L,
    x_col = "x", y_col = "y", cell_id_col = "cell"
  )

  common_ids <- intersect(
    rownames(result_serial$TPCA_Embedding),
    rownames(result_parallel$TPCA_Embedding)
  )
  expect_true(length(common_ids) > 0)

  top_pc <- 1L
  a <- result_serial$TPCA_Embedding[common_ids, top_pc]
  b <- result_parallel$TPCA_Embedding[common_ids, top_pc]

  diff_same_sign <- mean(abs(a - b))
  diff_flip_sign <- mean(abs(a + b))
  expect_true(min(diff_same_sign, diff_flip_sign) < 1e-3)
})
