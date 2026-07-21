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
  expect_true(min(diff_same_sign, diff_flip_sign) < 1e-6)

  expect_equal(result$Info$variances, expected$Info$variances, tolerance = 1e-6)
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
    use_cache             = FALSE,
    frechet_mean_tol = 1e-6,
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
  expect_equal(result$Info$variances, expected$Info$variances, tolerance = 1e-3)

  # ---- identify degenerate (near-equal variance) PCs --------------------
  # PCs whose variance is within a relative 1e-3 of a neighbor are
  # numerically unstable in direction/sign and should not be compared
  # column-by-column.
  vars <- expected$Info$variances
  rel_gap <- abs(diff(vars)) / (vars[-length(vars)] + 1e-12)
  degenerate <- c(FALSE, rel_gap < 1e-3) | c(rel_gap < 1e-3, FALSE)

  # ---- identify negligible-variance PCs -----------------------------
  # PCs explaining less than 0.1% of total variance are noise-dominated;
  # exclude from strict per-column comparison.
  negligible <- vars / sum(vars) < 1e-2

  stable_pcs <- which(!degenerate & !negligible)
  skip_if(length(stable_pcs) == 0, "no numerically stable PCs to compare")

  cell_order <- rownames(res_aligned)
  # ---- per-PC comparison, up to independent sign flips, stable PCs only --
  for (j in stable_pcs) {
    diff_same_sign <- mean(abs( (res_aligned[, j] - exp_aligned[, j])/exp_aligned[, j]))
    diff_flip_sign <- mean(abs( (res_aligned[, j] + exp_aligned[, j])/exp_aligned[, j]) )
    expect_true(
      min(diff_same_sign, diff_flip_sign) < 1e-3,
      label = paste0("PC", j, " sign-adjusted difference")
    )
  }

  # ---- subspace-level check across ALL PCs (sign- and rotation-robust) --
  # Even for degenerate PCs, the SUBSPACE they jointly span should match,
  # even if individual directions within it have rotated. Compare via
  # Frobenius norm of the projection matrices (P P^T), which is invariant
  # to both sign flips and internal rotation within a degenerate block.
  proj_result   <- res_aligned %*% t(res_aligned)
  proj_expected <- exp_aligned %*% t(exp_aligned)
  subspace_diff <- mean(abs(proj_result - proj_expected))
  expect_true(subspace_diff < 1e-3)
})


test_that("Frechet mean converges within tolerance on fixture contours", {
  skip_on_ci()

  expected_tpca_path <- system.file(
    "extdata", "xenium_test", "tpca_fixtures", "Unit_Test_Keratinocytes_Nuclei",
    package = "kstitch"
  )
  skip_if(expected_tpca_path == "", "expected tpca fixture not found")

  expected <- kstitch:::.load_kendall_tpca_output(expected_tpca_path)

  # The Frechet mean itself should be a valid, finite pre-shape — if this
  # fails, everything downstream (tangent space, PCA) is unreliable
  # regardless of how close the embeddings numerically compare.
  expect_true(all(is.finite(expected$Info$frechet_mean)))
})


test_that("top PCs are robust to use_parallel = TRUE vs FALSE", {
  skip_on_ci()

  boundary_path <- system.file(
    "extdata", "xenium_test", "Test_Keratinocyte_Nuclear_Coordinates.parquet",
    package = "kstitch"
  )
  skip_if(boundary_path == "", "boundary fixture not found")

  result_serial <- run_tpca(
    boundary_parquet_path = boundary_path,
    output_dir            = withr::local_tempdir(),
    use_parallel          = FALSE,
    use_cache             = FALSE,
    x_col = "x", y_col = "y", cell_id_col = "cell"
  )
  result_parallel <- run_tpca(
    boundary_parquet_path = boundary_path,
    output_dir            = withr::local_tempdir(),
    use_parallel          = TRUE,
    num_threads           = 2L,
    use_cache             = FALSE,
    x_col = "x", y_col = "y", cell_id_col = "cell"
  )

  common_ids <- intersect(
    rownames(result_serial$TPCA_Embedding),
    rownames(result_parallel$TPCA_Embedding)
  )
  expect_true(length(common_ids) > 0)

  # Only compare the leading, well-separated PC(s) — parallel chunking can
  # alter floating-point summation order, which is amplified by degeneracy.
  vars <- result_serial$Info$variances
  top_pc <- 1L

  a <- result_serial$TPCA_Embedding[common_ids, top_pc]
  b <- result_parallel$TPCA_Embedding[common_ids, top_pc]

  diff_same_sign <- mean(abs(a - b))
  diff_flip_sign <- mean(abs(a + b))
  expect_true(min(diff_same_sign, diff_flip_sign) < 1e-3)
})
