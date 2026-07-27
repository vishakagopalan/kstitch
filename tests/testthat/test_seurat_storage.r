make_fake_result <- function(cell_names, k = 3) {
  make_scores <- function(cells, prefix) {
    matrix(rnorm(length(cells) * k), length(cells), k,
           dimnames = list(cells, paste0(prefix, seq_len(k))))
  }
  list(
    CC_Corr_Coefs       = runif(k, 0.5, 0.9),
    CSP_Scores          = make_scores(cell_names, "CSP"),
    CEP_Scores          = make_scores(cell_names, "CEP"),
    CSP_Vectors         = matrix(rnorm(8 * k), 8, k),
    CEP_Vectors         = matrix(rnorm(6 * k), 6, k),
    Shape_Corr_With_CSP = matrix(rnorm(8 * k), 8, k),
    Exp_Corr_With_CEP   = matrix(rnorm(6 * k), 6, k),
    Anchor_Features     = setNames(paste0("ShapePC", seq_len(k)),
                                   paste0("CSP", seq_len(k))),
    Misc_CCA            = list()
  )
}

# ---------------------------------------------------------------------------
# store_kstitch_results
# ---------------------------------------------------------------------------

test_that("store_kstitch_results writes result under reduction_name in misc", {
  set.seed(1)
  cells  <- paste0("cell_", 1:80)
  obj    <- make_seurat_stub(cells)
  result <- make_fake_result(cells)

  obj2 <- store_kstitch_results(obj, result, reduction_name = "nuc_regressed")

  expect_true(!is.null(obj2@misc$kstitch))
  expect_named(obj2@misc$kstitch, "nuc_regressed")
  expect_true("Anchor_Features" %in% names(obj2@misc$kstitch[["nuc_regressed"]]))
})

test_that("store_kstitch_results creates correctly named DimReducs", {
  set.seed(2)
  cells  <- paste0("cell_", 1:80)
  obj    <- make_seurat_stub(cells)
  result <- make_fake_result(cells)

  obj2 <- store_kstitch_results(obj, result, reduction_name = "cell_regressed",
                                reduction_key_csp = "RCSP_",
                                reduction_key_cep = "RCEP_")

  expect_true("cell_regressed_csp" %in% names(obj2@reductions))
  expect_true("cell_regressed_cep" %in% names(obj2@reductions))
})

test_that("store_kstitch_results supports multiple independent reductions", {
  set.seed(3)
  cells   <- paste0("cell_", 1:80)
  obj     <- make_seurat_stub(cells)
  result1 <- make_fake_result(cells)
  result2 <- make_fake_result(cells)

  obj2 <- store_kstitch_results(obj,  result1, reduction_name = "cell_unregressed",
                                reduction_key_csp = "UCSP_", reduction_key_cep = "UCEP_")
  obj2 <- store_kstitch_results(obj2, result2, reduction_name = "cell_regressed",
                                reduction_key_csp = "RCSP_", reduction_key_cep = "RCEP_")

  expect_named(obj2@misc$kstitch, c("cell_unregressed", "cell_regressed"),
               ignore.order = TRUE)
  expect_true("cell_unregressed_csp" %in% names(obj2@reductions))
  expect_true("cell_regressed_csp"   %in% names(obj2@reductions))
})

test_that("store_kstitch_results fills NA for obj cells absent from result", {
  set.seed(4)
  scored_cells <- paste0("cell_", 1:60)
  all_cells    <- paste0("cell_", 1:80)   # 20 extra not in result
  obj          <- make_seurat_stub(all_cells)
  result       <- make_fake_result(scored_cells)

  expect_warning(
    obj2 <- store_kstitch_results(obj, result, reduction_name = "test"),
    "NA scores"
  )

  csp_emb <- Seurat::Embeddings(obj2, "test_csp")
  expect_equal(nrow(csp_emb), length(all_cells))
  expect_true(all(is.na(csp_emb[paste0("cell_", 61:80), ])))
  expect_true(all(!is.na(csp_emb[scored_cells, ])))
})

# ---------------------------------------------------------------------------
# get_kstitch_results
# ---------------------------------------------------------------------------

test_that("get_kstitch_results retrieves the stored result", {
  set.seed(5)
  cells  <- paste0("cell_", 1:60)
  obj    <- make_seurat_stub(cells)
  result <- make_fake_result(cells)

  obj2 <- store_kstitch_results(obj, result, reduction_name = "nuc_regressed")
  out  <- get_kstitch_results(obj2, reduction_name = "nuc_regressed")

  expect_true(is.list(out))
  expect_true("CC_Corr_Coefs" %in% names(out))
})

test_that("get_kstitch_results errors with available reductions when name is missing", {
  set.seed(6)
  cells  <- paste0("cell_", 1:60)
  obj    <- make_seurat_stub(cells)
  result <- make_fake_result(cells)
  obj2   <- store_kstitch_results(obj, result, reduction_name = "nuc_regressed")

  expect_error(
    get_kstitch_results(obj2, reduction_name = "nonexistent"),
    "nuc_regressed"   # available reductions shown in error
  )
})

test_that("get_kstitch_results errors when no kstitch results exist", {
  obj <- make_seurat_stub(paste0("cell_", 1:10))
  expect_error(get_kstitch_results(obj), "No kstitch results found")
})
