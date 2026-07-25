make_fake_results <- function(cell_names_a, cell_names_b, k = 3) {
  make_scores <- function(cells, prefix) {
    m <- matrix(rnorm(length(cells) * k), length(cells), k,
                dimnames = list(cells, paste0(prefix, seq_len(k))))
    m
  }
  list(
    TypeA = list(
      CC_Corr_Coefs         = runif(k, 0.5, 0.9),
      CSP_Scores            = make_scores(cell_names_a, "CSP"),
      CEP_Scores            = make_scores(cell_names_a, "CEP"),
      CSP_Vectors           = matrix(rnorm(8 * k), 8, k),
      CEP_Vectors           = matrix(rnorm(6 * k), 6, k),
      CSP_Self_Correlations = matrix(rnorm(8 * k), 8, k),
      CEP_Self_Correlations = matrix(rnorm(6 * k), 6, k),
      Anchor_Features       = setNames(paste0("ShapePC", seq_len(k)),
                                       paste0("CSP", seq_len(k))),
      Misc_CCA              = list()
    ),
    TypeB = list(
      CC_Corr_Coefs         = runif(k, 0.5, 0.9),
      CSP_Scores            = make_scores(cell_names_b, "CSP"),
      CEP_Scores            = make_scores(cell_names_b, "CEP"),
      CSP_Vectors           = matrix(rnorm(8 * k), 8, k),
      CEP_Vectors           = matrix(rnorm(6 * k), 6, k),
      CSP_Self_Correlations = matrix(rnorm(8 * k), 8, k),
      CEP_Self_Correlations = matrix(rnorm(6 * k), 6, k),
      Anchor_Features       = setNames(paste0("ShapePC", seq_len(k)),
                                       paste0("CSP", seq_len(k))),
      Misc_CCA              = list()
    )
  )
}

test_that("store_kstitch_results writes results under reduction_name", {
  set.seed(1)
  cells_a <- paste0("cell_", 1:80)
  cells_b <- paste0("cell_", 81:150)
  obj     <- make_seurat_stub(c(cells_a, cells_b))
  results <- make_fake_results(cells_a, cells_b)

  obj2 <- store_kstitch_results(obj, results, reduction_name = "cell_unregressed")

  expect_true(!is.null(obj2@misc$kstitch))
  expect_named(obj2@misc$kstitch, "cell_unregressed")
  expect_named(obj2@misc$kstitch[["cell_unregressed"]],
               c("TypeA", "TypeB"), ignore.order = TRUE)
  expect_true("Anchor_Features" %in%
                names(obj2@misc$kstitch[["cell_unregressed"]][["TypeA"]]))
})

test_that("store_kstitch_results creates correctly named DimReducs", {
  set.seed(2)
  cells_a <- paste0("cell_", 1:80)
  cells_b <- paste0("cell_", 81:150)
  obj     <- make_seurat_stub(c(cells_a, cells_b))
  results <- make_fake_results(cells_a, cells_b)

  obj2 <- store_kstitch_results(obj, results, reduction_name = "cell_unregressed")

  expect_true("cell_unregressed_csp" %in% names(obj2@reductions))
  expect_true("cell_unregressed_cep" %in% names(obj2@reductions))
})

test_that("store_kstitch_results supports multiple reductions", {
  set.seed(3)
  cells_a <- paste0("cell_", 1:80)
  cells_b <- paste0("cell_", 81:150)
  obj     <- make_seurat_stub(c(cells_a, cells_b))
  results <- make_fake_results(cells_a, cells_b)

  obj2 <- store_kstitch_results(obj, results, reduction_name = "cell_unregressed",
                                reduction_key_csp = "UCSP_", reduction_key_cep = "UCEP_")
  obj2 <- store_kstitch_results(obj2, results, reduction_name = "cell_regressed",
                                reduction_key_csp = "RCSP_", reduction_key_cep = "RCEP_")

  expect_named(obj2@misc$kstitch, c("cell_unregressed", "cell_regressed"),
               ignore.order = TRUE)
  expect_true("cell_unregressed_csp" %in% names(obj2@reductions))
  expect_true("cell_regressed_csp"   %in% names(obj2@reductions))
})

test_that("get_kstitch_results retrieves the correct group", {
  set.seed(4)
  cells_a <- paste0("cell_", 1:60)
  cells_b <- paste0("cell_", 61:120)
  obj     <- make_seurat_stub(c(cells_a, cells_b))
  results <- make_fake_results(cells_a, cells_b)

  obj2 <- store_kstitch_results(obj, results, reduction_name = "cell_unregressed")

  grp <- get_kstitch_results(obj2, reduction_name = "cell_unregressed", group = "TypeA")
  expect_true(is.list(grp))
  expect_true("CC_Corr_Coefs" %in% names(grp))

  full <- get_kstitch_results(obj2, reduction_name = "cell_unregressed")
  expect_named(full, c("TypeA", "TypeB"), ignore.order = TRUE)
})

test_that("get_kstitch_results errors informatively on missing reduction or group", {
  obj <- make_seurat_stub(paste0("cell_", 1:10))

  expect_error(get_kstitch_results(obj), "No kstitch results found")

  obj@misc$kstitch <- list(cell_unregressed = list(TypeA = list()))
  expect_error(get_kstitch_results(obj, reduction_name = "nonexistent"),
               "not found in obj@misc\\$kstitch")
  expect_error(get_kstitch_results(obj, reduction_name = "cell_unregressed",
                                   group = "TypeZ"),
               "not found in obj@misc\\$kstitch")
})
