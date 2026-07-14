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

test_that("store_kstitch_results writes misc entries for each group", {

  set.seed(1)
  cells_a <- paste0("cell_", 1:80)
  cells_b <- paste0("cell_", 81:150)
  obj     <- make_seurat_stub(c(cells_a, cells_b))
  results <- make_fake_results(cells_a, cells_b)

  obj2 <- store_kstitch_results(obj, results)

  expect_true(!is.null(obj2@misc$kstitch))
  expect_named(obj2@misc$kstitch, c("TypeA", "TypeB"), ignore.order = TRUE)
  expect_false("CSP_Scores" %in% names(obj2@misc$kstitch[["TypeA"]]))
  expect_false("CEP_Scores" %in% names(obj2@misc$kstitch[["TypeA"]]))
  expect_true("Anchor_Features" %in% names(obj2@misc$kstitch[["TypeA"]]))
})

test_that("get_kstitch_results retrieves the correct group", {

  set.seed(2)
  cells_a <- paste0("cell_", 1:60)
  cells_b <- paste0("cell_", 61:120)
  obj     <- make_seurat_stub(c(cells_a, cells_b))
  results <- make_fake_results(cells_a, cells_b)

  obj2 <- store_kstitch_results(obj, results)

  grp <- get_kstitch_results(obj2, "TypeA")
  expect_true(is.list(grp))
  expect_true("CC_Corr_Coefs" %in% names(grp))

  full <- get_kstitch_results(obj2)
  expect_named(full, c("TypeA", "TypeB"), ignore.order = TRUE)
})

test_that("get_kstitch_results errors informatively on missing group", {

  obj <- make_seurat_stub(paste0("cell_", 1:10))

  expect_error(get_kstitch_results(obj), "No kstitch results found")

  obj@misc$kstitch <- list(TypeA = list())
  expect_error(get_kstitch_results(obj, "TypeZ"),
               "not found in obj@misc\\$kstitch")
})
