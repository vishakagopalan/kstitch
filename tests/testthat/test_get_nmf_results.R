make_obj_with_nmf <- function(suffixes = c("Keratinocytes", "Fibroblasts")) {
  n     <- 100
  cells <- paste0("cell_", seq_len(n * length(suffixes)))
  obj   <- make_seurat_stub(cells)

  for (sfx in suffixes) {
    emb <- matrix(rnorm(n * 4), n, 4,
                  dimnames = list(paste0("cell_", seq_len(n)),
                                  paste0("NMF", 1:4)))
    safe <- gsub("[^A-Za-z0-9]", "", sfx)
    obj[[paste0("nmf_", sfx)]] <- Seurat::CreateDimReducObject(
      embeddings = emb,
      key        = paste0("NMF", safe, "_"),
      assay      = Seurat::DefaultAssay(obj)
    )
    obj[[paste0("nmf_", sfx)]]@misc <- list(
      Fit_Error        = 0.1,
      Fit_Summary      = data.frame(),
      Factor_Gene_List = list()
    )
  }
  obj
}

test_that("get_nmf_results returns misc for a specified reduction_suffix", {
  obj <- make_obj_with_nmf()
  res <- get_nmf_results(obj, reduction_suffix = "Keratinocytes")
  expect_true(is.list(res))
  expect_true("Fit_Error" %in% names(res))
})

test_that("get_nmf_results with reduction_suffix = NULL returns all as named list", {
  obj <- make_obj_with_nmf(suffixes = c("Keratinocytes", "Fibroblasts"))
  res <- get_nmf_results(obj, reduction_suffix = NULL)
  expect_true(is.list(res))
  expect_named(res, c("Keratinocytes", "Fibroblasts"), ignore.order = TRUE)
  expect_true("Fit_Error" %in% names(res$Keratinocytes))
})

test_that("get_nmf_results errors with available reductions when suffix is absent", {
  obj <- make_obj_with_nmf(suffixes = "Keratinocytes")
  expect_error(
    get_nmf_results(obj, reduction_suffix = "Fibroblasts"),
    "nmf_Keratinocytes"   # available reductions shown in error
  )
})

test_that("get_nmf_results errors with all available reductions when none match prefix", {
  obj <- make_seurat_stub(paste0("cell_", 1:20))
  expect_error(
    get_nmf_results(obj, reduction_suffix = NULL),
    "No reductions with prefix"
  )
})

test_that("get_nmf_results single-group 'all' suffix works (store_nmf_results default)", {
  n   <- 50
  obj <- make_seurat_stub(paste0("cell_", 1:n))
  emb <- matrix(rnorm(n * 4), n, 4,
                dimnames = list(paste0("cell_", 1:n), paste0("NMF", 1:4)))
  obj[["nmf_all"]] <- Seurat::CreateDimReducObject(
    embeddings = emb, key = "NMFall_", assay = Seurat::DefaultAssay(obj)
  )
  obj[["nmf_all"]]@misc <- list(Fit_Error = 0.2, Fit_Summary = data.frame(),
                                Factor_Gene_List = list())

  res <- get_nmf_results(obj, reduction_suffix = "all")
  expect_equal(res$Fit_Error, 0.2)
})
