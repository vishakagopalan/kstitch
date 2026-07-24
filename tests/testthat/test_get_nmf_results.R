make_obj_with_nmf <- function(groups = c("TypeA", "TypeB")) {
  n <- 100
  cells <- paste0("cell_", seq_len(n * length(groups)))
  obj   <- make_seurat_stub(cells)
  for (grp in groups) {
    emb <- matrix(rnorm(n * 4), n, 4,
                  dimnames = list(paste0("cell_", seq_len(n)), paste0("NMF", 1:4)))
    obj[[paste0("nmf_", grp)]] <- Seurat::CreateDimReducObject(
      embeddings = emb,
      key        = paste0("NMF", gsub("[^A-Za-z0-9]", "", grp), "_"),
      assay      = Seurat::DefaultAssay(obj)
    )
    obj[[paste0("nmf_", grp)]]@misc <- list(
      Fit_Error        = 0.1,
      Fit_Summary      = data.frame(),
      Factor_Gene_List = list()
    )
  }
  obj
}

test_that("get_nmf_results returns misc for a specified group", {
  obj <- make_obj_with_nmf()
  res <- get_nmf_results(obj, group = "TypeA")
  expect_true(is.list(res))
  expect_true("Fit_Error" %in% names(res))
})

test_that("get_nmf_results with group = NULL returns all groups as named list", {
  obj <- make_obj_with_nmf(groups = c("TypeA", "TypeB"))
  res <- get_nmf_results(obj, group = NULL)
  expect_true(is.list(res))
  expect_named(res, c("TypeA", "TypeB"), ignore.order = TRUE)
  expect_true("Fit_Error" %in% names(res$TypeA))
})

test_that("get_nmf_results errors when no NMF reductions exist", {
  obj <- make_seurat_stub(paste0("cell_", 1:20))
  expect_error(get_nmf_results(obj, group = NULL), "No reductions with prefix")
})

test_that("get_nmf_results errors when specified group is absent", {
  obj <- make_obj_with_nmf(groups = "TypeA")
  expect_error(get_nmf_results(obj, group = "TypeB"), "not found in obj")
})
