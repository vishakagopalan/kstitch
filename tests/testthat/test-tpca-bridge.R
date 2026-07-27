make_fake_tpca_result <- function(cell_names, k = 6, contour_type = "cell") {
  emb <- matrix(rnorm(length(cell_names) * k), length(cell_names), k,
                dimnames = list(cell_names, paste0("Shape_PC", seq_len(k))))
  list(
    TPCA_Embedding = emb,
    Info = list(
      variances           = sort(runif(k, 0.5, 5), decreasing = TRUE),
      v_matrix            = matrix(rnorm(100 * k), 100, k),
      frechet_mean        = matrix(rnorm(2 * 50), 2, 50),
      pre_shape_embedding = NULL
    ),
    contour_type = contour_type,
    output_dir   = tempdir()
  )
}


# ---- .load_kendall_tpca_output ----------------------------------------------

test_that(".load_kendall_tpca_output errors clearly when TPCA_Info.h5 is absent", {
  tmp <- file.path(tempdir(), "empty_tpca_dir")
  dir.create(tmp, showWarnings = FALSE)
  expect_error(kstitch:::.load_kendall_tpca_output(tmp), "TPCA_Info.h5 not found")
})

test_that(".load_kendall_tpca_output errors clearly when Shape_Metadata.csv.gz is absent", {
  tmp <- file.path(tempdir(), "missing_meta_dir")
  dir.create(tmp, showWarnings = FALSE)
  dummy <- file.path(tmp, "TPCA_Info.h5")
  file.create(dummy)
  expect_error(kstitch:::.load_kendall_tpca_output(tmp), "Shape_Metadata.csv.gz not found")
  file.remove(dummy)
})


# ---- store_tpca_results — single group --------------------------------------

test_that("store_tpca_results writes DimReduc and misc for single-group cell result", {
  set.seed(1)
  cells  <- paste0("cell_", 1:80)
  obj    <- make_seurat_stub(cells)
  result <- make_fake_tpca_result(cells, contour_type = "cell")

  obj2 <- store_tpca_results(obj, result)

  expect_true(!is.null(obj2[["tpca_cell"]]))
  expect_equal(nrow(Seurat::Embeddings(obj2, "tpca_cell")), length(cells))
  expect_true(!is.null(obj2@misc$kstitch$tpca$cell))
  expect_true(all(c("variances", "v_matrix", "frechet_mean") %in%
                    names(obj2@misc$kstitch$tpca$cell)))
})

test_that("store_tpca_results warns when obj cells are missing from single-group embedding", {
  set.seed(2)
  emb_cells <- paste0("cell_", 1:60)
  obj_cells <- paste0("cell_", 1:80)
  obj       <- make_seurat_stub(obj_cells)
  result    <- make_fake_tpca_result(emb_cells, contour_type = "nucleus")

  expect_warning(
    obj2 <- store_tpca_results(obj, result),
    "no TPCA embedding"
  )
  expect_equal(nrow(Seurat::Embeddings(obj2, "tpca_nucleus")), length(obj_cells))
})


# ---- get_tpca_results -------------------------------------------------------

test_that("get_tpca_results retrieves the correct contour type", {
  obj <- make_seurat_stub(paste0("cell_", 1:10))
  obj@misc$kstitch <- list(
    tpca = list(
      cell    = list(variances = 1:3),
      nucleus = list(variances = 4:6)
    )
  )

  expect_equal(get_tpca_results(obj, "cell")$variances,    1:3)
  expect_equal(get_tpca_results(obj, "nucleus")$variances, 4:6)
  expect_named(get_tpca_results(obj), c("cell", "nucleus"))
})

test_that("get_tpca_results errors informatively on missing results", {
  obj <- make_seurat_stub(paste0("cell_", 1:5))
  expect_error(get_tpca_results(obj), "No TPCA results found")

  obj@misc$kstitch <- list(tpca = list(cell = list()))
  expect_error(get_tpca_results(obj, "nucleus"),
               "contour_type 'nucleus' not found")
})

test_that("get_tpca_results shows available groups when group is missing", {
  obj <- make_seurat_stub(paste0("cell_", 1:5))
  obj@misc$kstitch <- list(tpca = list(cell = list(TypeA = list(), TypeB = list())))
  expect_error(
    get_tpca_results(obj, "cell", group = "TypeC"),
    "TypeA"
  )
})

# ---- load_kstitch_results type=tpca -----------------------------------------

test_that("load_kstitch_results type=tpca returns a flat result", {
  result <- load_kstitch_results(
    path = system.file(
      "extdata", "xenium_test", "tpca_fixtures", "Unit_Test_Keratinocytes_Nuclei",
      package = "kstitch"
    ),
    type = "tpca"
  )
  # flat result — TPCA_Embedding at top level, no "all" wrapper
  expect_true("TPCA_Embedding" %in% names(result))
  expect_true("Info"           %in% names(result))
  expect_true("Metadata"       %in% names(result))
  expect_false("all" %in% names(result))
})
