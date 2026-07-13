# Full integration tests for run_tpca() require Python + real parquet data
# and are left for manual testing. Tests here cover:
#   - .load_pga_output error handling (missing files)
#   - store_tpca_results storage logic (mocked Seurat calls)
#   - get_tpca_results accessor and error messages

make_tpca_stub_obj <- function(cell_names) {
  list(meta.data = data.frame(row.names = cell_names))
}

make_fake_tpca_result <- function(cell_names, k = 6,
                                  contour_type = "cell") {
  emb <- matrix(rnorm(length(cell_names) * k), length(cell_names), k,
                dimnames = list(cell_names, paste0("Shape_PC", seq_len(k))))
  list(
    PGA_Embedding = emb,
    Info = list(
      variances    = runif(k),
      v_matrix     = matrix(rnorm(100 * k), 100, k),
      frechet_mean = matrix(rnorm(2 * 50), 2, 50),
      pre_shape_embedding = NULL
    ),
    contour_type = contour_type,
    output_dir   = tempdir()
  )
}

# ---- .load_pga_output -------------------------------------------------------

test_that(".load_pga_output errors clearly when PGA_Info.h5 is absent", {
  tmp <- file.path(tempdir(), "empty_pga_dir")
  dir.create(tmp, showWarnings = FALSE)
  expect_error(
    kstitch:::.load_pga_output(tmp),
    "PGA_Info.h5 not found"
  )
})

test_that(".load_pga_output errors clearly when Shape_Metadata.csv.gz is absent", {
  tmp <- file.path(tempdir(), "missing_meta_dir")
  dir.create(tmp, showWarnings = FALSE)
  # create a dummy PGA_Info.h5 so the first check passes
  dummy <- file.path(tmp, "PGA_Info.h5")
  file.create(dummy)
  expect_error(
    kstitch:::.load_pga_output(tmp),
    "Shape_Metadata.csv.gz not found"
  )
  file.remove(dummy)
})

# ---- store_tpca_results -----------------------------------------------------

test_that("store_tpca_results writes DimReduc and misc for cell contours", {

  set.seed(1)
  cells  <- paste0("cell_", 1:80)
  obj    <- make_tpca_stub_obj(cells)
  result <- make_fake_tpca_result(cells, contour_type = "cell")

  with_mocked_bindings(
    CreateDimReducObject = function(embeddings, key, assay)
      list(embeddings = embeddings, key = key),
    DefaultAssay = function(obj) "RNA",
    .package = "Seurat",
    {
      obj2 <- store_tpca_results(obj, result)
    }
  )

  expect_true(!is.null(obj2[["tpca_cell"]]))
  expect_equal(nrow(obj2[["tpca_cell"]]$embeddings), length(cells))

  tpca_misc <- obj2[["misc"]][["kstitch"]][["tpca"]][["cell"]]
  expect_true(!is.null(tpca_misc))
  expect_true(all(c("variances", "v_matrix", "frechet_mean") %in% names(tpca_misc)))
})

test_that("store_tpca_results warns when obj cells are missing from embedding", {

  set.seed(2)
  emb_cells <- paste0("cell_", 1:60)
  obj_cells <- paste0("cell_", 1:80)   # 20 extra cells not in embedding

  obj    <- make_tpca_stub_obj(obj_cells)
  result <- make_fake_tpca_result(emb_cells, contour_type = "nucleus")

  with_mocked_bindings(
    CreateDimReducObject = function(embeddings, key, assay)
      list(embeddings = embeddings),
    DefaultAssay = function(obj) "RNA",
    .package = "Seurat",
    {
      expect_warning(
        obj2 <- store_tpca_results(obj, result),
        "no TPCA embedding"
      )
    }
  )

  expect_equal(nrow(obj2[["tpca_nucleus"]]$embeddings), length(obj_cells))
})

# ---- get_tpca_results -------------------------------------------------------

test_that("get_tpca_results retrieves the correct contour type", {

  obj <- make_tpca_stub_obj(paste0("cell_", 1:10))
  obj[["misc"]] <- list(
    kstitch = list(
      tpca = list(
        cell     = list(variances = 1:3),
        nucleus  = list(variances = 4:6)
      )
    )
  )

  expect_equal(get_tpca_results(obj, "cell")$variances,     1:3)
  expect_equal(get_tpca_results(obj, "nucleus")$variances,  4:6)

  full <- get_tpca_results(obj)
  expect_named(full, c("cell", "nucleus"))
})

test_that("get_tpca_results errors informatively on missing results", {

  obj <- make_tpca_stub_obj(paste0("cell_", 1:5))
  obj[["misc"]] <- list()

  expect_error(get_tpca_results(obj), "No TPCA results found")

  obj[["misc"]][["kstitch"]] <- list(tpca = list(cell = list()))
  expect_error(get_tpca_results(obj, "nucleus"),
               "contour_type 'nucleus' not found")
})
