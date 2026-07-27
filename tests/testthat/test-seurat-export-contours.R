# Tests for export_seurat_contours() and its internal helpers.
# Python / TPCA execution is not tested here; tests cover everything up to
# the reticulate boundary.

# ── Fixtures ──────────────────────────────────────────────────────────────────

make_seg_list <- function(cell_ids, contour_type = "cell", use_cell_col = FALSE) {
  n   <- length(cell_ids)
  # Each cell gets a small square contour (4 vertices)
  id_col <- if (use_cell_col) "cell" else "cell_id"
  df <- data.frame(
    x       = rep(c(0, 1, 1, 0), n),
    y       = rep(c(0, 0, 1, 1), n),
    cell_id = rep(cell_ids, each = 4)
  )
  names(df)[names(df) == "cell_id"] <- id_col

  slot_name <- if (contour_type == "nucleus") "nucleus_segmentations" else "segmentations"
  stats::setNames(list(df), slot_name)
}

CELL_IDS <- paste0("cell_", 1:5)

# ── .extract_contour_df ───────────────────────────────────────────────────────

test_that(".extract_contour_df returns x, y, cell_id for cell contours", {
  seg <- make_seg_list(CELL_IDS, "cell")
  df  <- kstitch:::.extract_contour_df(seg, "cell")
  expect_named(df, c("x", "y", "cell_id"), ignore.order = TRUE)
  expect_equal(sort(unique(df$cell_id)), sort(CELL_IDS))
})

test_that(".extract_contour_df accepts 'cell' column alias", {
  seg <- make_seg_list(CELL_IDS, "cell", use_cell_col = TRUE)
  df  <- kstitch:::.extract_contour_df(seg, "cell")
  expect_true("cell_id" %in% names(df))
})

test_that(".extract_contour_df returns x, y, cell_id for nucleus contours", {
  seg <- make_seg_list(CELL_IDS, "nucleus")
  df  <- kstitch:::.extract_contour_df(seg, "nucleus")
  expect_named(df, c("x", "y", "cell_id"), ignore.order = TRUE)
})

test_that(".extract_contour_df errors when slot is missing", {
  seg <- make_seg_list(CELL_IDS, "cell")   # has 'segmentations', not 'nucleus_segmentations'
  expect_error(
    kstitch:::.extract_contour_df(seg, "nucleus"),
    "nucleus_segmentations not found"
  )
})

test_that(".extract_contour_df errors when cell slot missing for cell type", {
  expect_error(
    kstitch:::.extract_contour_df(list(), "cell"),
    "'segmentations' not found"
  )
})

# ── .normalise_columns ────────────────────────────────────────────────────────

test_that(".normalise_columns renames candidate column to target", {
  df  <- data.frame(cell = 1:3, x = 1:3, y = 1:3)
  out <- kstitch:::.normalise_columns(df, list(cell_id = c("cell", "cell_id")))
  expect_true("cell_id" %in% names(out))
  expect_false("cell" %in% names(out))
})

test_that(".normalise_columns errors when no candidate found", {
  df <- data.frame(a = 1:3, x = 1:3, y = 1:3)
  expect_error(
    kstitch:::.normalise_columns(df, list(cell_id = c("cell", "cell_id"))),
    "Could not find column"
  )
})

# ── export_seurat_contours() input validation ─────────────────────────────────

test_that("contour_type is matched and rejects bad values", {
  obj <- make_seurat_stub(CELL_IDS)
  seg <- make_seg_list(CELL_IDS, "cell")
  expect_error(
    export_seurat_contours(obj, seg, contour_type = "mitochondria"),
    regexp = "arg"   # match.arg error
  )
})

test_that("errors when no contours remain after cell_id subsetting", {
  obj <- make_seurat_stub(CELL_IDS)
  seg <- make_seg_list(CELL_IDS, "cell")
  expect_error(
    export_seurat_contours(obj, seg, cell_ids = "nonexistent_cell"),
    "No contours remain"
  )
})

test_that("cell_ids defaults to all cells in obj", {
  # Verify subsetting logic: a seg_list with extra cells is trimmed to obj cells.
  # We can't run Python, so mock at the reticulate boundary.
  skip_if_not_installed("mockery")
  obj      <- make_seurat_stub(CELL_IDS[1:3])
  all_ids  <- CELL_IDS          # seg_list has 5 cells
  seg      <- make_seg_list(all_ids, "cell")

  # Capture the cell_ids actually passed into compute_pre_shape_embedding
  captured <- NULL
  mock_mod <- list(
    compute_pre_shape_embedding = function(...) { captured <<- list(...)$cell_ids_to_analyze },
    run_kendall_tpca            = function(...) list(NULL, NULL, numeric(0))
  )
  mockery::stub(export_seurat_contours, "reticulate::import_from_path", mock_mod)
  mockery::stub(export_seurat_contours, ".load_kendall_tpca_output",
                list(TPCA_Embedding = NULL, Info = NULL, Metadata = NULL))

  suppressMessages(export_seurat_contours(obj, seg))
  expect_setequal(as.character(reticulate::py_to_r(captured)), CELL_IDS[1:3])
})

# ── return_results = FALSE ────────────────────────────────────────────────────

test_that("return_results = FALSE returns list with output_path", {
  skip_if_not_installed("mockery")
  obj <- make_seurat_stub(CELL_IDS)
  seg <- make_seg_list(CELL_IDS, "cell")

  mock_mod <- list(
    compute_pre_shape_embedding = function(...) invisible(NULL),
    run_kendall_tpca            = function(...) list(NULL, NULL, numeric(0))
  )
  mockery::stub(export_seurat_contours, "reticulate::import_from_path", mock_mod)

  result <- suppressMessages(
    export_seurat_contours(obj, seg, return_results = FALSE)
  )
  expect_named(result, "output_path")
  expect_type(result$output_path, "character")
})
