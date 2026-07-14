#' Export contours from a Seurat object and run TPCA
#'
#' Accepts the raw segmentation list returned by \code{ReadXenium} or
#' \code{ReadNanostring}, normalises the contour data frame, passes it
#' directly to PGA (no intermediate parquet), runs TPCA, and returns a
#' result ready for \code{store_tpca_results()}.
#'
#' @param obj          A Seurat object. Used to optionally subset contours to
#'                     cells present in the object.
#' @param seg_list     The list returned by \code{ReadXenium} or
#'                     \code{ReadNanostring}.
#' @param contour_type One of \code{"cell"} or \code{"nucleus"}.
#' @param output_dir   Directory for TPCA intermediate outputs.
#' @param cell_ids     Optional character vector of cell IDs to analyse.
#'                     Defaults to all cells in \code{obj}.
#' @param num_vertices Integer. Vertices to sample per contour (default 50).
#' @param eta          Numeric. TPCA eta parameter (default 1).
#' @param use_parallel Logical. Use parallel processing (default FALSE).
#' @param num_threads  Integer. Threads when \code{use_parallel = TRUE}.
#' @param frechet_mean_tol Numeric. Convergence tolerance for Frechet mean.
#' @param max_frechet_iter Integer. Maximum Frechet mean iterations.
#'
#' @return A list in the format expected by \code{store_tpca_results()}.
#' @export
export_seurat_contours <- function(obj,
                                   seg_list,
                                   contour_type      = c("cell", "nucleus"),
                                   output_dir,
                                   cell_ids          = NULL,
                                   num_vertices      = 50L,
                                   eta               = 1,
                                   use_parallel      = FALSE,
                                   num_threads       = 8L,
                                   frechet_mean_tol  = 1e-4,
                                   max_frechet_iter  = 1000L) {

  contour_type <- match.arg(contour_type)

  # ── 1. Extract and normalise the contour data frame ─────────────────────────
  df <- .extract_contour_df(seg_list, contour_type)

  # ── 2. Optionally subset to cells present in obj ────────────────────────────
  if (is.null(cell_ids)) {
    cell_ids <- rownames(obj@meta.data)
  }
  df <- df[df$cell_id %in% cell_ids, , drop = FALSE]
  if (nrow(df) == 0L)
    stop("No contours remain after subsetting to cells in obj.")

  # ── 3. Load PGA ──────────────────────────────────────────────────────────────
  pga_py_dir <- system.file("python", package = "kstitch")
  if (!nzchar(pga_py_dir))
    stop("Could not locate inst/python/ inside the kstitch package.")
  kendall_tpca <- reticulate::import_from_path("kendall_tpca", path = pga_py_dir)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  py_cell_ids <- reticulate::r_to_py(as.list(cell_ids))

  # ── 4. Pre-shape embedding (in-memory, no parquet) ───────────────────────────
  message("Computing pre-shape embedding ...")
  kendall_tpca$compute_pre_shape_embedding(
    pre_shape_output_dir   = output_dir,
    df                     = reticulate::r_to_py(df),
    num_vertices_to_sample = as.integer(num_vertices),
    cell_ids_to_analyze    = py_cell_ids,
    x_vertex_col           = "x",
    y_vertex_col           = "y",
    cell_id_col            = "cell_id"
  )

  # ── 5. TPCA ──────────────────────────────────────────────────────────────────
  message("Running TPCA ...")
  kendall_tpca$run_kendall_tpca(
    pre_shape_input_dir   = output_dir,
    pga_output_dir        = output_dir,
    cell_ids_to_analyze   = py_cell_ids,
    cell_id_col           = "cell_id",
    max_frechet_mean_iter = as.integer(max_frechet_iter),
    eta                   = eta,
    use_parallel          = use_parallel,
    num_threads           = as.integer(num_threads),
    frechet_mean_tol      = frechet_mean_tol
  )

  # ── 6. Read outputs from disk ─────────────────────────────────────────────
  message("Reading TPCA results ...")
  result              <- .load_kendall_tpca_output(output_dir)
  result$contour_type <- contour_type
  result$output_dir   <- output_dir
  result
}


# ── Internal helpers ───────────────────────────────────────────────────────────

#' @keywords internal
.extract_contour_df <- function(seg_list, contour_type) {

  # Xenium: ReadXenium returns nucleus_segmentations / segmentations
  # CosMx:  ReadNanostring returns segmentations only (no nucleus)
  slot_name <- switch(contour_type,
    nucleus = "nucleus_segmentations",
    cell    = "segmentations"
  )

  if (!slot_name %in% names(seg_list)) {
    if (contour_type == "nucleus")
      stop("nucleus_segmentations not found in seg_list. ",
           "Xenium objects may lack nucleus contours depending on the run. ",
           "Try contour_type = 'cell'. ",
           "Note: CosMx nucleus contour support will be added in a future ",
           "release depending on platform availability.")
    else
      stop(sprintf("'%s' not found in seg_list.", slot_name))
  }

  df <- seg_list[[slot_name]]

  # Xenium columns: cell, x, y, <NA integer col>
  # CosMx columns:  x, y, cell
  # Normalise to: x, y, cell_id
  col_map <- list(
    cell_id = c("cell", "cell_id"),
    x       = "x",
    y       = "y"
  )

  df <- .normalise_columns(df, col_map)
  df[, c("x", "y", "cell_id")]   # drop cell_numeric_id (Xenium-only integer column)
}

#' @keywords internal
.normalise_columns <- function(df, col_map) {
  for (target in names(col_map)) {
    candidates <- col_map[[target]]
    found <- intersect(candidates, names(df))
    if (length(found) == 0L)
      stop(sprintf(
        "Could not find column '%s' in segmentation data frame. ",
        "Available columns: %s.",
        target, paste(names(df), collapse = ", ")
      ))
    if (!target %in% names(df))
      names(df)[names(df) == found[1L]] <- target
  }
  df
}
