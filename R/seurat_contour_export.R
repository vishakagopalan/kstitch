#' Export contours from a Seurat object and run TPCA
#'
#' Accepts the raw segmentation list returned by \code{ReadXenium} or
#' \code{ReadNanostring}, normalises the contour data frame, passes it
#' directly to kendall_tpca (no intermediate parquet), runs TPCA, and returns a
#' result ready for \code{store_tpca_results()}.
#'
#' @param obj          A Seurat object. Used to optionally subset contours to
#'                     cells present in the object.
#' @param seg_list     The list returned by \code{ReadXenium} or
#'                     \code{ReadNanostring}.
#' @param contour_type One of \code{"cell"} or \code{"nucleus"}.
#' @param output_dir   Directory for TPCA intermediate outputs. Optional. If
#'                     \code{NULL} (the default), a temporary directory is
#'                     used. When \code{return_results = TRUE} the temp
#'                     directory is deleted after the call; when
#'                     \code{return_results = FALSE} it is kept so saved
#'                     files remain accessible. An explicit path is created
#'                     if absent and always left on disk.
#' @param cell_ids     Optional character vector of cell IDs to analyse.
#'                     Defaults to all cells in \code{obj}.
#' @param num_vertices Integer. Vertices to sample per contour (default 50).
#' @param eta          Numeric. TPCA eta parameter (default 1).
#' @param use_parallel Logical. Use parallel processing (default FALSE).
#' @param num_threads  Integer. Threads when \code{use_parallel = TRUE}.
#' @param frechet_mean_tol Numeric. Convergence tolerance for Frechet mean.
#' @param max_frechet_iter Integer. Maximum Frechet mean iterations.
#' @param return_results Logical. If \code{TRUE} (default), results are
#'   loaded into R and returned. If \code{FALSE}, results are written to
#'   disk and only the output path is returned, with a message showing
#'   how to load them via \code{\link{load_kstitch_results}}.
#'
#' @return When \code{return_results = TRUE}: a named list with element
#'   \code{"all"} containing \code{TPCA_Embedding}, \code{Info},
#'   \code{Metadata}, \code{contour_type}, \code{group}, and
#'   \code{is_groupwise}.
#'   When \code{return_results = FALSE}: a list with \code{output_path}
#'   (character scalar) giving the directory written by
#'   \code{load_kstitch_results(..., type = "tpca")}.
#' @seealso \code{\link{load_kstitch_results}}, \code{\link{store_tpca_results}}
#' @export
export_seurat_contours <- function(obj,
                                   seg_list,
                                   contour_type      = c("cell", "nucleus"),
                                   output_dir        = NULL,
                                   cell_ids          = NULL,
                                   num_vertices      = 50L,
                                   eta               = 1,
                                   use_parallel      = FALSE,
                                   num_threads       = 8L,
                                   frechet_mean_tol  = 1e-4,
                                   max_frechet_iter  = 1000L,
                                   return_results    = TRUE) {

  contour_type <- match.arg(contour_type)

  # ── 1. Extract and normalise the contour data frame ─────────────────────────
  df <- .extract_contour_df(seg_list, contour_type)

  # ── 2. Resolve cell_ids ──────────────────────────────────────────────────────
  if (is.null(cell_ids))
    cell_ids <- rownames(obj@meta.data)
  df <- df[df$cell_id %in% cell_ids, , drop = FALSE]
  if (nrow(df) == 0L)
    stop("No contours remain after subsetting to cells in obj.")

  # ── 3. Resolve output directory ──────────────────────────────────────────────
  if (!return_results && is.null(output_dir)) {
    work_dir <- file.path(tempdir(), paste0("kstitch_tpca_", .random_id()))
    dir.create(work_dir, recursive = TRUE)
    resolved <- list(path = work_dir, is_temp = FALSE)
  } else {
    resolved <- .resolve_tpca_output_dir(output_dir)
  }
  work_dir <- resolved$path

  # ── 4. Load Python module ────────────────────────────────────────────────────
  kendall_tpca_py_dir <- system.file("python", package = "kstitch")
  if (!nzchar(kendall_tpca_py_dir))
    stop("Could not locate inst/python/ inside the kstitch package.")
  kendall_tpca <- reticulate::import_from_path(
    "kendall_tpca", path = kendall_tpca_py_dir
  )

  # ── 5. Run TPCA ──────────────────────────────────────────────────────────────
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  py_cell_ids <- reticulate::r_to_py(as.list(cell_ids))

  message("Computing pre-shape embedding ...")
  kendall_tpca$compute_pre_shape_embedding(
    pre_shape_output_dir   = work_dir,
    df                     = reticulate::r_to_py(df),
    num_vertices_to_sample = as.integer(num_vertices),
    cell_ids_to_analyze    = py_cell_ids,
    x_vertex_col           = "x",
    y_vertex_col           = "y",
    cell_id_col            = "cell_id"
  )

  message("Running TPCA ...")
  kendall_tpca$run_kendall_tpca(
    pre_shape_input_dir   = work_dir,
    output_dir            = work_dir,
    cell_ids_to_analyze   = py_cell_ids,
    cell_id_col           = "cell_id",
    max_frechet_mean_iter = as.integer(max_frechet_iter),
    eta                   = eta,
    use_parallel          = use_parallel,
    num_threads           = as.integer(num_threads),
    frechet_mean_tol      = frechet_mean_tol
  )

  if (!return_results) {
    message("Result not loaded into memory (return_results = FALSE).")
    message(sprintf("  TPCA output written to: %s", work_dir))
    message(sprintf(
      "  Load with load_kstitch_results(\"%s\", type = \"tpca\")", work_dir
    ))
    return(list(output_path = work_dir))
  }

  # ── 6. Load and return results ───────────────────────────────────────────────
  message("Reading TPCA results ...")
  result              <- .load_kendall_tpca_output(work_dir)
  result$contour_type <- contour_type
  result$output_dir   <- if (resolved$is_temp) NULL else work_dir
  result$group        <- "all"
  result$is_groupwise <- FALSE

  list(all = result)
}


# ── Internal helpers ───────────────────────────────────────────────────────────

#' @keywords internal
.extract_contour_df <- function(seg_list, contour_type) {

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

  col_map <- list(
    cell_id = c("cell", "cell_id"),
    x       = "x",
    y       = "y"
  )

  df <- .normalise_columns(df, col_map)
  df[, c("x", "y", "cell_id")]
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
