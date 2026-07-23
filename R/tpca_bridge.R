# TODO (PR #6b): Add contour export functions that extract cell/nucleus polygon
# coordinates directly from Seurat v5 objects created by ReadXenium /
# ReadCosMx, write them to a parquet file, and call run_tpca() internally.
# For now, the user supplies a pre-built boundary parquet file.

# ---- private helpers --------------------------------------------------------

#' @noRd
.load_kendall_tpca_output <- function(output_dir) {

  kendall_tpca_path       <- file.path(output_dir, "TPCA_Info.h5")
  meta_path      <- file.path(output_dir, "Shape_Metadata.csv.gz")
  pre_shape_path <- file.path(output_dir, "Pre_Shape_Space_Embedding.h5")

  if (!file.exists(kendall_tpca_path))
    stop(sprintf("TPCA_Info.h5 not found in '%s'. Has run_kendall_tpca() completed?", output_dir))
  if (!file.exists(meta_path))
    stop(sprintf("Shape_Metadata.csv.gz not found in '%s'.", output_dir))

  kendall_tpca_mat      <- t(rhdf5::h5read(kendall_tpca_path, name = "embedding"))
  kendall_tpca_idxes    <- rhdf5::h5read(kendall_tpca_path, name = "processed_idxes")
  frechet_mean <- rhdf5::h5read(kendall_tpca_path, name = "frechet_mean")
  v_matrix     <- rhdf5::h5read(kendall_tpca_path, name = "v_matrix")
  variances    <- rhdf5::h5read(kendall_tpca_path, name = "variances")

  dimnames(kendall_tpca_mat) <- list(
    as.character(kendall_tpca_idxes),
    paste0("Shape_PC", seq_len(ncol(kendall_tpca_mat)))
  )

  meta <- readr::read_csv(meta_path, show_col_types = FALSE)
  meta <- meta |>
    dplyr::mutate(numpy_idx = as.character(numpy_idx))

  kendall_tpca_df <- as.data.frame(kendall_tpca_mat) |>
    tibble::rownames_to_column("numpy_idx") |>
    dplyr::mutate(numpy_idx = as.character(numpy_idx)) |>
    dplyr::inner_join(dplyr::select(meta,numpy_idx, cell_id), by = "numpy_idx") |>
    dplyr::select(-numpy_idx) |>
    tibble::column_to_rownames("cell_id")

  kendall_tpca_mat <- as.matrix(kendall_tpca_df)

  pre_shape <- NULL
  if (file.exists(pre_shape_path)) {
    raw       <- rhdf5::h5read(pre_shape_path, name = "pre_shape_space_embedding")
    pre_shape <- aperm(raw, perm = c(3L, 2L, 1L))
  }

  list(
    TPCA_Embedding = kendall_tpca_mat,
    Info = list(
      variances           = variances,
      v_matrix            = v_matrix,
      frechet_mean        = frechet_mean,
      pre_shape_embedding = pre_shape
    ),
    Metadata = meta
  )
}

#' Resolve output_dir, falling back to an auto-cleaned temp directory
#'
#' If \code{output_dir} is \code{NULL}, a fresh temporary directory is
#' created and scheduled for deletion when the calling function exits
#' (success or error), via \code{withr::local_tempdir()} in the caller's
#' environment. If \code{output_dir} is supplied, it is expanded and created
#' if absent, and is never auto-deleted -- the caller is responsible for it.
#'
#' @param output_dir Character or \code{NULL}.
#' @param env Environment to scope the temp directory's lifetime to
#'   (normally the calling function's environment, via
#'   \code{parent.frame()}).
#' @return A list with \code{path} (character, expanded/created directory)
#'   and \code{is_temp} (logical, whether it was auto-generated).
#' @keywords internal
.resolve_tpca_output_dir <- function(output_dir, env = parent.frame()) {
  if (is.null(output_dir)) {
    path <- withr::local_tempdir(.local_envir = env)
    return(list(path = path, is_temp = TRUE))
  }
  path <- path.expand(output_dir)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  list(path = path, is_temp = FALSE)
}

# ---- exported functions -----------------------------------------------------

#' Run Tangent PCA on cell or nucleus boundary contours
#'
#' Calls \code{compute_pre_shape_embedding()} and \code{run_kendall_tpca()} from
#' \code{kendall_tpca.py} via reticulate, then reads the results back into R.
#'
#' @param boundary_parquet_path Path to a parquet file containing boundary
#'   coordinates.
#' @param output_dir Directory where intermediate/output files will be
#'   written. Optional. If \code{NULL} (the default), a temporary directory
#'   is created for the duration of the call and deleted automatically when
#'   \code{run_tpca()} returns -- the returned list already contains
#'   everything downstream code needs, so nothing is lost. If supplied, the
#'   directory is created if absent, expanded (\code{~} is resolved), and
#'   left on disk after the call returns for later inspection or reuse.
#' @param cell_ids Optional character vector of cell IDs to analyse.
#' @param contour_type One of \code{"cell"} or \code{"nucleus"}.
#' @param num_vertices Integer. Vertices to resample each contour to. Default
#'   50L.
#' @param cell_id_col Column name for cell IDs. Default \code{"cell_id"}.
#' @param x_col Column name for x coordinates. Default \code{"vertex_x"}.
#' @param y_col Column name for y coordinates. Default \code{"vertex_y"}.
#' @param eta Learning rate for Fréchet mean iteration. Default 1.
#' @param use_parallel Logical. Use multiprocessing in Python. Default
#'   \code{FALSE}.
#' @param num_threads Integer. Parallel threads when \code{use_parallel = TRUE}.
#'   Default 8L.
#' @param frechet_mean_tol Convergence tolerance. Default 1e-4.
#' @param max_frechet_iter Maximum iterations. Default 1000L.
#'
#' @return A list with \code{TPCA_Embedding}, \code{Info}, \code{contour_type},
#'   and \code{output_dir} (the directory actually used; \code{NULL} if a
#'   temporary directory was used and has already been deleted).
#' @seealso \code{\link{store_tpca_results}}
#' @export
run_tpca <- function(boundary_parquet_path,
                     output_dir       = NULL,
                     cell_ids         = NULL,
                     contour_type     = c("cell", "nucleus"),
                     num_vertices     = 50L,
                     cell_id_col      = "cell_id",
                     x_col            = "vertex_x",
                     y_col            = "vertex_y",
                     eta              = 1,
                     use_parallel     = FALSE,
                     num_threads      = 8L,
                     frechet_mean_tol = 1e-4,
                     max_frechet_iter = 1000L) {

  contour_type          <- match.arg(contour_type)
  boundary_parquet_path <- path.expand(boundary_parquet_path)

  if (!file.exists(boundary_parquet_path))
    stop(sprintf("boundary_parquet_path '%s' does not exist.", boundary_parquet_path))

  resolved   <- .resolve_tpca_output_dir(output_dir)
  work_dir   <- resolved$path

  kendall_tpca_py_dir <- system.file("python", package = "kstitch")
  if (!nzchar(kendall_tpca_py_dir))
    stop("Could not locate inst/python/ inside the kstitch package.")

  kendall_tpca <- reticulate::import_from_path("kendall_tpca", path = kendall_tpca_py_dir)

  py_cell_ids <- if (!is.null(cell_ids))
    reticulate::r_to_py(as.list(cell_ids)) else NULL

  message("Computing pre-shape embedding ...")
  kendall_tpca$compute_pre_shape_embedding(
    boundary_parquet_path  = boundary_parquet_path,
    pre_shape_output_dir   = work_dir,
    num_vertices_to_sample = as.integer(num_vertices),
    cell_ids_to_analyze    = py_cell_ids,
    x_vertex_col           = x_col,
    y_vertex_col           = y_col,
    cell_id_col            = cell_id_col
  )

  message("Running TPCA ...")
  kendall_tpca$run_kendall_tpca(
    pre_shape_input_dir   = work_dir,
    output_dir            = work_dir,
    cell_ids_to_analyze   = py_cell_ids,
    cell_id_col           = cell_id_col,
    max_frechet_mean_iter = as.integer(max_frechet_iter),
    eta                   = eta,
    use_parallel          = use_parallel,
    num_threads           = as.integer(num_threads),
    frechet_mean_tol      = frechet_mean_tol
  )

  message("Reading TPCA results ...")
  result              <- .load_kendall_tpca_output(work_dir)
  result$contour_type <- contour_type
  # If a temp dir was used, it will be deleted (via withr) when this
  # function returns -- report NULL rather than a path that no longer
  # exists by the time the caller can act on it.
  result$output_dir   <- if (resolved$is_temp) NULL else work_dir

  result
}


#' Store TPCA results in a Seurat object
#'
#' Writes the TPCA embedding as a \code{DimReduc} slot and fit metadata to
#' \code{obj@misc$kstitch$tpca[[contour_type]]}.
#'
#' @param obj A Seurat v5 object.
#' @param tpca_result The list returned by \code{\link{run_tpca}}.
#' @param reduction_key_prefix Key prefix for the \code{DimReduc}. Default
#'   \code{"ShapePC_"}.
#'
#' @return The Seurat object with TPCA results stored.
#' @seealso \code{\link{run_tpca}}, \code{\link{get_tpca_results}}
#' @export
store_tpca_results <- function(obj,
                               tpca_result,
                               reduction_key_prefix = "ShapePC_") {

  contour_type   <- tpca_result$contour_type
  reduction_name <- paste0("tpca_", contour_type)

  emb       <- tpca_result$TPCA_Embedding
  obj_cells <- rownames(obj@meta.data)
  missing   <- setdiff(obj_cells, rownames(emb))

  if (length(missing) > 0) {
    warning(sprintf(
      "%d cells in obj have no TPCA embedding and will have NA scores: %s%s",
      length(missing),
      paste(head(missing, 5), collapse = ", "),
      if (length(missing) > 5) " ..." else ""
    ))
    na_rows <- matrix(NA_real_, nrow = length(missing), ncol = ncol(emb),
                      dimnames = list(missing, colnames(emb)))
    emb     <- rbind(emb, na_rows)[obj_cells, , drop = FALSE]
  } else {
    emb <- emb[intersect(obj_cells, rownames(emb)), , drop = FALSE]
  }

  obj[[reduction_name]] <- Seurat::CreateDimReducObject(
    embeddings = emb,
    key        = reduction_key_prefix,
    assay      = Seurat::DefaultAssay(obj)
  )

  if (is.null(obj@misc$kstitch))       obj@misc$kstitch       <- list()
  if (is.null(obj@misc$kstitch$tpca))  obj@misc$kstitch$tpca  <- list()

  obj@misc$kstitch$tpca[[contour_type]] <- tpca_result$Info

  obj
}


#' Retrieve TPCA fit metadata from a Seurat object
#'
#' @param obj A Seurat v5 object.
#' @param contour_type \code{"cell"} or \code{"nucleus"}. \code{NULL} returns
#'   the full tpca misc list.
#'
#' @return A list of TPCA fit metadata.
#' @seealso \code{\link{store_tpca_results}}
#' @export
get_tpca_results <- function(obj, contour_type = NULL) {

  tpca_misc <- obj@misc$kstitch$tpca

  if (is.null(tpca_misc))
    stop("No TPCA results found in obj. Run store_tpca_results() first.")

  if (is.null(contour_type)) return(tpca_misc)

  if (!contour_type %in% names(tpca_misc))
    stop(sprintf(
      "contour_type '%s' not found. Available: %s.",
      contour_type,
      paste(names(tpca_misc), collapse = ", ")
    ))

  tpca_misc[[contour_type]]
}
