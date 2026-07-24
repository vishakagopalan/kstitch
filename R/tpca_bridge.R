# ---- private helpers --------------------------------------------------------

#' @noRd
.load_kendall_tpca_output <- function(output_dir, load_pre_shape = TRUE) {

  kendall_tpca_path <- file.path(output_dir, "TPCA_Info.h5")
  meta_path         <- file.path(output_dir, "Shape_Metadata.csv.gz")
  pre_shape_path    <- file.path(output_dir, "Pre_Shape_Space_Embedding.h5")

  if (!file.exists(kendall_tpca_path))
    stop(sprintf("TPCA_Info.h5 not found in '%s'. Has run_kendall_tpca() completed?", output_dir))
  if (!file.exists(meta_path))
    stop(sprintf("Shape_Metadata.csv.gz not found in '%s'.", output_dir))

  kendall_tpca_mat   <- t(rhdf5::h5read(kendall_tpca_path, name = "embedding"))
  kendall_tpca_idxes <- rhdf5::h5read(kendall_tpca_path, name = "processed_idxes")
  frechet_mean       <- rhdf5::h5read(kendall_tpca_path, name = "frechet_mean")
  v_matrix           <- rhdf5::h5read(kendall_tpca_path, name = "v_matrix")
  variances          <- rhdf5::h5read(kendall_tpca_path, name = "variances")

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
    dplyr::inner_join(dplyr::select(meta, numpy_idx, cell_id), by = "numpy_idx") |>
    dplyr::select(-numpy_idx) |>
    tibble::column_to_rownames("cell_id")

  kendall_tpca_mat <- as.matrix(kendall_tpca_df)

  pre_shape <- NULL
  if (load_pre_shape && file.exists(pre_shape_path)) {
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

#' @noRd
.resolve_tpca_output_dir <- function(output_dir, env = parent.frame()) {
  if (is.null(output_dir)) {
    path <- withr::local_tempdir(.local_envir = env)
    return(list(path = path, is_temp = TRUE))
  }
  path <- path.expand(output_dir)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  list(path = path, is_temp = FALSE)
}

#' @noRd
.random_id <- function() {
  paste0(sample(c(letters, 0:9), 8, replace = TRUE), collapse = "")
}

#' @noRd
.tpca_fields_msg <- c(
  "  $TPCA_Embedding  matrix, cells x shape PCs",
  "  $Info            list: variances, v_matrix, frechet_mean, pre_shape_embedding",
  "  $Metadata        data frame: per-cell shape metadata including Frobenius norm (scale)",
  "  $contour_type    character: 'cell' or 'nucleus'"
)

#' @noRd
.emit_tpca_deferred_msg <- function(output_paths) {
  first_name <- names(output_paths)[1]
  paths_lines <- paste0(
    "  ", format(names(output_paths), justify = "left"),
    ": ", output_paths,
    collapse = "\n"
  )
  message(paste(c(
    "Results saved to disk. Load with:",
    sprintf("  load_kstitch_results(output_paths[[\"%s\"]], type = \"tpca\")", first_name),
    "",
    "Results written to:",
    paths_lines,
    "",
    "Returned list contains:",
    .tpca_fields_msg
  ), collapse = "\n"))
}

#' @noRd
.run_tpca_single <- function(boundary_parquet_path,
                             work_dir,
                             cell_ids,
                             contour_type,
                             num_vertices,
                             cell_id_col,
                             x_col,
                             y_col,
                             eta,
                             use_parallel,
                             num_threads,
                             frechet_mean_tol,
                             max_frechet_iter,
                             kendall_tpca) {
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

  # run_kendall_tpca reads Shape_Metadata.csv.gz which always has a column
  # named "cell_id" (hardcoded in compute_pre_shape_embedding), regardless
  # of what cell_id_col was in the source parquet. Pass cell_ids_to_analyze
  # with cell_id_col fixed to "cell_id" so subsetting works correctly.
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
}

# ---- exported functions -----------------------------------------------------

#' Load kstitch results from disk
#'
#' Unified loader for results previously saved by \code{\link{run_tpca}},
#' \code{\link{compute_nmf}}, or \code{\link{link_shape_and_factors}} when
#' called with \code{return_results = FALSE}.
#'
#' @param path For \code{type = "tpca"}: path to the output directory
#'   containing \code{TPCA_Info.h5} and \code{Shape_Metadata.csv.gz}.
#'   For \code{type = "nmf"} or \code{type = "cca"}: path to an RDS file.
#' @param type One of \code{"tpca"}, \code{"nmf"}, or \code{"cca"}.
#'
#' @return For \code{"tpca"}: a list with \code{TPCA_Embedding},
#'   \code{Info} (variances, v_matrix, frechet_mean, pre_shape_embedding),
#'   \code{Metadata}, and \code{contour_type}.
#'   For \code{"nmf"}: a list with \code{H}, \code{W}, and \code{stability}.
#'   For \code{"cca"}: a list with \code{CSP_Scores}, \code{CEP_Scores},
#'   \code{CC_Corr_Coefs}, \code{corr.shape.with.csp},
#'   \code{corr.shape.with.cep}, \code{corr.exp.with.csp}, and
#'   \code{corr.exp.with.cep}.
#'
#' @param load_pre_shape Logical. For \code{type = "tpca"} only: whether to
#'   load the pre-shape embedding from \code{Pre_Shape_Space_Embedding.h5}
#'   into \code{Info$pre_shape_embedding}. Default \code{TRUE}. Set to
#'   \code{FALSE} to skip loading this (potentially large) array when it is
#'   not needed — e.g. before storing results in a Seurat object.
#'
#' @seealso \code{\link{run_tpca}}, \code{\link{compute_nmf}},
#'   \code{\link{link_shape_and_factors}}
#' @export
load_kstitch_results <- function(path, type = c("tpca", "nmf", "cca"),
                                 load_pre_shape = TRUE) {
  type <- match.arg(type)
  if (type == "tpca") {
    result <- .load_kendall_tpca_output(path, load_pre_shape = load_pre_shape)
    return(list(all = result))

  } else {
    if (!file.exists(path))
      stop(sprintf("File not found: '%s'", path))
    readRDS(path)
  }
}

#' Run Tangent PCA on cell or nucleus boundary contours
#'
#' Calls \code{compute_pre_shape_embedding()} and \code{run_kendall_tpca()}
#' from \code{kendall_tpca.py} via reticulate, then reads the results back
#' into R.
#'
#' @param boundary_parquet_path Path to a parquet file containing boundary
#'   coordinates.
#' @param output_dir Directory where intermediate/output files will be
#'   written. Optional. If \code{NULL} (the default), a temporary directory
#'   is used. When \code{return_results = TRUE} the temp directory is deleted
#'   after the call; when \code{return_results = FALSE} it is kept so the
#'   saved files remain accessible. An explicit path is created if absent and
#'   always left on disk.
#' @param cell_ids Optional character vector of cell IDs to analyse. Ignored
#'   when \code{cell_groups} is supplied.
#' @param cell_groups Optional named character vector mapping cell IDs to
#'   group labels (e.g. \code{c(cell1 = "typeA", cell2 = "typeB")}). When
#'   supplied, TPCA is run separately per group and a named list of per-group
#'   results is returned. Each group is written to
#'   \code{<output_dir>/<group>_<contour_type>/}.
#' @param contour_type One of \code{"cell"} or \code{"nucleus"}.
#' @param num_vertices Integer. Vertices to resample each contour to.
#'   Default 50L.
#' @param cell_id_col Column name for cell IDs. Default \code{"cell_id"}.
#' @param x_col Column name for x coordinates. Default \code{"vertex_x"}.
#' @param y_col Column name for y coordinates. Default \code{"vertex_y"}.
#' @param eta Learning rate for Fréchet mean iteration. Default 1.
#' @param use_parallel Logical. Use multiprocessing in Python. Default
#'   \code{FALSE}.
#' @param num_threads Integer. Parallel threads when
#'   \code{use_parallel = TRUE}. Default 8L.
#' @param frechet_mean_tol Convergence tolerance. Default 1e-4.
#' @param max_frechet_iter Maximum iterations. Default 1000L.
#' @param load_pre_shape Logical. Whether to load the pre-shape embedding
#'   into \code{Info$pre_shape_embedding}. Default \code{TRUE}. Set to
#'   \code{FALSE} to skip loading this (potentially large) array — e.g.
#'   when results will be stored in a Seurat object via
#'   \code{\link{store_tpca_results}}.
#' @param return_results Logical. If \code{TRUE} (default), results are
#'   loaded into R and returned. If \code{FALSE}, results are written to
#'   disk and only the output paths are returned, with a message showing
#'   how to load them via \code{\link{load_kstitch_results}}.
#'
#' @return When \code{return_results = TRUE}: a list (or named list of lists
#'   when \code{cell_groups} is supplied) with \code{TPCA_Embedding},
#'   \code{Info}, \code{Metadata}, \code{contour_type}, and
#'   \code{output_dir}.
#'   When \code{return_results = FALSE}: a list with \code{output_paths}
#'   (named character vector) and \code{fields} describing the structure
#'   returned by \code{load_kstitch_results(..., type = "tpca")}.
#' @seealso \code{\link{load_kstitch_results}}, \code{\link{store_tpca_results}}
#' @export
run_tpca <- function(boundary_parquet_path,
                     output_dir       = NULL,
                     cell_ids         = NULL,
                     cell_groups      = NULL,
                     contour_type     = c("cell", "nucleus"),
                     num_vertices     = 50L,
                     cell_id_col      = "cell_id",
                     x_col            = "vertex_x",
                     y_col            = "vertex_y",
                     eta              = 1,
                     use_parallel     = FALSE,
                     num_threads      = 8L,
                     frechet_mean_tol = 1e-4,
                     max_frechet_iter = 1000L,
                     return_results   = TRUE,
                     load_pre_shape   = TRUE) {

  contour_type          <- match.arg(contour_type)
  boundary_parquet_path <- path.expand(boundary_parquet_path)

  if (!file.exists(boundary_parquet_path))
    stop(sprintf("boundary_parquet_path '%s' does not exist.", boundary_parquet_path))

  # When return_results = FALSE we must keep the output dir alive after the
  # call returns (the files are the point). Bypass withr's auto-cleanup by
  # allocating a plain temp dir ourselves; the OS will reclaim it eventually.
  if (!return_results && is.null(output_dir)) {
    work_dir    <- file.path(tempdir(), paste0("kstitch_tpca_", .random_id()))
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    resolved    <- list(path = work_dir, is_temp = FALSE)
  } else {
    resolved <- .resolve_tpca_output_dir(output_dir)
  }
  work_dir <- resolved$path

  kendall_tpca_py_dir <- system.file("python", package = "kstitch")
  if (!nzchar(kendall_tpca_py_dir))
    stop("Could not locate inst/python/ inside the kstitch package.")

  kendall_tpca <- reticulate::import_from_path(
    "kendall_tpca", path = kendall_tpca_py_dir
  )

  shared_args <- list(
    boundary_parquet_path = boundary_parquet_path,
    contour_type          = contour_type,
    num_vertices          = num_vertices,
    cell_id_col           = cell_id_col,
    x_col                 = x_col,
    y_col                 = y_col,
    eta                   = eta,
    use_parallel          = use_parallel,
    num_threads           = num_threads,
    frechet_mean_tol      = frechet_mean_tol,
    max_frechet_iter      = max_frechet_iter,
    kendall_tpca          = kendall_tpca
  )

  # ---- groupwise -----------------------------------------------------------
  if (!is.null(cell_groups)) {
    groups       <- unique(cell_groups)
    output_paths <- stats::setNames(
      file.path(work_dir, paste0(groups, "_", contour_type)),
      groups
    )

    results <- lapply(groups, function(grp) {
      grp_cells <- names(cell_groups)[cell_groups == grp]
      grp_dir   <- output_paths[[grp]]
      dir.create(grp_dir, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("--- Group: %s ---", grp))
      do.call(.run_tpca_single, c(shared_args, list(
        cell_ids = grp_cells,
        work_dir = grp_dir
      )))
      if (return_results) {
        message("Reading TPCA results ...")
        res               <- .load_kendall_tpca_output(grp_dir,
                                                       load_pre_shape = load_pre_shape)
        res$contour_type  <- contour_type
        res$output_dir    <- grp_dir
        res$group         <- grp
        res$is_groupwise  <- TRUE
        res
      } else {
        NULL
      }
    })

    if (!return_results) {
      .emit_tpca_deferred_msg(output_paths)
      return(list(
        output_paths = output_paths,
        fields = list(
          TPCA_Embedding = "matrix, cells x shape PCs",
          Info           = "list: variances, v_matrix, frechet_mean, pre_shape_embedding",
          Metadata       = "data frame: per-cell shape metadata including Frobenius norm (scale)",
          contour_type   = "character: 'cell' or 'nucleus'"
        )
      ))
    }

    return(stats::setNames(results, groups))
  }

  # ---- single group --------------------------------------------------------
  do.call(.run_tpca_single, c(shared_args, list(
    cell_ids = cell_ids,
    work_dir = work_dir
  )))

  if (!return_results) {
    output_paths <- stats::setNames(work_dir, "all")
    .emit_tpca_deferred_msg(output_paths)
    return(list(
      output_paths = output_paths,
      fields = list(
        TPCA_Embedding = "matrix, cells x shape PCs",
        Info           = "list: variances, v_matrix, frechet_mean, pre_shape_embedding",
        Metadata       = "data frame: per-cell shape metadata including Frobenius norm (scale)",
        contour_type   = "character: 'cell' or 'nucleus'"
      )
    ))
  }

  message("Reading TPCA results ...")
  result              <- .load_kendall_tpca_output(work_dir,
                                                   load_pre_shape = load_pre_shape)
  result$contour_type <- contour_type
  result$output_dir   <- if (resolved$is_temp) NULL else work_dir
  result$group        <- "all"
  result$is_groupwise <- FALSE
  list(all = result)
}


#' Store TPCA results in a Seurat object
#'
#' Writes the TPCA embedding as a \code{DimReduc} slot and fit metadata to
#' \code{obj@misc$kstitch$tpca}.
#'
#' For a single-group result, the reduction is named
#' \code{tpca_<contour_type>} and metadata stored under
#' \code{obj@misc$kstitch$tpca[[contour_type]]}.
#'
#' For a groupwise result (named list of per-group results), the reduction
#' for each group is named \code{tpca_<contour_type>_<group>} and metadata
#' stored under \code{obj@misc$kstitch$tpca[[contour_type]][[group]]}.
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

  for (grp in names(tpca_result)) {
    grp_result <- tpca_result[[grp]]
    reduction_name <- if (isTRUE(grp_result$is_groupwise)) {
      .sanitize_reduction_key(
        paste0("tpca_", grp_result$contour_type, "_", grp)
      )
    } else {
      paste0("tpca_", grp_result$contour_type)
    }
    obj <- .store_tpca_single(
      obj            = obj,
      result         = grp_result,
      reduction_name = reduction_name,
      key_prefix     = reduction_key_prefix,
      misc_group     = if (isTRUE(grp_result$is_groupwise)) grp else NULL
    )
  }
  obj
}

#' @noRd
.sanitize_reduction_key <- function(key) {
  gsub("_+", "_", gsub("[^A-Za-z0-9_]", "_", key))
}

#' @noRd
.store_tpca_single <- function(obj, result, reduction_name, key_prefix,
                               misc_group = NULL) {
  emb       <- result$TPCA_Embedding
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
    key        = key_prefix,
    assay      = Seurat::DefaultAssay(obj)
  )

  if (is.null(obj@misc$kstitch))      obj@misc$kstitch      <- list()
  if (is.null(obj@misc$kstitch$tpca)) obj@misc$kstitch$tpca <- list()

  contour_type <- result$contour_type
  if (is.null(misc_group)) {
    obj@misc$kstitch$tpca[[contour_type]] <- result$Info
  } else {
    if (is.null(obj@misc$kstitch$tpca[[contour_type]]))
      obj@misc$kstitch$tpca[[contour_type]] <- list()
    obj@misc$kstitch$tpca[[contour_type]][[misc_group]] <- result$Info
  }

  obj
}


#' Retrieve TPCA fit metadata from a Seurat object
#'
#' @param obj A Seurat v5 object.
#' @param contour_type \code{"cell"} or \code{"nucleus"}. \code{NULL}
#'   returns the full tpca misc list.
#' @param group Optional group label. When supplied, returns the per-group
#'   fit metadata for the given contour type.
#'
#' @return A list of TPCA fit metadata.
#' @seealso \code{\link{store_tpca_results}}
#' @export
get_tpca_results <- function(obj, contour_type = NULL, group = NULL) {

  tpca_misc <- obj@misc$kstitch$tpca

  if (is.null(tpca_misc))
    stop("No TPCA results found in obj. Run store_tpca_results() first.")

  if (is.null(contour_type)) return(tpca_misc)

  if (!contour_type %in% names(tpca_misc))
    stop(sprintf(
      "contour_type '%s' not found. Available: %s.",
      contour_type, paste(names(tpca_misc), collapse = ", ")
    ))

  ct_misc <- tpca_misc[[contour_type]]

  if (is.null(group)) return(ct_misc)

  if (!group %in% names(ct_misc))
    stop(sprintf(
      "group '%s' not found for contour_type '%s'. Available: %s.",
      group, contour_type, paste(names(ct_misc), collapse = ", ")
    ))

  ct_misc[[group]]
}
