# ── Cache path helpers ─────────────────────────────────────────────────────────

#' Return the root kstitch TPCA cache directory
#' @keywords internal
.tpca_cache_dir <- function() {
  tools::R_user_dir("kstitch", "cache")
}

#' Build the cache key from contour data and run parameters
#'
#' The key is an SHA-256 hash of the contour data (either raw bytes from a
#' parquet file via \code{readBin}, or an in-memory data frame) combined with
#' the run parameters that affect the output.
#'
#' @param contour_data Either a raw vector (from \code{readBin} on a parquet
#'   file) or a data frame of contour coordinates.
#' @keywords internal
.tpca_cache_key <- function(contour_data, contour_type, num_vertices, eta,
                             frechet_mean_tol, max_frechet_iter) {
  payload <- list(
    contour_data     = contour_data,
    contour_type     = contour_type,
    num_vertices     = num_vertices,
    eta              = eta,
    frechet_mean_tol = frechet_mean_tol,
    max_frechet_iter = max_frechet_iter
  )
  digest::digest(payload, algo = "sha256")
}

#' Return the cache directory for a specific hash key
#' @keywords internal
.tpca_cache_entry_dir <- function(key) {
  file.path(.tpca_cache_dir(), key)
}

#' Check whether a cache entry exists and is complete
#' @keywords internal
.tpca_cache_exists <- function(key) {
  entry_dir <- .tpca_cache_entry_dir(key)
  all(file.exists(file.path(entry_dir, c(
    "TPCA_Info.h5",
    "Shape_Metadata.csv.gz",
    "Pre_Shape_Space_Embedding.h5",
    "cache_meta.rds"
  ))))
}


# ── Cache read / write ─────────────────────────────────────────────────────────

#' Write TPCA outputs to the persistent cache
#'
#' @param key        Cache key (SHA-256 string).
#' @param output_dir Directory containing the TPCA run outputs to cache.
#' @param cache_meta Named list of human-readable metadata to store alongside
#'                   the outputs (contour_type, parameters, timestamp, etc.).
#' @keywords internal
.tpca_cache_write <- function(key, output_dir, cache_meta) {
  entry_dir <- .tpca_cache_entry_dir(key)
  dir.create(entry_dir, recursive = TRUE, showWarnings = FALSE)

  files_to_cache <- c(
    "TPCA_Info.h5",
    "Shape_Metadata.csv.gz",
    "Pre_Shape_Space_Embedding.h5"
  )
  for (f in files_to_cache) {
    src <- file.path(output_dir, f)
    if (file.exists(src))
      file.copy(src, file.path(entry_dir, f), overwrite = TRUE)
  }

  cache_meta$cached_at <- Sys.time()
  saveRDS(cache_meta, file.path(entry_dir, "cache_meta.rds"))
  invisible(entry_dir)
}

#' Read TPCA outputs from the persistent cache
#'
#' @param key Cache key (SHA-256 string).
#' @return A list in the format returned by \code{store_tpca_results()}.
#' @keywords internal
.tpca_cache_read <- function(key) {
  entry_dir  <- .tpca_cache_entry_dir(key)
  cache_meta <- readRDS(file.path(entry_dir, "cache_meta.rds"))
  result     <- .load_kendall_tpca_output(entry_dir)
  result$contour_type <- cache_meta$contour_type
  result$output_dir   <- entry_dir
  result
}


# ── Public API ─────────────────────────────────────────────────────────────────

#' List cached TPCA runs
#'
#' Displays all entries currently held in the kstitch TPCA cache, with
#' their contour type, key parameters, cache date, and disk size.
#'
#' @return A data frame with one row per cache entry, invisibly. Also prints
#'   a formatted summary to the console.
#' @export
list_tpca_cache <- function() {
  cache_root <- .tpca_cache_dir()
  if (!dir.exists(cache_root)) {
    message("No TPCA cache found at ", cache_root)
    return(invisible(data.frame()))
  }

  entries <- list.dirs(cache_root, full.names = TRUE, recursive = FALSE)
  if (length(entries) == 0L) {
    message("TPCA cache is empty.")
    return(invisible(data.frame()))
  }

  rows <- lapply(entries, function(entry_dir) {
    meta_path <- file.path(entry_dir, "cache_meta.rds")
    if (!file.exists(meta_path)) return(NULL)
    meta <- readRDS(meta_path)
    size_mb <- sum(
      file.size(list.files(entry_dir, full.names = TRUE)), na.rm = TRUE
    ) / 1024^2

    data.frame(
      key              = basename(entry_dir),
      contour_type     = meta$contour_type     %||% NA_character_,
      num_vertices     = meta$num_vertices     %||% NA_integer_,
      eta              = meta$eta              %||% NA_real_,
      frechet_mean_tol = meta$frechet_mean_tol %||% NA_real_,
      max_frechet_iter = meta$max_frechet_iter %||% NA_integer_,
      cached_at        = format(meta$cached_at),
      size_mb          = round(size_mb, 2L),
      stringsAsFactors = FALSE
    )
  })

  rows <- do.call(rbind, Filter(Negate(is.null), rows))
  print(rows, row.names = FALSE)
  invisible(rows)
}

#' Clear one or all TPCA cache entries
#'
#' @param key Character. One or more SHA-256 cache keys to remove, as shown
#'   by \code{list_tpca_cache()}. If \code{NULL} (the default), \emph{all}
#'   cache entries are removed after a confirmation prompt.
#' @param confirm Logical. If \code{TRUE} (the default), prompts for
#'   confirmation before deleting all entries. Set to \code{FALSE} for
#'   non-interactive use.
#'
#' @return Invisibly returns the paths that were removed.
#' @export
clear_tpca_cache <- function(key = NULL, confirm = TRUE) {
  cache_root <- .tpca_cache_dir()

  if (!dir.exists(cache_root)) {
    message("No TPCA cache found — nothing to clear.")
    return(invisible(character()))
  }

  if (is.null(key)) {
    entries <- list.dirs(cache_root, full.names = TRUE, recursive = FALSE)
    if (length(entries) == 0L) {
      message("TPCA cache is already empty.")
      return(invisible(character()))
    }
    if (confirm && interactive()) {
      answer <- readline(sprintf(
        "Delete all %d cache entries in %s? [y/N]: ",
        length(entries), cache_root
      ))
      if (!tolower(trimws(answer)) %in% c("y", "yes")) {
        message("Aborted.")
        return(invisible(character()))
      }
    }
    targets <- entries
  } else {
    targets <- file.path(cache_root, key)
    missing <- targets[!dir.exists(targets)]
    if (length(missing) > 0L)
      warning("Cache entries not found: ",
              paste(basename(missing), collapse = ", "))
    targets <- targets[dir.exists(targets)]
  }

  unlink(targets, recursive = TRUE)
  message("Removed ", length(targets), " cache ",
          if (length(targets) == 1L) "entry." else "entries.")
  invisible(targets)
}


# ── Null-coalescing helper (internal) ─────────────────────────────────────────
`%||%` <- function(x, y) if (is.null(x)) y else x
