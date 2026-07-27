#' Store NMF results in a Seurat object
#'
#' Takes the list returned by \code{\link{compute_nmf}} and writes the fit into
#' the Seurat object as a \code{DimReduc} slot named
#' \code{obj[["<reduction_prefix><reduction_suffix>"]]}. Feature loadings are
#' stored in the reduction's \code{@@feature.loadings}; fit metadata
#' (\code{Fit_Error}, \code{Fit_Summary}, \code{Factor_Gene_List}) is stored
#' in the reduction's \code{@@misc}.
#'
#' For multi-group workflows, call this function once per group in an external
#' loop, passing a distinct \code{reduction_suffix} each time (e.g. the cell
#' type name), so each group gets its own \code{DimReduc}.
#'
#' Cells in \code{obj} not present in the NMF result receive \code{NA}
#' embeddings — this guards against partial or mismatched inputs.
#'
#' If \code{NMF_Matrix} is \code{NULL} (e.g. the group was skipped for having
#' too few cells), a message is emitted and no \code{DimReduc} is created.
#'
#' @param obj A Seurat v5 object.
#' @param nmf_result List as returned by \code{\link{compute_nmf}}.
#' @param reduction_suffix Character. Appended to \code{reduction_prefix} to
#'   form the full reduction name, e.g. \code{"Keratinocytes"} →
#'   \code{"nmf_Keratinocytes"}. Default \code{"all"}.
#' @param reduction_prefix Character. Prefix for the reduction name.
#'   Default \code{"nmf_"}.
#' @param reduction_key_prefix Character. Key prefix passed to
#'   \code{Seurat::CreateDimReducObject()}. The suffix (with non-alphanumeric
#'   characters removed) is appended to keep keys unique across groups, e.g.
#'   \code{"NMFKeratinocytes_"}. Default \code{"NMF"}.
#'
#' @return The Seurat object with one new \code{DimReduc} added.
#' @seealso \code{\link{compute_nmf}}, \code{\link{get_nmf_results}}
#' @export
store_nmf_results <- function(obj,
                              nmf_result,
                              reduction_suffix     = "all",
                              reduction_prefix     = "nmf_",
                              reduction_key_prefix = "NMF") {

  if (is.null(nmf_result$NMF_Matrix)) {
    message(sprintf(
      "Skipping '%s%s': no NMF fit available (NMF_Matrix is NULL).",
      reduction_prefix, reduction_suffix
    ))
    return(obj)
  }

  obj_cells <- rownames(obj@meta.data)
  emb       <- nmf_result$NMF_Matrix
  missing   <- setdiff(obj_cells, rownames(emb))

  if (length(missing) > 0L) {
    warning(sprintf(
      "%d cells in obj have no NMF embedding and will have NA scores: %s%s",
      length(missing),
      paste(head(missing, 5L), collapse = ", "),
      if (length(missing) > 5L) " ..." else ""
    ))
    na_rows <- matrix(NA_real_, nrow = length(missing), ncol = ncol(emb),
                      dimnames = list(missing, colnames(emb)))
    emb <- rbind(emb, na_rows)[obj_cells, , drop = FALSE]
  } else {
    emb <- emb[intersect(obj_cells, rownames(emb)), , drop = FALSE]
  }

  reduction_name <- paste0(reduction_prefix, reduction_suffix)
  safe_suffix    <- gsub("[^A-Za-z0-9]", "", reduction_suffix)
  reduction_key  <- paste0(reduction_key_prefix, safe_suffix, "_")

  dimreduc <- Seurat::CreateDimReducObject(
    embeddings = emb,
    loadings   = if (!is.null(nmf_result$NMF_Loading)) nmf_result$NMF_Loading
    else new("matrix"),
    key        = reduction_key,
    assay      = Seurat::DefaultAssay(obj)
  )

  actual_key <- Seurat::Key(dimreduc)
  if (!identical(actual_key, reduction_key)) {
    message(sprintf(
      "Note: Seurat adjusted the reduction key for '%s' from '%s' to '%s'. ",
      reduction_name, reduction_key, actual_key
    ), sprintf(
      "Use Seurat::Key(obj[[\"%s\"]]) to confirm feature names for plotting.",
      reduction_name
    ))
  }

  dimreduc@misc <- list(
    Fit_Error        = nmf_result$Fit_Error,
    Fit_Summary      = nmf_result$Fit_Summary,
    Factor_Gene_List = nmf_result$Factor_Gene_List
  )

  obj[[reduction_name]] <- dimreduc
  obj
}

#' Retrieve NMF results from a Seurat object
#'
#' Extracts the \code{@@misc} slot of one or all NMF \code{DimReduc} objects
#' written by \code{\link{store_nmf_results}}.
#'
#' @param obj A Seurat v5 object.
#' @param reduction_suffix Character or NULL. Suffix used when storing the
#'   reduction (e.g. \code{"Keratinocytes"} for \code{"nmf_Keratinocytes"}).
#'   When NULL, all reductions with \code{reduction_prefix} are returned as a
#'   named list. If the requested suffix is not found, available NMF reductions
#'   are shown.
#' @param reduction_prefix Character. Prefix used when storing reductions.
#'   Default \code{"nmf_"}.
#'
#' @return When \code{reduction_suffix} is specified: the \code{@@misc} list
#'   for that reduction, containing \code{Fit_Error}, \code{Fit_Summary}, and
#'   \code{Factor_Gene_List}.
#'
#'   When \code{reduction_suffix = NULL}: a named list of such \code{@@misc}
#'   lists, one per matching NMF reduction, keyed by suffix.
#'
#' @seealso \code{\link{store_nmf_results}}, \code{\link{compute_nmf}}
#' @export
get_nmf_results <- function(obj, reduction_suffix = NULL,
                            reduction_prefix = "nmf_") {

  all_reductions <- Seurat::Reductions(obj) %||% character(0)
  nmf_reductions <- all_reductions[startsWith(all_reductions, reduction_prefix)]

  if (length(nmf_reductions) == 0L) {
    stop(sprintf(
      "No reductions with prefix '%s' found in obj. ",
      reduction_prefix
    ), "Available reductions: ",
    if (length(all_reductions) == 0L) "(none)"
    else paste(all_reductions, collapse = ", "), ".")
  }

  if (is.null(reduction_suffix)) {
    suffixes <- sub(paste0("^", reduction_prefix), "", nmf_reductions)
    return(stats::setNames(
      lapply(nmf_reductions, function(r) obj[[r]]@misc),
      suffixes
    ))
  }

  reduction_name <- paste0(reduction_prefix, reduction_suffix)

  if (!reduction_name %in% nmf_reductions) {
    stop(sprintf(
      "NMF reduction '%s' not found in obj.\nAvailable NMF reductions: %s.",
      reduction_name,
      paste(nmf_reductions, collapse = ", ")
    ))
  }

  obj[[reduction_name]]@misc
}
