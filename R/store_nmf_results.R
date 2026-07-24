#' Store per-group NMF results in a Seurat object
#'
#' Takes the named list returned by \code{\link{compute_nmf}} (one entry per
#' group, or a single \code{"all"} entry when \code{group.by = NULL}) and
#' writes each group's fit into the Seurat object as its own \code{DimReduc}
#' slot, named \code{obj[["nmf_<group>"]]}. Feature loadings are stored in the
#' reduction's \code{@feature.loadings}; fit metadata (\code{Fit_Error},
#' \code{Fit_Summary}, \code{Factor_Gene_List}) is stored in the reduction's
#' \code{@misc}. Mirrors \code{\link{store_tpca_results}}'s pattern for the
#' TPCA bridge.
#'
#' Cells not present in a given group's \code{NMF_Matrix} receive \code{NA}
#' embeddings in that group's reduction (a cell only ever belongs to one
#' group's reduction under normal use, so this mainly guards against partial
#' or mismatched inputs).
#'
#' Groups whose fit was empty (e.g. skipped for having too few cells --
#' \code{NMF_Matrix} is \code{NULL}) are skipped with a message; no
#' \code{DimReduc} is created for them.
#'
#' @param obj A Seurat v5 object.
#' @param nmf_results Named list as returned by \code{\link{compute_nmf}}.
#' @param reduction_prefix Character. Prefix for each group's reduction name.
#'   Default \code{"nmf_"}, giving \code{obj[["nmf_<group>"]]}.
#' @param reduction_key_prefix Character. Key prefix passed to
#'   \code{Seurat::CreateDimReducObject()} for each reduction (Seurat requires
#'   reduction keys to be unique across the object, and to consist of
#'   alphanumeric characters with a single trailing underscore -- Seurat
#'   silently strips any other underscores). Default \code{"NMF"}; the group
#'   name is appended (with non-alphanumeric characters removed) to keep keys
#'   unique across groups, e.g. \code{"NMFTypeA_"}, \code{"NMFTypeB_"}, giving
#'   feature names \code{"NMFTypeA_1"}, \code{"NMFTypeA_2"}, ... Use
#'   \code{Seurat::Key(obj[["nmf_<group>"]])} to confirm the exact key
#'   actually assigned, since Seurat may adjust it further.
#'
#' @return The Seurat object with one \code{obj[["nmf_<group>"]]} reduction
#'   per non-empty group.
#' @seealso \code{\link{compute_nmf}}, \code{\link{store_tpca_results}}
#' @export
store_nmf_results <- function(obj,
                              nmf_results,
                              reduction_prefix     = "nmf_",
                              reduction_key_prefix = "NMF") {

  obj_cells <- rownames(obj@meta.data)

  for (grp in names(nmf_results)) {

    fit <- nmf_results[[grp]]

    if (is.null(fit$NMF_Matrix)) {
      message(sprintf("Skipping group '%s': no NMF fit available (NMF_Matrix is NULL).", grp))
      next
    }

    emb <- fit$NMF_Matrix
    missing <- setdiff(obj_cells, rownames(emb))

    if (length(missing) > 0) {
      warning(sprintf(
        "%d cells in obj have no NMF embedding for group '%s' and will have NA scores: %s%s",
        length(missing), grp,
        paste(head(missing, 5), collapse = ", "),
        if (length(missing) > 5) " ..." else ""
      ))
      na_rows <- matrix(NA_real_, nrow = length(missing), ncol = ncol(emb),
                        dimnames = list(missing, colnames(emb)))
      emb <- rbind(emb, na_rows)[obj_cells, , drop = FALSE]
    } else {
      emb <- emb[intersect(obj_cells, rownames(emb)), , drop = FALSE]
    }

    reduction_name <- paste0(reduction_prefix, grp)

    # Seurat keys must be alphanumeric with a single trailing underscore --
    # strip anything else out of the group name before appending, and warn
    # if the requested key had to be altered so the caller isn't surprised
    # by a mismatch between what they passed and what Seurat actually used.
    safe_grp <- gsub("[^A-Za-z0-9]", "", grp)
    reduction_key <- paste0(reduction_key_prefix, safe_grp, "_")

    dimreduc <- Seurat::CreateDimReducObject(
      embeddings      = emb,
      loadings        = if (!is.null(fit$NMF_Loading)) fit$NMF_Loading else new("matrix"),
      key             = reduction_key,
      assay           = Seurat::DefaultAssay(obj)
    )

    actual_key <- Seurat::Key(dimreduc)
    if (!identical(actual_key, reduction_key)) {
      message(sprintf(
        "Note: Seurat adjusted the reduction key for group '%s' from '%s' to '%s'. ",
        grp, reduction_key, actual_key
      ), "Use Seurat::Key(obj[[\"", reduction_name, "\"]]) to confirm feature names for plotting.")
    }

    dimreduc@misc <- list(
      Fit_Error        = fit$Fit_Error,
      Fit_Summary      = fit$Fit_Summary,
      Factor_Gene_List = fit$Factor_Gene_List
    )

    obj[[reduction_name]] <- dimreduc
  }

  obj
}

#' Retrieve NMF results from a Seurat object
#'
#' Extracts the \code{@@misc} slot of one or all NMF DimReduc objects written
#' by \code{\link{store_nmf_results}}.
#'
#' @param obj A Seurat v5 object.
#' @param group Character or NULL. Name of the group to retrieve (e.g.
#'   \code{"keratinocyte"} or \code{"all"} for single-group results). When
#'   NULL, all reductions with the given prefix are returned as a named list.
#' @param reduction_prefix Character. Prefix used when storing reductions
#'   (default \code{"nmf_"}). Must match the prefix passed to
#'   \code{store_nmf_results()}.
#'
#' @return When \code{group} is specified: the \code{@@misc} list for that
#'   group's NMF reduction, containing \code{Fit_Error}, \code{Fit_Summary},
#'   and \code{Factor_Gene_List}.
#'
#'   When \code{group = NULL}: a named list of such \code{@@misc} lists, one
#'   per NMF reduction found in \code{obj}, keyed by group name.
#'
#' @seealso \code{\link{store_nmf_results}}, \code{\link{compute_nmf}}
#' @export
get_nmf_results <- function(obj, group = NULL, reduction_prefix = "nmf_") {
  if (is.null(group)) {
    all_reductions <- Seurat::Reductions(obj) %||% character(0)
    nmf_reductions <- all_reductions[startsWith(all_reductions, reduction_prefix)]
    if (length(nmf_reductions) == 0) {
      stop(sprintf(
        "No reductions with prefix '%s' found in obj. Run store_nmf_results() first.",
        reduction_prefix
      ))
    }
    grp_names <- sub(paste0("^", reduction_prefix), "", nmf_reductions)
    out <- stats::setNames(
      lapply(nmf_reductions, function(r) obj[[r]]@misc),
      grp_names
    )
    return(out)
  }

  reduction_name <- paste0(reduction_prefix, group)
  if (!reduction_name %in% Seurat::Reductions(obj)) {
    stop(sprintf(
      "Reduction '%s' not found in obj. Available reductions: %s.",
      reduction_name,
      paste(Seurat::Reductions(obj), collapse = ", ")
    ))
  }
  obj[[reduction_name]]@misc
}

