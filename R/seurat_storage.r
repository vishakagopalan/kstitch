#' Store kstitch CCA results in a Seurat object
#'
#' Takes the named list returned by \code{\link{link_shape_and_factors}} and
#' writes results into the Seurat object in two places:
#' \enumerate{
#'   \item CSP and CEP scores are assembled across groups into object-spanning
#'     \code{DimReduc} slots named \code{<reduction_name>_csp} and
#'     \code{<reduction_name>_cep}.
#'   \item All other per-group fit information is written to
#'     \code{obj@misc$kstitch[[reduction_name]][[group_level]]}.
#' }
#'
#' @param obj A Seurat v5 object.
#' @param results Named list as returned by \code{\link{link_shape_and_factors}}.
#' @param reduction_name Character. Base name used to construct the
#'   \code{DimReduc} slot names (\code{<reduction_name>_csp} and
#'   \code{<reduction_name>_cep}) and the \code{obj@misc$kstitch} key.
#'   Must be unique across calls if storing multiple CCA results.
#'   Default \code{"kstitch"}.
#' @param reduction_key_csp Character. Key prefix for the CSP \code{DimReduc}.
#'   Default \code{"CSP_"}.
#' @param reduction_key_cep Character. Key prefix for the CEP \code{DimReduc}.
#'   Default \code{"CEP_"}.
#'
#' @return The Seurat object with \code{obj[["<reduction_name>_csp"]]},
#'   \code{obj[["<reduction_name>_cep"]]}, and
#'   \code{obj@misc$kstitch[[reduction_name]]} populated.
#'
#' @seealso \code{\link{link_shape_and_factors}}, \code{\link{get_kstitch_results}}
#' @export
store_kstitch_results <- function(obj,
                                  results,
                                  reduction_name    = "kstitch",
                                  reduction_key_csp = "CSP_",
                                  reduction_key_cep = "CEP_") {
  all_cells <- rownames(obj@meta.data)
  k         <- ncol(results[[1]]$CSP_Scores)
  csp_mat <- matrix(NA_real_, nrow = length(all_cells), ncol = k,
                    dimnames = list(all_cells, paste0("CSP", seq_len(k))))
  cep_mat <- matrix(NA_real_, nrow = length(all_cells), ncol = k,
                    dimnames = list(all_cells, paste0("CEP", seq_len(k))))
  for (grp in names(results)) {
    grp_cells <- rownames(results[[grp]]$CSP_Scores)
    csp_mat[grp_cells, ] <- results[[grp]]$CSP_Scores
    cep_mat[grp_cells, ] <- results[[grp]]$CEP_Scores
  }
  if (any(is.na(csp_mat))) {
    warning("Some cells in the Seurat object have no CSP scores (NA). ",
            "These cells were not present in any group in `results`.")
  }
  csp_name <- paste0(reduction_name, "_csp")
  cep_name <- paste0(reduction_name, "_cep")
  obj[[csp_name]] <- Seurat::CreateDimReducObject(
    embeddings = csp_mat,
    key        = reduction_key_csp,
    assay      = Seurat::DefaultAssay(obj)
  )
  obj[[cep_name]] <- Seurat::CreateDimReducObject(
    embeddings = cep_mat,
    key        = reduction_key_cep,
    assay      = Seurat::DefaultAssay(obj)
  )
  if (is.null(obj@misc$kstitch)) obj@misc$kstitch <- list()
  obj@misc$kstitch[[reduction_name]] <- results
  obj
}

#' Retrieve kstitch per-group fit results from a Seurat object
#'
#' Convenience accessor for \code{obj@misc$kstitch[[reduction_name]][[group]]}.
#'
#' @param obj A Seurat v5 object with kstitch results stored via
#'   \code{\link{store_kstitch_results}}.
#' @param reduction_name Character. The reduction name passed to
#'   \code{\link{store_kstitch_results}}. Default \code{"kstitch"}.
#' @param group Character. The group level to retrieve. If \code{NULL}
#'   (default), returns the full \code{obj@misc$kstitch[[reduction_name]]} list.
#'
#' @return A list of fit objects for the requested group, or all groups for the
#'   requested reduction if \code{group = NULL}.
#'
#' @seealso \code{\link{store_kstitch_results}}
#' @export
get_kstitch_results <- function(obj, reduction_name = "kstitch", group = NULL) {
  if (is.null(obj@misc$kstitch)) {
    stop("No kstitch results found in obj@misc. Run store_kstitch_results() first.")
  }
  if (!reduction_name %in% names(obj@misc$kstitch)) {
    stop(sprintf(
      "Reduction '%s' not found in obj@misc$kstitch. Available reductions: %s.",
      reduction_name,
      paste(names(obj@misc$kstitch), collapse = ", ")
    ))
  }
  reduction <- obj@misc$kstitch[[reduction_name]]
  if (is.null(group)) return(reduction)
  if (!group %in% names(reduction)) {
    stop(sprintf(
      "Group '%s' not found in obj@misc$kstitch[[\"%s\"]]. Available groups: %s.",
      group,
      reduction_name,
      paste(names(reduction), collapse = ", ")
    ))
  }
  reduction[[group]]
}
