#' Store kstitch CCA results in a Seurat object
#'
#' Takes the list returned by \code{\link{link_shape_and_factors}} and writes
#' results into the Seurat object in two places:
#' \enumerate{
#'   \item CSP and CEP scores are written into object-spanning \code{DimReduc}
#'     slots named \code{<reduction_name>_csp} and \code{<reduction_name>_cep}.
#'     Cells in the object that were not present in the CCA (e.g. cells
#'     excluded due to missing shape polygons) receive \code{NA} scores.
#'   \item All other fit information is written to
#'     \code{obj@misc$kstitch[[reduction_name]]}.
#' }
#'
#' To store multiple CCA results (e.g. cell vs nucleus, regressed vs
#' unregressed, or different cell types), call this function once per result
#' with a distinct \code{reduction_name}.
#'
#' @param obj A Seurat v5 object.
#' @param result List as returned by \code{\link{link_shape_and_factors}}.
#' @param reduction_name Character. Base name used to construct the
#'   \code{DimReduc} slot names (\code{<reduction_name>_csp} and
#'   \code{<reduction_name>_cep}) and the \code{obj@misc$kstitch} key.
#'   Must be unique across calls. Default \code{"kstitch"}.
#' @param reduction_key_csp Character. Key prefix for the CSP
#'   \code{DimReduc}. Default \code{"CSP_"}.
#' @param reduction_key_cep Character. Key prefix for the CEP
#'   \code{DimReduc}. Default \code{"CEP_"}.
#'
#' @return The Seurat object with \code{obj[["<reduction_name>_csp"]]},
#'   \code{obj[["<reduction_name>_cep"]]}, and
#'   \code{obj@misc$kstitch[[reduction_name]]} populated.
#'
#' @seealso \code{\link{link_shape_and_factors}}, \code{\link{get_kstitch_results}}
#' @export
store_kstitch_results <- function(obj,
                                  result,
                                  reduction_name    = "kstitch",
                                  reduction_key_csp = "CSP_",
                                  reduction_key_cep = "CEP_") {

  all_cells <- rownames(obj@meta.data)
  k         <- ncol(result$CSP_Scores)

  csp_mat <- matrix(NA_real_, nrow = length(all_cells), ncol = k,
                    dimnames = list(all_cells, paste0("CSP", seq_len(k))))
  cep_mat <- matrix(NA_real_, nrow = length(all_cells), ncol = k,
                    dimnames = list(all_cells, paste0("CEP", seq_len(k))))

  scored_cells <- rownames(result$CSP_Scores)
  present      <- intersect(scored_cells, all_cells)

  if (length(present) < length(scored_cells)) {
    warning(sprintf(
      "%d cells in the CCA result are not present in obj and will be ignored.",
      length(scored_cells) - length(present)
    ))
  }

  missing <- setdiff(all_cells, scored_cells)
  if (length(missing) > 0L) {
    warning(sprintf(
      "%d cells in obj have no CCA scores and will have NA scores: %s%s",
      length(missing),
      paste(head(missing, 5L), collapse = ", "),
      if (length(missing) > 5L) " ..." else ""
    ))
  }

  csp_mat[present, ] <- result$CSP_Scores[present, ]
  cep_mat[present, ] <- result$CEP_Scores[present, ]

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
  obj@misc$kstitch[[reduction_name]] <- result

  obj
}

#' Retrieve kstitch CCA results from a Seurat object
#'
#' Convenience accessor for \code{obj@misc$kstitch[[reduction_name]]}.
#'
#' @param obj A Seurat v5 object with results stored via
#'   \code{\link{store_kstitch_results}}.
#' @param reduction_name Character. The reduction name passed to
#'   \code{\link{store_kstitch_results}}. Default \code{"kstitch"}.
#'   If not found, available reduction names are shown.
#'
#' @return The CCA result list for the requested reduction.
#'
#' @seealso \code{\link{store_kstitch_results}}
#' @export
get_kstitch_results <- function(obj, reduction_name = "kstitch") {

  if (is.null(obj@misc$kstitch)) {
    stop("No kstitch results found in obj@misc. Run store_kstitch_results() first.")
  }

  available <- names(obj@misc$kstitch)

  if (!reduction_name %in% available) {
    stop(sprintf(
      "Reduction '%s' not found in obj@misc$kstitch.\nAvailable reductions: %s.",
      reduction_name,
      paste(available, collapse = ", ")
    ))
  }

  obj@misc$kstitch[[reduction_name]]
}
