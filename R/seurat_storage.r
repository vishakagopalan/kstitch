#' Store kstitch CCA results in a Seurat object
#'
#' Takes the named list returned by \code{\link{link_shape_and_factors}} and
#' writes results into the Seurat object in two places:
#' \enumerate{
#'   \item CSP and CEP scores are assembled across groups into object-spanning
#'     \code{DimReduc} slots (\code{obj[["csp"]]} and \code{obj[["cep"]]}).
#'   \item All other per-group fit information is written to
#'     \code{obj@misc$kstitch[[group_level]]}.
#' }
#'
#' @param obj A Seurat v5 object.
#' @param results Named list as returned by \code{\link{link_shape_and_factors}}.
#' @param reduction_key_csp Character. Key prefix for the CSP \code{DimReduc}.
#'   Default \code{"CSP_"}.
#' @param reduction_key_cep Character. Key prefix for the CEP \code{DimReduc}.
#'   Default \code{"CEP_"}.
#'
#' @return The Seurat object with \code{obj[["csp"]]}, \code{obj[["cep"]]}, and
#'   \code{obj@misc$kstitch} populated.
#'
#' @seealso \code{\link{link_shape_and_factors}}, \code{\link{get_kstitch_results}}
#' @export
store_kstitch_results <- function(obj,
                                  results,
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

  obj[["csp"]] <- Seurat::CreateDimReducObject(
    embeddings = csp_mat,
    key        = reduction_key_csp,
    assay      = Seurat::DefaultAssay(obj)
  )

  obj[["cep"]] <- Seurat::CreateDimReducObject(
    embeddings = cep_mat,
    key        = reduction_key_cep,
    assay      = Seurat::DefaultAssay(obj)
  )

  if (is.null(obj@misc$kstitch)) obj@misc$kstitch <- list()

  score_fields <- c("CSP_Scores", "CEP_Scores")

  for (grp in names(results)) {
    obj@misc$kstitch[[grp]] <- results[[grp]][
      setdiff(names(results[[grp]]), score_fields)
    ]
  }

  obj
}


#' Retrieve kstitch per-group fit results from a Seurat object
#'
#' Convenience accessor for \code{obj@misc$kstitch[[group]]}.
#'
#' @param obj A Seurat v5 object with kstitch results stored via
#'   \code{\link{store_kstitch_results}}.
#' @param group Character. The group level to retrieve. If \code{NULL}
#'   (default), returns the full \code{obj@misc$kstitch} list.
#'
#' @return A list of fit objects for the requested group, or the full kstitch
#'   misc list if \code{group = NULL}.
#'
#' @seealso \code{\link{store_kstitch_results}}
#' @export
get_kstitch_results <- function(obj, group = NULL) {

  if (is.null(obj@misc$kstitch)) {
    stop("No kstitch results found in obj@misc. Run store_kstitch_results() first.")
  }

  if (is.null(group)) return(obj@misc$kstitch)

  if (!group %in% names(obj@misc$kstitch)) {
    stop(sprintf(
      "Group '%s' not found in obj@misc$kstitch. Available groups: %s.",
      group,
      paste(names(obj@misc$kstitch), collapse = ", ")
    ))
  }

  obj@misc$kstitch[[group]]
}
