#' Link shape and gene-expression embeddings via CCA
#'
#' Runs canonical correlation analysis (CCA) between a morphology embedding
#' (e.g. TPCA shape PCs) and a gene-expression dimensionality reduction
#' (e.g. NMF factors, PCA scores, or any cell x feature matrix), optionally
#' stratified by a grouping variable. Results are returned as a named list
#' intended to be stored in \code{obj@misc$kstitch} via
#' \code{\link{store_kstitch_results}}.
#'
#' @param obj A Seurat v5 object. Used to resolve \code{group.by} via
#'   \code{@meta.data}; cell names are taken from \code{colnames(obj)}.
#' @param expr_mat Numeric matrix of gene-expression embeddings (e.g. NMF
#'   factors, PCA scores), cells x features, with cell names as row names.
#'   Any dimensionality reduction method is supported; NMF-specific
#'   post-processing (e.g. \code{\link{prioritize_cca_components}}) should
#'   be applied separately after this function returns.
#' @param shape_mat Numeric matrix of shape embeddings (e.g. TPCA scores),
#'   cells x features, with cell names as row names.
#' @param group.by Character scalar naming a column in \code{obj@meta.data} to
#'   stratify by. CCA is run once per unique value. \code{NULL} (default) pools
#'   all cells into a single group labelled \code{"all"}.
#' @param min_cells Integer. Groups with fewer cells (after intersection) are
#'   skipped with a message. Default 100.
#' @param scale Logical. Whether to z-score both matrices before CCA.
#'   Default \code{TRUE}.
#' @param test_significance Logical. Whether to compute deflation-based
#'   permutation p-values via \code{\link{cca_pvalues}} for each group.
#'   Default \code{FALSE}.
#' @param nperm Integer. Number of permutations passed to \code{cca_pvalues}
#'   when \code{test_significance = TRUE}. Default 1000L.
#' @param perm_seed Optional integer seed for reproducibility of permutations.
#' @param verbose Logical. Print progress messages. Default \code{FALSE}.
#'
#' @return A named list, one element per group, each containing:
#'   \describe{
#'     \item{CC_Corr_Coefs}{Canonical correlations, length k.}
#'     \item{CSP_Scores}{Cell x k matrix of shape canonical variates (CSP).}
#'     \item{CEP_Scores}{Cell x k matrix of expression canonical variates (CEP).}
#'     \item{CSP_Vectors}{Feature x k matrix of shape canonical weight vectors.}
#'     \item{CEP_Vectors}{Feature x k matrix of expression canonical weight vectors.}
#'     \item{CSP_Self_Correlations}{Structure correlations of shape features with CSP scores.}
#'     \item{CEP_Self_Correlations}{Structure correlations of expression features with CEP scores.}
#'     \item{Anchor_Features}{Named character vector recording the anchor feature per component.}
#'     \item{Misc_CCA}{Full \code{run_cca()} output for downstream use.}
#'   }
#'
#' @seealso \code{\link{run_cca}}, \code{\link{cca_pvalues}},
#'   \code{\link{anchor_cca_signs}}, \code{\link{prioritize_cca_components}}
#' @export
link_shape_and_factors <- function(obj,
                                   expr_mat,
                                   shape_mat,
                                   group.by          = NULL,
                                   min_cells         = 100,
                                   scale             = TRUE,
                                   test_significance = FALSE,
                                   nperm             = 1000L,
                                   perm_seed         = NULL,
                                   verbose           = FALSE) {

  # --- input checks -----------------------------------------------------------

  stopifnot(
    is.matrix(expr_mat) || is.data.frame(expr_mat),
    is.matrix(shape_mat) || is.data.frame(shape_mat),
    !is.null(rownames(expr_mat)),
    !is.null(rownames(shape_mat))
  )

  if (!is.null(group.by)) {
    if (!group.by %in% colnames(obj@meta.data)) {
      stop(sprintf(
        "`group.by` column '%s' not found in obj@meta.data. Available columns: %s.",
        group.by,
        paste(colnames(obj@meta.data), collapse = ", ")
      ))
    }
  }

  # --- resolve cells and groups -----------------------------------------------

  cells_all <- Reduce(intersect, list(
    rownames(obj@meta.data),
    rownames(expr_mat),
    rownames(shape_mat)
  ))

  if (length(cells_all) == 0) {
    stop("No cells remain after intersecting rownames(obj@meta.data), rownames(expr_mat), ",
         "and rownames(shape_mat). Check that row names match Seurat cell names.")
  }

  if (is.null(group.by)) {
    groups <- list(all = cells_all)
  } else {
    meta_vec <- obj@meta.data[cells_all, group.by]
    groups   <- split(cells_all, meta_vec)
  }

  # --- per-group CCA ----------------------------------------------------------

  results <- list()

  for (grp in names(groups)) {

    cells <- groups[[grp]]
    n     <- length(cells)

    if (n < min_cells) {
      message(sprintf(
        "Skipping group '%s': %d cells available, minimum is %d.",
        grp, n, min_cells
      ))
      next
    }

    if (verbose) message(sprintf("Running CCA for group '%s' (%d cells).", grp, n))

    X <- shape_mat[cells, , drop = FALSE]
    Y <- expr_mat[cells,  , drop = FALSE]

    if (scale) {
      X <- scale(X)
      Y <- scale(Y)
    }

    cca <- run_cca(X, Y, scale = FALSE)

    k <- length(cca$cor)

    csp_scores <- cca$scores$xscores
    colnames(csp_scores) <- paste0("CSP", seq_len(k))
    rownames(csp_scores) <- cells

    csp_vectors <- cca$xcoef
    colnames(csp_vectors) <- paste0("CSP", seq_len(k))

    cep_scores <- cca$scores$yscores
    colnames(cep_scores) <- paste0("CEP", seq_len(k))
    rownames(cep_scores) <- cells

    cep_vectors <- cca$ycoef
    colnames(cep_vectors) <- paste0("CEP", seq_len(k))

    csp_self_cor <- cca$scores$corr.X.xscores
    cep_self_cor <- cca$scores$corr.Y.yscores

    results[[grp]] <- list(
      CC_Corr_Coefs         = cca$cor,
      CSP_Scores            = csp_scores,
      CEP_Scores            = cep_scores,
      CSP_Vectors           = csp_vectors,
      CEP_Vectors           = cep_vectors,
      CSP_Self_Correlations = csp_self_cor,
      CEP_Self_Correlations = cep_self_cor,
      Misc_CCA              = cca
    )

    if (test_significance) {
      if (verbose) message(sprintf(
        "  Computing permutation p-values for group '%s' (%d perms) ...", grp, nperm
      ))
      results[[grp]][["P_Value_Info"]] <- cca_pvalues(
        X       = X,
        Y       = Y,
        nperm   = nperm,
        seed    = perm_seed,
        verbose = verbose
      )
    }

    results[[grp]] <- anchor_cca_signs(results[[grp]])
  }

  if (length(results) == 0) {
    warning("All groups were skipped (all had fewer than `min_cells` cells). ",
            "Returning an empty list.")
  }

  results
}
