#' Link cell shape and gene expression via Canonical Correlation Analysis
#'
#' Runs CCA between a shape PC matrix and an expression factor matrix.
#' Operates on a single pre-filtered cell set. For multi-group workflows,
#' call this function once per group in an external loop, passing pre-subsetted
#' matrices each time.
#'
#' @param obj A Seurat v5 object. Used only to validate that matrix row names
#'   match known cell names; never mutated.
#' @param expr_mat Numeric matrix (cells x features). Expression-side input,
#'   e.g. NMF factor scores. Row names must match Seurat cell names.
#' @param shape_mat Numeric matrix (cells x shape PCs). Shape-side input from
#'   TPCA. Row names must match Seurat cell names.
#' @param min_cells Integer. If fewer cells remain after intersecting the three
#'   cell sets, the function errors (default 100).
#' @param scale Logical. Whether to z-score \code{shape_mat} and
#'   \code{expr_mat} before CCA (default TRUE).
#' @param test_significance Logical. Whether to compute deflation-based
#'   permutation p-values via \code{cca_pvalues()} (default FALSE).
#' @param nperm Integer. Number of permutations for significance testing.
#' @param perm_seed Integer or NULL. RNG seed for permutation testing.
#' @param verbose Logical. Whether to emit progress messages (default FALSE).
#' @param return_results Logical. When TRUE (default) the result list is
#'   returned. When FALSE, the result is serialized to
#'   \code{<output_dir>/all.rds} and the file path is returned.
#' @param output_dir Character or NULL. Directory for serialized output when
#'   \code{return_results = FALSE}. When NULL a temporary directory is used
#'   and not cleaned up. Ignored when \code{return_results = TRUE}.
#'
#' @return When \code{return_results = TRUE}: a list containing
#'   \code{CC_Corr_Coefs}, \code{Shape_Corr_With_CSP}, \code{Exp_Corr_With_CEP},
#'   \code{CSP_Scores}, \code{CEP_Scores}, \code{CSP_Vectors},
#'   \code{CEP_Vectors}, \code{Misc_CCA}, and \code{Anchor_Features}. If
#'   \code{test_significance = TRUE}, also \code{P_Value_Info}.
#'
#'   When \code{return_results = FALSE}: a list with element
#'   \code{output_path} (character scalar, path to the serialized RDS).
#'
#' @seealso \code{\link{store_kstitch_results}}, \code{\link{load_kstitch_results}}
#' @export
link_shape_and_factors <- function(obj,
                                   expr_mat,
                                   shape_mat,
                                   min_cells         = 100,
                                   scale             = TRUE,
                                   test_significance = FALSE,
                                   nperm             = 1000L,
                                   perm_seed         = NULL,
                                   verbose           = FALSE,
                                   return_results    = TRUE,
                                   output_dir        = NULL) {

  stopifnot(
    is.matrix(expr_mat)  || is.data.frame(expr_mat),
    is.matrix(shape_mat) || is.data.frame(shape_mat),
    !is.null(rownames(expr_mat)),
    !is.null(rownames(shape_mat))
  )

  if (!return_results && is.null(output_dir)) {
    output_dir <- file.path(tempdir(), paste0("kstitch_cca_", .random_id()))
    dir.create(output_dir, recursive = TRUE)
  } else if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  }

  cells <- Reduce(intersect, list(
    rownames(obj@meta.data),
    rownames(expr_mat),
    rownames(shape_mat)
  ))

  if (length(cells) == 0L) {
    stop(
      "No cells remain after intersecting rownames(obj@meta.data), ",
      "rownames(expr_mat), and rownames(shape_mat). ",
      "Check that row names match Seurat cell names."
    )
  }

  if (length(cells) < min_cells) {
    stop(sprintf(
      "%d cells available after intersection, minimum is %d.",
      length(cells), min_cells
    ))
  }

  X <- shape_mat[cells, , drop = FALSE]
  Y <- expr_mat[cells,  , drop = FALSE]

  complete_idx <- complete.cases(X) & complete.cases(Y)
  n_dropped    <- sum(!complete_idx)
  if (n_dropped > 0L) {
    message(sprintf(
      "Dropping %d cells with NA in shape_mat or expr_mat (%d remain).",
      n_dropped, sum(complete_idx)
    ))
    X <- X[complete_idx, , drop = FALSE]
    Y <- Y[complete_idx, , drop = FALSE]
  }

  if (nrow(X) < min_cells) {
    stop(sprintf(
      "Only %d complete cells remain after NA removal, minimum is %d.",
      nrow(X), min_cells
    ))
  }

  cells <- rownames(X)

  if (scale) {
    X <- scale(X)
    Y <- scale(Y)
  }

  if (verbose) message(sprintf("Running CCA on %d cells ...", length(cells)))

  cca <- run_cca(X, Y, scale = FALSE)
  k   <- length(cca$cor)

  csp_scores           <- cca$scores$xscores
  colnames(csp_scores) <- paste0("CSP", seq_len(k))
  rownames(csp_scores) <- cells

  csp_vectors           <- cca$xcoef
  colnames(csp_vectors) <- paste0("CSP", seq_len(k))

  cep_scores           <- cca$scores$yscores
  colnames(cep_scores) <- paste0("CEP", seq_len(k))
  rownames(cep_scores) <- cells

  cep_vectors           <- cca$ycoef
  colnames(cep_vectors) <- paste0("CEP", seq_len(k))

  csp_self_cor <- cca$scores$corr.X.xscores
  cep_self_cor <- cca$scores$corr.Y.yscores

  names(cca$scores)[match(
    c("corr.X.xscores", "corr.X.yscores", "corr.Y.xscores", "corr.Y.yscores"),
    names(cca$scores)
  )] <- c("corr.shape.with.csp", "corr.shape.with.cep",
          "corr.exp.with.csp",   "corr.exp.with.cep")

  result <- list(
    CC_Corr_Coefs       = cca$cor,
    Shape_Corr_With_CSP = csp_self_cor,
    Exp_Corr_With_CEP   = cep_self_cor,
    CSP_Scores          = csp_scores,
    CEP_Scores          = cep_scores,
    CSP_Vectors         = csp_vectors,
    CEP_Vectors         = cep_vectors,
    Misc_CCA            = cca
  )

  if (test_significance) {
    if (verbose) message(sprintf("Computing permutation p-values (%d perms) ...", nperm))
    result[["P_Value_Info"]] <- cca_pvalues(
      X       = X,
      Y       = Y,
      nperm   = nperm,
      seed    = perm_seed,
      verbose = verbose
    )
  }

  result <- anchor_cca_signs(result)

  if (!return_results) {
    path <- file.path(output_dir, "all.rds")
    saveRDS(result, file = path)
    message("CCA result written to: ", path)
    message(sprintf("Load with load_kstitch_results(\"%s\", type = \"cca\")", path))
    return(list(output_path = path))
  }

  result
}
