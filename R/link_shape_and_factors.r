#' Link cell shape and gene expression via Canonical Correlation Analysis
#'
#' Runs CCA between a shape PC matrix and an expression factor matrix,
#' optionally per cell group. Results are returned as a named list of per-group
#' CCA outputs, or serialized to disk when \code{return_results = FALSE}.
#'
#' @param obj A Seurat v5 object. Used only for \code{@@meta.data} (group
#'   resolution and cell-name intersection); never mutated.
#' @param expr_mat Numeric matrix (cells x features). Expression-side input,
#'   e.g. NMF factor scores. Row names must match Seurat cell names.
#' @param shape_mat Numeric matrix (cells x shape PCs). Shape-side input from
#'   TPCA. Row names must match Seurat cell names.
#' @param group.by Character or NULL. Column in \code{obj@@meta.data} defining
#'   groups. When NULL, all cells are treated as a single group returned under
#'   the key \code{"all"}.
#' @param min_cells Integer. Groups with fewer cells are skipped (default 100).
#' @param scale Logical. Whether to z-score \code{shape_mat} and
#'   \code{expr_mat} before CCA (default TRUE).
#' @param test_significance Logical. Whether to compute deflation-based
#'   permutation p-values via \code{cca_pvalues()} (default FALSE).
#' @param nperm Integer. Number of permutations for significance testing.
#' @param perm_seed Integer or NULL. RNG seed for permutation testing.
#' @param verbose Logical. Whether to emit per-group progress messages
#'   (default FALSE).
#' @param return_results Logical. When TRUE (default) results are returned as a
#'   named list. When FALSE, each group's result is serialized to
#'   \code{<output_dir>/<group>.rds} and a named list of file paths is returned
#'   instead.
#' @param output_dir Character or NULL. Directory for serialized output when
#'   \code{return_results = FALSE}. When NULL a temporary directory is created
#'   automatically and is \emph{not} cleaned up on exit. Ignored when
#'   \code{return_results = TRUE}.
#'
#' @return When \code{return_results = TRUE}: a named list of per-group result
#'   lists. Each entry contains \code{CC_Corr_Coefs}, \code{Shape_Corr_With_CSP},
#'   \code{Exp_Corr_With_CEP}, \code{CSP_Scores}, \code{CEP_Scores},
#'   \code{CSP_Vectors}, \code{CEP_Vectors}, \code{Misc_CCA},
#'   \code{Anchor_Features}, \code{group}, and \code{is_groupwise}. If
#'   \code{test_significance = TRUE}, also \code{P_Value_Info}. Single-group
#'   results are returned under the key \code{"all"}.
#'
#'   When \code{return_results = FALSE}: a list with element
#'   \code{output_paths}, a named character vector mapping group names to RDS
#'   file paths.
#'
#' @seealso \code{\link{store_kstitch_results}}, \code{\link{load_kstitch_results}}
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
                                   verbose           = FALSE,
                                   return_results    = TRUE,
                                   output_dir        = NULL) {

  # ---- input checks --------------------------------------------------------
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

  # ---- resolve output directory --------------------------------------------
  if (!return_results && is.null(output_dir)) {
    output_dir <- file.path(tempdir(), paste0("kstitch_cca_", .random_id()))
    dir.create(output_dir, recursive = TRUE)
  } else if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  }

  # ---- resolve cells and groups --------------------------------------------
  cells_all <- Reduce(intersect, list(
    rownames(obj@meta.data),
    rownames(expr_mat),
    rownames(shape_mat)
  ))

  if (length(cells_all) == 0) {
    stop("No cells remain after intersecting rownames(obj@meta.data), rownames(expr_mat), ",
         "and rownames(shape_mat). Check that row names match Seurat cell names.")
  }

  is_groupwise <- !is.null(group.by)

  if (is.null(group.by)) {
    groups <- list(all = cells_all)
  } else {
    meta_vec <- obj@meta.data[cells_all, group.by]
    groups   <- split(cells_all, meta_vec)
  }

  # ---- per-group CCA -------------------------------------------------------
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

    # Drop rows where either matrix has any NA
    complete_idx <- complete.cases(X) & complete.cases(Y)
    n_dropped    <- sum(!complete_idx)
    if (n_dropped > 0L)
      message(sprintf(
        "  Group '%s': dropping %d cells with NA in shape_mat or expr_mat (%d remain).",
        grp, n_dropped, sum(complete_idx)
      ))
      X <- X[complete_idx, , drop = FALSE]
      Y <- Y[complete_idx, , drop = FALSE]

      if (nrow(X) < min_cells) {
        message(sprintf(
          "Skipping group '%s': only %d complete cells remain after NA removal, minimum is %d.",
          grp, nrow(X), min_cells
        ))
        next
      }

    cells <- rownames(X)  # update cells to complete set for score rowname assignment
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

    names(cca$scores)[match(c("corr.X.xscores", "corr.X.yscores",
                              "corr.Y.xscores", "corr.Y.yscores"),
                            names(cca$scores))] <-
      c("corr.shape.with.csp", "corr.shape.with.cep",
        "corr.exp.with.csp",   "corr.exp.with.cep")

    res <- list(
      CC_Corr_Coefs       = cca$cor,
      Shape_Corr_With_CSP = csp_self_cor,
      Exp_Corr_With_CEP   = cep_self_cor,
      CSP_Scores          = csp_scores,
      CEP_Scores          = cep_scores,
      CSP_Vectors         = csp_vectors,
      CEP_Vectors         = cep_vectors,
      Misc_CCA            = cca,
      group               = grp,
      is_groupwise        = is_groupwise
    )

    if (test_significance) {
      if (verbose) message(sprintf(
        "  Computing permutation p-values for group '%s' (%d perms) ...", grp, nperm
      ))
      res[["P_Value_Info"]] <- cca_pvalues(
        X       = X,
        Y       = Y,
        nperm   = nperm,
        seed    = perm_seed,
        verbose = verbose
      )
    }

    results[[grp]] <- anchor_cca_signs(res)
  }

  if (length(results) == 0) {
    warning("All groups were skipped (all had fewer than `min_cells` cells). ",
            "Returning an empty list.")
  }

  # ---- return_results = FALSE ----------------------------------------------
  if (!return_results) {
    output_paths <- .serialize_group_results(results, output_dir, type = "cca")
    return(list(output_paths = output_paths))
  }

  results
}
