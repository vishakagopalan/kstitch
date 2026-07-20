#' Prioritize CCA components by cross-loading stability
#'
#' Identifies the most stable features on each canonical component via
#' repeated subsampling. On each replicate, \code{\link{run_cca}} is run on a
#' random fraction of cells and features are ranked by their absolute
#' cross-loading magnitude. Features that consistently appear among the top
#' ranks across replicates are retained.
#'
#' @param shape_mat Numeric matrix, cells x shape features, rownames = cell
#'   names. Should be the same (optionally scaled) matrix passed to
#'   \code{\link{link_shape_and_factors}}.
#' @param expr_mat Numeric matrix, cells x NMF factors, rownames = cell names.
#' @param ccs_to_consider Integer vector of component indices to evaluate
#'   (e.g. \code{1:3}).
#' @param shape_side_rank Integer. A shape feature must rank within the top
#'   \code{shape_side_rank} positions to count toward its tally. Default 5.
#' @param factor_side_rank Integer. Same threshold for NMF factors. Default 5.
#' @param num_replicates Integer. Number of subsampling replicates. Default 100.
#' @param subsample_frac Numeric in (0, 1]. Fraction of cells to use per
#'   replicate. Default 0.9.
#' @param cell_subsample Optional character vector of cell names to restrict
#'   the pool before subsampling (e.g. already-intersected cells for a group).
#'   Distinct from \code{group.by} filtering — this is purely for replicate
#'   subsampling.
#' @param num_times_among_top_ranks Integer. Minimum number of replicates in
#'   which a feature must appear within the top-rank threshold to be retained.
#'   Default 75.
#'
#' @return A list with two data frames:
#'   \describe{
#'     \item{Shape}{Shape features passing the stability filter, with columns
#'       \code{feature_name}, \code{CC}, \code{cross_loading_mag}, \code{n}.}
#'     \item{Factor}{NMF factors passing the stability filter, same columns.}
#'   }
#'
#' @seealso \code{\link{link_shape_and_factors}}, \code{\link{run_cca}}
#' @export
prioritize_cca_components <- function(shape_mat,
                                      expr_mat,
                                      ccs_to_consider,
                                      shape_side_rank          = 5,
                                      factor_side_rank         = 5,
                                      num_replicates           = 100,
                                      subsample_frac           = 0.9,
                                      cell_subsample           = NULL,
                                      num_times_among_top_ranks = 75) {

  if (!is.matrix(shape_mat)) shape_mat <- as.matrix(shape_mat)
  if (!is.matrix(expr_mat))   expr_mat   <- as.matrix(expr_mat)

  num_factors        <- ncol(expr_mat)
  num_shape_features <- ncol(shape_mat)

  cell_list <- intersect(rownames(shape_mat), rownames(expr_mat))
  if (!is.null(cell_subsample)) {
    cell_list <- intersect(cell_list, cell_subsample)
  }
  if (length(cell_list) == 0) stop("No cells remain after intersecting row names.")

  combined_rank_df <- data.frame()

  for (rep_i in seq_len(num_replicates)) {

    sub_cells <- sample(cell_list, floor(subsample_frac * length(cell_list)))

    fit <- run_cca(
      shape_mat[sub_cells, , drop = FALSE],
      expr_mat[sub_cells,   , drop = FALSE],
      scale = FALSE   # caller is responsible for scaling upstream
    )

    for (cc_idx in ccs_to_consider) {

      if (cc_idx > length(fit$cor)) next

      factor_abs_cor <- abs(fit$scores$corr.Y.xscores[, cc_idx])
      shape_abs_cor  <- abs(fit$scores$corr.X.yscores[, cc_idx])

      factor_rank_df <- tibble::enframe(
        sort(factor_abs_cor, decreasing = TRUE),
        name = "feature_name", value = "loading_mag"
      ) |> dplyr::mutate(
        rank         = seq_len(num_factors),
        CC           = cc_idx,
        feature_type = "Factor",
        replicate    = rep_i
      )

      shape_rank_df <- tibble::enframe(
        sort(shape_abs_cor, decreasing = TRUE),
        name = "feature_name", value = "loading_mag"
      ) |> dplyr::mutate(
        rank         = seq_len(num_shape_features),
        CC           = cc_idx,
        feature_type = "Shape",
        replicate    = rep_i
      )

      combined_rank_df <- rbind(combined_rank_df, factor_rank_df, shape_rank_df)
    }
  }

  top_factors_df <- dplyr::filter(combined_rank_df, feature_type == "Factor") |>
    dplyr::group_by(feature_name, CC) |>
    dplyr::reframe(
      cross_loading_mag = median(loading_mag),
      n                 = sum(rank < factor_side_rank)
    ) |>
    dplyr::filter(n > num_times_among_top_ranks)

  top_shape_df <- dplyr::filter(combined_rank_df, feature_type == "Shape") |>
    dplyr::group_by(feature_name, CC) |>
    dplyr::reframe(
      cross_loading_mag = median(loading_mag),
      n                 = sum(rank < shape_side_rank)      # fix: was factor_side_rank
    ) |>
    dplyr::filter(n > num_times_among_top_ranks)

  list(Shape = top_shape_df, Factor = top_factors_df)
}
