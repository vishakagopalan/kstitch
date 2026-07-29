# ---- Figure 1: TPCA shape-mode grid ----------------------------------------

#' Plot TPCA shape modes
#'
#' Reconstructs and plots cell/nucleus shapes along each principal geodesic
#' component at a range of standard deviation values, producing a PC x SD
#' faceted grid.
#'
#' Calls \code{reconstruct_shapes_from_pca()} from \code{kendall_tpca.py} via
#' reticulate.
#'
#' @param tpca_result A single group's TPCA result — i.e. one element of the
#'   named list returned by \code{\link{run_tpca}} or
#'   \code{\link{run_tpca_from_seurat}} (e.g. \code{tpca_result$all} or
#'   \code{tpca_result[["keratinocyte"]]}). Do not pass the top-level named
#'   list directly.
#' @param sds_to_plot Numeric vector of SD values to plot along each PC.
#'   Default \code{seq(-3, 3, by = 1)}.
#' @param num_pcs Integer. Number of PCs to show. Default 4L.
#' @param line_colour Character. Colour for shape outlines. Default
#'   \code{"steelblue"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_shape_modes <- function(tpca_result,
                             sds_to_plot  = seq(-3, 3, by = 1),
                             num_pcs      = 4L,
                             line_colour  = "steelblue") {

  info <- tpca_result$Info

  kendall_tpca_py_dir <- system.file("python", package = "kstitch")
  kendall_tpca        <- reticulate::import_from_path("kendall_tpca", path = kendall_tpca_py_dir)

  raw <- kendall_tpca$reconstruct_shapes_from_pca(
    v           = t(info$v_matrix),
    mu          = t(info$frechet_mean),
    lambdas     = info$variances,
    sds_to_plot = sds_to_plot,
    num_pcs     = as.integer(num_pcs)
  )
  df <- as.data.frame(lapply(raw, as.vector))

  variances  <- info$variances
  pct_var    <- round(100 * variances / sum(variances), 2)
  pct_var_df <- data.frame(
    PC      = paste0("PC", seq_len(length(pct_var))),
    disp_pc = paste0("PC", seq_len(length(pct_var)),
                     "\n(", pct_var, "%)")
  )

  df <- merge(df, pct_var_df, by = "PC")
  df$disp_pc <- factor(df$disp_pc,
                       levels = pct_var_df$disp_pc[seq_len(num_pcs)])
  df$sd_val  <- factor(df$sd_val)

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, group = interaction(PC, sd_val))) +
    ggplot2::geom_polygon(fill = NA, colour = line_colour, linewidth = 0.6) +
    ggplot2::facet_grid(disp_pc ~ sd_val, switch = "y") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title    = "Shape modes (TPCA)",
      subtitle = paste0("Contour type: ", tpca_result$contour_type),
      x        = expression(sigma),
      y        = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 9),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}


# ---- Figure 2: CEP loading bar plot -----------------------------------------

#' Plot CEP (expression canonical variate) loadings
#'
#' Bar chart of NMF factor loadings on a given CEP component.
#'
#' @param group_result A single group's entry from the list returned by
#'   \code{\link{link_shape_and_factors}}.
#' @param cc_idx Integer. Which canonical component to plot. Default 1L.
#' @param bar_colour_pos Character. Fill colour for positive loadings.
#'   Default \code{"steelblue"}.
#' @param bar_colour_neg Character. Fill colour for negative loadings.
#'   Default \code{"tomato"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_cep_loadings <- function(group_result,
                              cc_idx         = 1L,
                              bar_colour_pos = "steelblue",
                              bar_colour_neg = "tomato") {

  vec <- group_result$CEP_Vectors[, cc_idx]

  df <- data.frame(
    feature = names(vec),
    loading = as.numeric(vec)
  )
  df$feature <- factor(df$feature, levels = df$feature[order(df$loading)])
  df$sign    <- ifelse(df$loading >= 0, "positive", "negative")

  ggplot2::ggplot(df, ggplot2::aes(x = feature, y = loading, fill = sign)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c(positive = bar_colour_pos,
                                          negative = bar_colour_neg)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Loading on CEP", cc_idx),
      x     = NULL,
      y     = "Canonical weight"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}


# ---- Figure 3: CSP loading bar plot -----------------------------------------

#' Plot CSP (shape canonical variate) loadings
#'
#' Bar chart of shape feature loadings on a given CSP component.
#'
#' @param group_result A single group's entry from the list returned by
#'   \code{\link{link_shape_and_factors}}.
#' @param cc_idx Integer. Which canonical component to plot. Default 1L.
#' @param bar_colour_pos Character. Fill colour for positive loadings.
#'   Default \code{"steelblue"}.
#' @param bar_colour_neg Character. Fill colour for negative loadings.
#'   Default \code{"tomato"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_csp_loadings <- function(group_result,
                              cc_idx         = 1L,
                              bar_colour_pos = "steelblue",
                              bar_colour_neg = "tomato") {

  vec <- group_result$CSP_Vectors[, cc_idx]

  df <- data.frame(
    feature = names(vec),
    loading = as.numeric(vec)
  )
  df$feature <- factor(df$feature, levels = df$feature[order(df$loading)])
  df$sign    <- ifelse(df$loading >= 0, "positive", "negative")

  anchor <- group_result$Anchor_Features[[paste0("CSP", cc_idx)]]

  ggplot2::ggplot(df, ggplot2::aes(x = feature, y = loading, fill = sign)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c(positive = bar_colour_pos,
                                          negative = bar_colour_neg)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = paste0("Loading on CSP", cc_idx),
      subtitle = paste0("Anchor feature: ", anchor),
      x        = NULL,
      y        = "Canonical weight"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}


# ---- Figure 4: CSP-binned boundary montage ----------------------------------

#' Boundary montage along a CSP or CEP axis
#'
#' Bins cells by their score on a chosen canonical variate (CSP or CEP),
#' selects medoid cells within each bin, and plots their Procrustes-aligned
#' contours in a faceted montage. When the pre-shape embedding is unavailable,
#' shapes are reconstructed from TPCA PCA coordinates.
#'
#' @param group_result A kstitch group result object containing \code{CSP_Scores}
#'   and \code{CEP_Scores} matrices (cells x canonical variates).
#' @param tpca_result A kstitch TPCA result object containing \code{Info},
#'   \code{TPCA_Embedding}, and optionally \code{Metadata} and
#'   \code{contour_type}.
#' @param cc_idx Integer. Index of the canonical variate to bin cells by.
#'   Default \code{1L}.
#' @param num_bins Integer. Number of quantile bins along the score axis.
#'   Default \code{10L}.
#' @param num_cells_for_boundary Integer. Number of medoid cells to display per
#'   bin. Default \code{6L}.
#' @param pcs_to_use Integer vector. TPCA PCs used for medoid selection and,
#'   when reconstructing shapes, for shape reconstruction. PCs exceeding the
#'   embedding rank are silently dropped. Default \code{1:10}.
#' @param norm_rescale Logical. If \code{TRUE}, rescales each reconstructed or
#'   retrieved contour by the per-cell Frobenius norm stored in
#'   \code{tpca_result$Metadata$scale}, restoring approximate original contour
#'   size. Disabled with a message if the \code{scale} column is absent.
#'   Default \code{FALSE}.
#' @param reconstruction_pcs Integer or \code{NULL}. When shapes are
#'   reconstructed (i.e. \code{pre_shape_embedding} is \code{NULL}), caps the
#'   number of PCs from \code{pcs_to_use} passed to
#'   \code{reconstruct_shape_from_pca_coords}. \code{NULL} uses all of
#'   \code{pcs_to_use}. Default \code{NULL}.
#' @param score_type Character. Which canonical variate scores to bin cells by.
#'   One of \code{"CSP"} (cell shape program) or \code{"CEP"} (cell expression
#'   program). Default \code{"CSP"}.
#' @param line_colour Character. Colour passed to \code{geom_polygon} for
#'   contour outlines. Default \code{"steelblue"}.
#'
#' @return A \code{ggplot} object. The assembled coordinate data frame is
#'   attached as \code{attr(p, "coord_data")} with columns \code{x},
#'   \code{y}, \code{cell}, and \code{bin}.
#'
#' @details
#' Medoid cells are identified as those with the smallest mean pairwise
#' Euclidean distance to other cells in the same bin, computed in the subspace
#' defined by \code{pcs_to_use}.
#'
#' When \code{tpca_result$Info$pre_shape_embedding} is \code{NULL}, shapes are
#' reconstructed via \code{kendall_tpca$reconstruct_shape_from_pca_coords}
#' using the Fréchet mean and eigenvectors stored in \code{tpca_result$Info}.
#' A message is emitted indicating how many PCs were used; fine morphological
#' detail may be lost relative to the stored pre-shape path.
#'
#' All contours are Procrustes-aligned to the Fréchet mean via
#' \code{kendall_tpca$reparam_OPA} before plotting.
#'
#' @seealso \code{\link{plot_mu_history}}, \code{\link{plot_frechet_convergence}}
#'
#' @export
plot_csp_boundary_montage <- function(group_result,
                                      tpca_result,
                                      cc_idx                 = 1L,
                                      num_bins               = 10L,
                                      num_cells_for_boundary = 6L,
                                      pcs_to_use             = 1:10,
                                      norm_rescale           = FALSE,
                                      reconstruction_pcs     = NULL,
                                      score_type             = c("CSP", "CEP"),
                                      line_colour            = "steelblue") {

  score_type <- match.arg(score_type)

  kendall_tpca_py_dir <- system.file("python", package = "kstitch")
  kendall_tpca        <- reticulate::import_from_path("kendall_tpca", path = kendall_tpca_py_dir)

  info       <- tpca_result$Info
  mean_shape <- t(info$frechet_mean)   # 2 x L
  emb        <- if (score_type == "CSP") group_result$CSP_Scores else group_result$CEP_Scores
  shape_emb  <- tpca_result$TPCA_Embedding
  valid_idx  <- which(complete.cases(shape_emb))
  shape_emb  <- shape_emb[valid_idx, , drop = FALSE]

  pre_shape          <- info$pre_shape_embedding   # cells x 2 x L, or NULL
  use_reconstruction <- is.null(pre_shape)

  pcs_to_use <- pcs_to_use[pcs_to_use <= ncol(shape_emb)]

  if (use_reconstruction) {
    k_rec <- if (is.null(reconstruction_pcs)) length(pcs_to_use) else
      min(reconstruction_pcs, length(pcs_to_use))
    recon_msg <- sprintf(
      paste0("Pre-shape embedding not available \u2014 shapes reconstructed ",
             "from top %d PCA component%s. Fine morphological detail may be lost."),
      k_rec, if (k_rec == 1L) "" else "s"
    )
    message(recon_msg)
  } else {
    recon_msg <- NULL
    k_rec     <- NULL
  }

  # Validate norm_rescale — requires scale column in Metadata
  meta <- tpca_result$Metadata
  if (norm_rescale && (is.null(meta) || !"scale" %in% names(meta))) {
    message("scale column not found in tpca_result$Metadata; norm rescaling disabled.")
    norm_rescale <- FALSE
  }

  # bin cells by score
  cells_with_both <- intersect(rownames(emb), rownames(shape_emb))
  scores          <- emb[cells_with_both, cc_idx]

  bin_labels <- paste0("Bin ", seq_len(num_bins))
  bin_vec    <- dplyr::ntile(scores, num_bins)
  bin_df     <- data.frame(
    cell  = cells_with_both,
    score = scores,
    bin   = factor(paste0("Bin ", bin_vec), levels = bin_labels)
  )

  coord_df <- data.frame()

  for (bin_ in levels(bin_df$bin)) {

    cells_in_bin <- dplyr::filter(bin_df, bin == bin_)[["cell"]]
    if (length(cells_in_bin) < 1L) next

    # medoid selection: smallest mean pairwise TPCA distance
    sub_emb   <- shape_emb[cells_in_bin, pcs_to_use, drop = FALSE]
    dist_mat  <- as.matrix(dist(sub_emb))
    mean_dist <- colMeans(dist_mat)
    sorted    <- sort(mean_dist)
    boundary_cells <- names(sorted)[seq_len(min(num_cells_for_boundary,
                                                length(sorted)))]

    for (cell in boundary_cells) {

      if (use_reconstruction) {
        pca_coords <- shape_emb[cell, pcs_to_use[seq_len(k_rec)]]
        shape_mat  <- kendall_tpca$reconstruct_shape_from_pca_coords(
          pca_coords = reticulate::r_to_py(as.numeric(pca_coords)),
          v          = t(info$v_matrix),
          mu         = t(info$frechet_mean),
          lambdas    = info$variances,
          return_df  = FALSE
        )
        mat <- t(shape_mat)   # L x 2
      } else {
        cell_idx <- which(rownames(shape_emb) == cell)
        mat      <- t(pre_shape[cell_idx, , ])   # L x 2
      }

      # Procrustes-align to Fréchet mean
      rotated <- t(kendall_tpca$reparam_OPA(t(mat), mean_shape))   # L x 2

      # Frobenius norm rescaling — restores original contour scale
      if (norm_rescale && !is.null(meta)) {
        cell_scale <- meta$scale[meta$cell_id == cell]
        if (length(cell_scale) == 1L && is.finite(cell_scale) && cell_scale > 0) {
          rotated <- rotated * cell_scale
        }
      }

      coord_df <- rbind(coord_df, data.frame(
        x    = rotated[, 1],
        y    = rotated[, 2],
        cell = cell,
        bin  = bin_
      ))
    }
  }

  if (nrow(coord_df) == 0L)
    stop("No coordinates assembled — check that scores and TPCA embeddings share cell names.")

  coord_df$bin <- factor(coord_df$bin, levels = bin_labels)

  subtitle <- if (use_reconstruction) recon_msg else
    paste0("Contour type: ", tpca_result$contour_type,
           if (norm_rescale) " \u2022 Frobenius norm-rescaled" else "")

  n_cols <- min(num_bins, 5L)
  lim <- max(abs(c(coord_df$x, coord_df$y)), na.rm = TRUE)

  p <- ggplot2::ggplot(
    coord_df,
    ggplot2::aes(x = x, y = y, group = cell)
  ) +
    ggplot2::geom_polygon(fill = NA, colour = line_colour, linewidth = 0.5) +
    ggplot2::facet_wrap(~ bin, ncol = n_cols) +
   # ggplot2::coord_equal() +
    ggplot2::labs(
      title    = paste0(score_type, cc_idx, " boundary montage"),
      subtitle = subtitle,
      x        = NULL,
      y        = NULL
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::coord_cartesian(
      xlim = c(-lim, lim),
      ylim = c(-lim, lim),
      expand = FALSE
    ) +
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8),
      aspect.ratio = 1
    )

  attr(p, "coord_data") <- coord_df
  p
}
#' Plot Fréchet mean shape history
#'
#' Visualises the evolution of the Fréchet mean shape across iterations of the
#' Kendall TPCA algorithm. Produces two side-by-side panels: a grid of
#' evenly-spaced mean shape outlines (left) and an overlay of the same shapes
#' with alpha scaled by iteration order (right).
#'
#' @param tpca_result A flat list returned by \code{\link{run_tpca}} with
#'   \code{store_history = TRUE}. Must have a \code{"frechet_history"}
#'   attribute containing \code{$mu_history} (a 3-D array,
#'   iterations \eqn{\times} 2 \eqn{\times} landmarks).
#' @param n Integer. Number of evenly-spaced iterawtions to display.
#'   Default \code{9L}.
#' @param colour Fill/outline colour for contours. Default \code{"steelblue"}.
#' @param linewidth Line width passed to \code{geom_path}. Default \code{0.6}.
#'
#' @return A \code{ggplot} object (assembled with \code{ggpubr::ggarrange}).
#' @seealso \code{\link{run_tpca}}, \code{\link{plot_frechet_convergence}}
#' @export
plot_mu_history <- function(tpca_result,
                            n         = 9L,
                            colour    = "steelblue",
                            linewidth = 0.6) {

  history <- attr(tpca_result, "frechet_history")
  if (is.null(history))
    stop(
      "No frechet_history attribute found. ",
      "Re-run run_tpca() with store_history = TRUE."
    )

  mu_history <- history$mu_history          # [iter, landmark, 2]
  n_iter     <- dim(mu_history)[1]
  n          <- min(as.integer(n), n_iter)

  idx   <- unique(round(seq(1, n_iter, length.out = n)))
  alpha <- seq(0.15, 1, length.out = length(idx))

  # Build a long data frame for all selected iterations
  iter_dfs <- Map(function(i, a) {
    coords <- t(mu_history[i, , ])  # landmarks x 2

    coords <- rbind(coords, coords[1L, ])

    data.frame(
      x         = coords[, 1],
      y         = coords[, 2],
      iteration = i,
      iter_lab  = paste0("iter ", i),
      alpha     = a
    )
  }, idx, alpha)

  df <- do.call(rbind, iter_dfs)
  df$iter_lab <- factor(df$iter_lab, levels = unique(df$iter_lab))

  # ---- grid panel -----------------------------------------------------------
  p_grid <- ggplot2::ggplot(
    df,
    ggplot2::aes(x, y, group = iteration)
  ) +
    ggplot2::geom_path(
      colour = colour,
      linewidth = linewidth
    ) +
    ggplot2::facet_wrap(
      ~ iter_lab,
      scales = "fixed"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 7),
      panel.spacing = grid::unit(4, "pt"),
      plot.title = ggplot2::element_text(
        size = 9,
        hjust = 0.5
      )
    ) +
    ggplot2::labs(title = "Mean shape per iteration")

  # ---- overlay panel --------------------------------------------------------
  # geom_path doesn't vectorise alpha per-group cleanly, so draw layer by layer
  p_overlay <- ggplot2::ggplot() +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9, hjust = 0.5)
    ) +
    ggplot2::labs(title = "Overlay (early \u2192 late)")

  for (i in seq_along(idx)) {
    sub <- df[df$iteration == idx[i], ]
    p_overlay <- p_overlay +
      ggplot2::geom_path(
        data      = sub,
        mapping   = ggplot2::aes(x, y),
        colour    = colour,
        alpha     = alpha[i],
        linewidth = linewidth
      )
  }

  ggpubr::ggarrange(p_grid, p_overlay, ncol = 2L, nrow = 1L)
}

#' Plot Fréchet mean convergence
#'
#' Plots the gradient norm at each iteration of the Kendall TPCA Fréchet mean
#' computation. Requires \code{run_tpca()} or \code{run_tpca_from_seurat()}
#' to have been called with \code{store_history = TRUE}.
#'
#' @param result A TPCA result object carrying a \code{"frechet_history"}
#'   attribute (i.e. the return value of \code{run_tpca()} or
#'   \code{load_kstitch_results()} when \code{store_history = TRUE}).
#' @param log_scale Logical. If \code{TRUE} (default), the y-axis is
#'   log10-transformed, which makes exponential decay easier to read.
#' @param point_size Numeric. Size of points. Default \code{1.5}.
#' @param line_width Numeric. Width of the connecting line. Default \code{0.7}.
#' @param title Character. Plot title. Defaults to
#'   \code{"Fréchet mean convergence"}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot_mu_history}}
#' @export
plot_frechet_convergence <- function(
    result,
    log_scale  = TRUE,
    point_size = 1.5,
    line_width = 0.7,
    title      = "Fr\u00e9chet mean convergence"
) {
  history <- attr(result, "frechet_history")
  if (is.null(history)) {
    stop(
      "No 'frechet_history' attribute found. ",
      "Re-run with store_history = TRUE."
    )
  }

  grad_norm <- history[["grad_norm"]]
  if (is.null(grad_norm) || length(grad_norm) == 0) {
    stop("'frechet_history$grad_norm' is missing or empty.")
  }

  df <- data.frame(
    iteration = seq_along(grad_norm),
    grad_norm = grad_norm
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = grad_norm)) +
    ggplot2::geom_line(linewidth = line_width, colour = "grey40") +
    ggplot2::geom_point(size = point_size, colour = "steelblue") +
    ggplot2::labs(
      title = title,
      x     = "Iteration",
      y     = if (log_scale) "Gradient norm (log10)" else "Gradient norm"
    ) +
    ggplot2::theme_bw()

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }

  p
}
