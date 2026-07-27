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
#'   \code{\link{export_seurat_contours}} (e.g. \code{tpca_result$all} or
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

#' Plot CSP-binned cell boundary montage
#'
#' Bins cells by their CSP score into deciles, selects bin medoids (cells with
#' smallest mean pairwise TPCA distance), Procrustes-aligns their contours to
#' the Fréchet mean, and plots them as a montage faceted by bin.
#'
#' When \code{area_rescale = TRUE} (default), normalised contours are rescaled
#' by \code{sqrt(area)} so the displayed size reflects true cell area. Area
#' values are read from the \code{Metadata} element of \code{tpca_result},
#' which is returned directly by \code{\link{run_tpca}} and
#' \code{\link{export_seurat_contours}} — no disk access is required.
#'
#' If the pre-shape embedding is not available in \code{tpca_result}, shapes
#' are reconstructed from the top \code{reconstruction_pcs} PCA components.
#' This is clearly flagged in the plot subtitle and via an R message.
#'
#' @param group_result A single group's entry from the list returned by
#'   \code{\link{link_shape_and_factors}}.
#' @param tpca_result A single group's TPCA result — i.e. one element of the
#'   named list returned by \code{\link{run_tpca}} or
#'   \code{\link{export_seurat_contours}} (e.g. \code{tpca_result$all} or
#'   \code{tpca_result[["keratinocyte"]]}). Do not pass the top-level named
#'   list directly.
#' @param cc_idx Integer. Which CSP component to use for binning. Default 1L.
#' @param num_bins Integer. Number of bins. Default 10L.
#' @param num_cells_for_boundary Integer. Number of medoid cells to display
#'   per bin. Default 6L.
#' @param pcs_to_use Integer vector. Which shape PCs to use for medoid
#'   distance computation. Default \code{1:10}.
#' @param area_rescale Logical. Rescale contours by \code{sqrt(area)} to
#'   reflect true cell size. Default \code{TRUE}.
#' @param reconstruction_pcs Integer. Number of PCA components to use for
#'   shape reconstruction when pre-shape embedding is unavailable.
#'   Default = all available.
#' @param line_colour Character. Colour for contour outlines. Default
#'   \code{"steelblue"}.
#'
#' @return A \code{ggplot} object. The underlying coordinate data frame is
#'   attached as attribute \code{"coord_data"} for downstream use.
#' @export
plot_csp_boundary_montage <- function(group_result,
                                      tpca_result,
                                      cc_idx                 = 1L,
                                      num_bins               = 10L,
                                      num_cells_for_boundary = 6L,
                                      pcs_to_use             = 1:10,
                                      norm_rescale           = FALSE,
                                      reconstruction_pcs     = NULL,
                                      line_colour            = "steelblue") {

  kendall_tpca_py_dir <- system.file("python", package = "kstitch")
  kendall_tpca        <- reticulate::import_from_path("kendall_tpca", path = kendall_tpca_py_dir)

  info       <- tpca_result$Info
  mean_shape <- t(info$frechet_mean)   # 2 x L
  emb        <- group_result$CSP_Scores
  shape_emb  <- tpca_result$TPCA_Embedding
  valid_idx  <- which(complete.cases(shape_emb))
  shape_emb  <- shape_emb[valid_idx, , drop = FALSE]

  pre_shape          <- info$pre_shape_embedding   # cells x 2 x L, or NULL
  use_reconstruction <- is.null(pre_shape)

  if (use_reconstruction) {
    k_rec <- if (is.null(reconstruction_pcs)) ncol(shape_emb) else
      min(reconstruction_pcs, ncol(shape_emb))
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

  # bin cells by CSP score
  cells_with_both <- intersect(rownames(emb), rownames(shape_emb))
  csp_scores      <- emb[cells_with_both, cc_idx]

  bin_labels <- paste0("CSP", cc_idx, " Bin ", seq_len(num_bins))
  bin_vec    <- dplyr::ntile(csp_scores, num_bins)
  bin_df     <- data.frame(
    cell  = cells_with_both,
    score = csp_scores,
    bin   = factor(paste0("CSP", cc_idx, " Bin ", bin_vec), levels = bin_labels)
  )

  pcs_to_use <- pcs_to_use[pcs_to_use <= ncol(shape_emb)]
  coord_df   <- data.frame()

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
        pca_coords <- shape_emb[cell, seq_len(k_rec)]
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
    stop("No coordinates assembled — check that CSP scores and TPCA embeddings share cell names.")

  coord_df$bin <- factor(coord_df$bin, levels = bin_labels)

  subtitle <- if (use_reconstruction) recon_msg else
    paste0("Contour type: ", tpca_result$contour_type,
           if (norm_rescale) " \u2022 Frobenius norm-rescaled" else "")

  p <- ggplot2::ggplot(
    coord_df,
    ggplot2::aes(x = x, y = y, group = cell)
  ) +
    ggplot2::geom_polygon(fill = NA, colour = line_colour, linewidth = 0.5) +
    ggplot2::facet_wrap(~ bin, nrow = 1) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title    = paste0("CSP", cc_idx, " boundary montage"),
      subtitle = subtitle,
      x        = NULL,
      y        = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8)
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
#' @param n Integer. Number of evenly-spaced iterations to display.
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
