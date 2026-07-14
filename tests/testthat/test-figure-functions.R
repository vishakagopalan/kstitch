# Tests for figure functions focus on input validation, return types, and
# correct use of group_result fields. Reticulate/Python calls are mocked.
# Visual correctness is verified manually.

make_group_result_fig <- function(n = 120, k_shape = 8, k_nmf = 6,
                                  k_cc = 4, seed = 1) {
  set.seed(seed)
  cells <- paste0("cell_", seq_len(n))

  csp_scores <- matrix(rnorm(n * k_cc), n, k_cc,
                       dimnames = list(cells, paste0("CSP", seq_len(k_cc))))
  cep_scores <- matrix(rnorm(n * k_cc), n, k_cc,
                       dimnames = list(cells, paste0("CEP", seq_len(k_cc))))
  csp_vec    <- matrix(rnorm(k_shape * k_cc), k_shape, k_cc,
                       dimnames = list(paste0("ShapePC", seq_len(k_shape)),
                                       paste0("CSP", seq_len(k_cc))))
  cep_vec    <- matrix(rnorm(k_nmf * k_cc), k_nmf, k_cc,
                       dimnames = list(paste0("Factor", seq_len(k_nmf)),
                                       paste0("CEP", seq_len(k_cc))))

  list(
    CC_Corr_Coefs         = runif(k_cc, 0.5, 0.9),
    CSP_Scores            = csp_scores,
    CEP_Scores            = cep_scores,
    CSP_Vectors           = csp_vec,
    CEP_Vectors           = cep_vec,
    CSP_Self_Correlations = matrix(rnorm(k_shape * k_cc), k_shape, k_cc,
                                   dimnames = dimnames(csp_vec)),
    CEP_Self_Correlations = matrix(rnorm(k_nmf * k_cc), k_nmf, k_cc,
                                   dimnames = dimnames(cep_vec)),
    Anchor_Features       = setNames(paste0("ShapePC", seq_len(k_cc)),
                                     paste0("CSP", seq_len(k_cc))),
    Misc_CCA              = list()
  )
}

make_tpca_result_fig <- function(cell_names, L = 50, k = 8,
                                 contour_type = "cell",
                                 include_pre_shape = TRUE) {
  emb <- matrix(rnorm(length(cell_names) * k), length(cell_names), k,
                dimnames = list(cell_names, paste0("Shape_PC", seq_len(k))))

  pre_shape <- if (include_pre_shape) {
    array(rnorm(length(cell_names) * 2 * L),
          dim = c(length(cell_names), 2, L))
  } else NULL

  list(
    TPCA_Embedding = emb,
    Info = list(
      variances      = sort(runif(k, 0.5, 5), decreasing = TRUE),
      v_matrix       = matrix(rnorm(2 * L * k), 2 * L, k),
      frechet_mean   = matrix(rnorm(2 * L), 2, L),
      pre_shape_embedding = pre_shape
    ),
    contour_type = contour_type,
    output_dir   = tempdir()
  )
}

# ---- plot_cep_loadings ------------------------------------------------------

test_that("plot_cep_loadings returns a ggplot with correct labels", {

  grp <- make_group_result_fig()
  p   <- plot_cep_loadings(grp, cc_idx = 1L)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Loading on CEP1")
  expect_equal(nrow(p$data), nrow(grp$CEP_Vectors))
})

test_that("plot_cep_loadings works for cc_idx > 1", {

  grp <- make_group_result_fig()
  expect_s3_class(plot_cep_loadings(grp, cc_idx = 3L), "ggplot")
})

# ---- plot_csp_loadings ------------------------------------------------------

test_that("plot_csp_loadings returns a ggplot with anchor in subtitle", {

  grp <- make_group_result_fig()
  p   <- plot_csp_loadings(grp, cc_idx = 1L)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Loading on CSP1")
  expect_true(grepl(grp$Anchor_Features[["CSP1"]], p$labels$subtitle))
  expect_equal(nrow(p$data), nrow(grp$CSP_Vectors))
})

# ---- plot_shape_modes -------------------------------------------------------

test_that("plot_shape_modes returns a ggplot (reticulate mocked)", {

  cells      <- paste0("cell_", 1:80)
  tpca_res   <- make_tpca_result_fig(cells)
  sds        <- c(-2, 0, 2)
  num_pcs    <- 3L

  # mock reticulate calls
  fake_df <- do.call(rbind, lapply(paste0("PC", seq_len(num_pcs)), function(pc) {
    do.call(rbind, lapply(sds, function(sd) {
      data.frame(x = rnorm(50), y = rnorm(50), PC = pc, sd_val = sd)
    }))
  }))

  with_mocked_bindings(
    import_from_path = function(...) {
      list(reconstruct_shapes_from_pca = function(...) fake_df)
    },
    .package = "reticulate",
    {
      p <- plot_shape_modes(tpca_res, sds_to_plot = sds, num_pcs = num_pcs)
    }
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Shape modes (TPCA)")
})

# ---- plot_csp_boundary_montage ----------------------------------------------

test_that("plot_csp_boundary_montage returns ggplot with coord_data attribute", {

  set.seed(5)
  cells    <- paste0("cell_", 1:120)
  grp      <- make_group_result_fig(n = 120, seed = 5)
  tpca_res <- make_tpca_result_fig(cells, include_pre_shape = TRUE)

  L          <- dim(tpca_res$Info$pre_shape_embedding)[3]
  mean_shape <- tpca_res$Info$frechet_mean   # 2 x L

  with_mocked_bindings(
    import_from_path = function(...) {
      list(
        reparam_OPA = function(A, B) A   # identity rotation
      )
    },
    .package = "reticulate",
    {
      p <- plot_csp_boundary_montage(
        grp, tpca_res,
        cc_idx = 1L, num_bins = 5L,
        num_cells_for_boundary = 3L,
        area_rescale = FALSE
      )
    }
  )

  expect_s3_class(p, "ggplot")
  expect_true(!is.null(attr(p, "coord_data")))
  expect_true(all(c("x", "y", "cell", "bin") %in%
                    colnames(attr(p, "coord_data"))))
})

test_that("plot_csp_boundary_montage messages and uses reconstruction when pre_shape is NULL", {

  set.seed(6)
  cells    <- paste0("cell_", 1:120)
  grp      <- make_group_result_fig(n = 120, seed = 6)
  tpca_res <- make_tpca_result_fig(cells, include_pre_shape = FALSE)

  L <- ncol(tpca_res$Info$frechet_mean)

  with_mocked_bindings(
    import_from_path = function(...) {
      list(
        reparam_OPA = function(A, B) A,
        reconstruct_shape_from_pca_coords = function(...) {
          matrix(rnorm(2 * L), 2, L)
        }
      )
    },
    .package = "reticulate",
    {
      expect_message(
        p <- plot_csp_boundary_montage(
          grp, tpca_res,
          cc_idx = 1L, num_bins = 5L,
          num_cells_for_boundary = 3L,
          area_rescale = FALSE
        ),
        "reconstructed"
      )
    }
  )

  expect_s3_class(p, "ggplot")
  expect_true(grepl("reconstructed", p$labels$subtitle))
})
