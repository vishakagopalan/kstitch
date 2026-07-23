make_group_result <- function(n = 150, k_shape = 8, k_nmf = 6, seed = 1) {
  set.seed(seed)
  Z         <- matrix(rnorm(n * 3), n, 3)
  X         <- scale(Z %*% matrix(rnorm(3 * k_shape), 3, k_shape) + matrix(rnorm(n * k_shape, sd = 0.3), n, k_shape))
  Y         <- scale(Z %*% matrix(rnorm(3 * k_nmf),   3, k_nmf)   + matrix(rnorm(n * k_nmf,   sd = 0.3), n, k_nmf))
  colnames(X) <- paste0("ShapePC", seq_len(k_shape))
  colnames(Y) <- paste0("Factor",  seq_len(k_nmf))
  rownames(X) <- rownames(Y) <- paste0("cell_", seq_len(n))

  cca <- run_cca(X, Y, scale = FALSE)
  k   <- length(cca$cor)

  csp_scores <- cca$scores$xscores;  colnames(csp_scores) <- paste0("CSP", seq_len(k)); rownames(csp_scores) <- rownames(X)
  cep_scores <- cca$scores$yscores;  colnames(cep_scores) <- paste0("CEP", seq_len(k)); rownames(cep_scores) <- rownames(Y)
  csp_vectors <- cca$xcoef;          colnames(csp_vectors) <- paste0("CSP", seq_len(k))
  cep_vectors <- cca$ycoef;          colnames(cep_vectors) <- paste0("CEP", seq_len(k))

  names(cca$scores)[match(c("corr.X.xscores", "corr.X.yscores",
                          "corr.Y.xscores", "corr.Y.yscores"),
                        names(cca$scores))] <-
      c("corr.shape.with.csp", "corr.shape.with.cep",
        "corr.exp.with.csp", "corr.exp.with.cep")


  list(
    CC_Corr_Coefs         = cca$cor,
    CSP_Scores            = csp_scores,
    CEP_Scores            = cep_scores,
    CSP_Vectors           = csp_vectors,
    CEP_Vectors           = cep_vectors,
    Shape_Corr_With_CSP  = cca$scores$corr.shape.with.csp,
    Exp_Corr_With_CEP  = cca$scores$corr.exp.with.cep,
    Misc_CCA              = cca
  )
}

test_that("anchor_cca_signs returns Anchor_Features with correct names and values", {

  grp <- make_group_result()
  out <- anchor_cca_signs(grp)

  k <- length(out$CC_Corr_Coefs)
  expect_named(out, c(names(grp), "Anchor_Features"), ignore.order = TRUE)
  expect_length(out$Anchor_Features, k)
  expect_equal(names(out$Anchor_Features), paste0("CSP", seq_len(k)))

  # each anchor feature name must be a real shape feature
  shape_features <- rownames(out$Shape_Corr_With_CSP)
  expect_true(all(out$Anchor_Features %in% shape_features))
})

test_that("after anchoring, anchor feature loads positively on each component", {

  grp <- make_group_result(seed = 7)
  out <- anchor_cca_signs(grp)

  k <- length(out$CC_Corr_Coefs)
  for (cc_idx in seq_len(k)) {
    anchor <- out$Anchor_Features[cc_idx]
    expect_gte(out$Shape_Corr_With_CSP[anchor, cc_idx], 0)
  }
})

test_that("anchor feature is the one with highest absolute structure correlation", {

  grp <- make_group_result(seed = 3)
  out <- anchor_cca_signs(grp)

  k <- length(out$CC_Corr_Coefs)
  for (cc_idx in seq_len(k)) {
    anchor      <- out$Anchor_Features[cc_idx]
    abs_cors    <- abs(out$Shape_Corr_With_CSP[, cc_idx])
    expect_equal(unname(as.character(anchor)), names(which.max(abs_cors)))
  }
})

test_that("anchor_cca_signs is idempotent: running twice gives the same result", {

  grp  <- make_group_result(seed = 5)
  out1 <- anchor_cca_signs(grp)
  out2 <- anchor_cca_signs(out1)

  expect_equal(out1$CSP_Scores,  out2$CSP_Scores)
  expect_equal(out1$CEP_Scores,  out2$CEP_Scores)
  expect_equal(out1$Anchor_Features, out2$Anchor_Features)
})
