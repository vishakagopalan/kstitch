#' Anchor CCA component signs by highest-loading shape feature
#'
#' Ensures canonical component signs are consistent across groups and runs by
#' flipping all score, vector, and structure-correlation matrices so that the
#' shape feature with the highest absolute structure correlation loads
#' positively on each CSP component. The anchor feature is allowed to differ
#' per component and is recorded in the output for auditability.
#'
#' @param group_result A single group's entry from the list returned by
#'   \code{\link{link_shape_and_factors}}, i.e.
#'   \code{link_shape_and_factors(...)[[group]]}.
#'
#' @return The same list with signs corrected in-place across all score,
#'   vector, and structure-correlation matrices, plus a new element:
#'   \describe{
#'     \item{Anchor_Features}{Named character vector, one entry per component
#'       (names: \code{CSP1}, \code{CSP2}, ...). Each value is the name of the
#'       shape feature whose absolute structure correlation was largest for that
#'       component and was used to determine the sign.}
#'   }
#'
#' @seealso \code{\link{link_shape_and_factors}}
#' @export
anchor_cca_signs <- function(group_result) {
  csp_self_cor  <- group_result$Shape_Corr_With_CSP   # shape features x k
  cep_self_cor  <- group_result$Exp_Corr_With_CEP   # NMF factors   x k
  csp_scores    <- group_result$CSP_Scores              # cells x k
  cep_scores    <- group_result$CEP_Scores              # cells x k
  csp_vectors   <- group_result$CSP_Vectors             # shape features x k
  cep_vectors   <- group_result$CEP_Vectors             # NMF factors   x k

  # cross-loading blocks from Misc_CCA
  cep_cor_with_csp <- group_result$Misc_CCA$scores$corr.exp.with.csp  # NMF x k
  csp_cor_with_cep <- group_result$Misc_CCA$scores$corr.shape.with.cep  # shape x k

  k <- ncol(csp_scores)
  anchor_features <- character(k)
  names(anchor_features) <- colnames(csp_scores)   # CSP1, CSP2, ...

  for (cc_idx in seq_len(k)) {

    # feature with highest absolute structure correlation on this component
    abs_cors          <- abs(csp_self_cor[, cc_idx])
    anchor_feature    <- names(which.max(abs_cors))
    anchor_features[cc_idx] <- anchor_feature

    # flip if that feature's raw (signed) correlation is negative
    if (csp_self_cor[anchor_feature, cc_idx] < 0) {
      csp_self_cor[, cc_idx]    <- -csp_self_cor[, cc_idx]
      cep_self_cor[, cc_idx]    <- -cep_self_cor[, cc_idx]
      csp_scores[, cc_idx]      <- -csp_scores[, cc_idx]
      cep_scores[, cc_idx]      <- -cep_scores[, cc_idx]
      csp_vectors[, cc_idx]     <- -csp_vectors[, cc_idx]
      cep_vectors[, cc_idx]     <- -cep_vectors[, cc_idx]
      cep_cor_with_csp[, cc_idx] <- -cep_cor_with_csp[, cc_idx]
      csp_cor_with_cep[, cc_idx] <- -csp_cor_with_cep[, cc_idx]
    }
  }

  group_result$Shape_Corr_With_CSP                    <- csp_self_cor
  group_result$Exp_Corr_With_CEP                    <- cep_self_cor
  group_result$CSP_Scores                               <- csp_scores
  group_result$CEP_Scores                               <- cep_scores
  group_result$CSP_Vectors                              <- csp_vectors
  group_result$CEP_Vectors                              <- cep_vectors
  group_result$Misc_CCA$scores$corr.exp.with.csp           <- cep_cor_with_csp
  group_result$Misc_CCA$scores$corr.shape.with.cep           <- csp_cor_with_cep
  group_result$Anchor_Features                          <- anchor_features

  group_result
}
