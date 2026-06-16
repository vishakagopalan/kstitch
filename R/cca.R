#' Canonical Correlation Analysis via SVD
#'
#' Computes CCA between two matrices using the standard whitening + SVD
#' formulation, implemented in base R. Intended as a drop-in replacement for
#' \code{CCA::rcc()} within kstitch: the return value mirrors the structure of
#' \code{CCA::cc()} / \code{CCA::rcc()} (\code{$cor}, \code{$xcoef},
#' \code{$ycoef}, and a \code{$scores} list holding the canonical variates and
#' the structure-correlation blocks), so downstream code that consumes those
#' fields needs no other changes.
#'
#' @param X Numeric matrix (cells x features), e.g. the shape / morphology
#'   matrix. Rownames should be cell IDs and colnames feature names.
#' @param Y Numeric matrix (cells x features), e.g. the NMF factor matrix.
#'   Must have the same rows (cells) as \code{X}.
#' @param scale Logical; if \code{TRUE} (default) the columns of \code{X} and
#'   \code{Y} are standardized (centered and scaled to unit variance) before the
#'   analysis. If \code{FALSE} the matrices are used as supplied (still centered
#'   internally for the covariance computation). Classical CCA is invariant to
#'   per-column rescaling, so this mainly affects numerical conditioning and the
#'   interpretation of the weight vectors, not the canonical correlations.
#' @param tol Eigenvalues below this threshold are clamped before inversion, to
#'   guard against numerical blow-up during whitening of a near-singular
#'   covariance matrix.
#'
#' @return A list with components:
#'   \describe{
#'     \item{cor}{Canonical correlations, length K = min(ncol(X), ncol(Y)).}
#'     \item{xcoef, ycoef}{Canonical weight vectors (features x K).}
#'     \item{scores}{List with \code{xscores}, \code{yscores} (canonical
#'       variates, cells x K) plus structure correlations
#'       \code{corr.X.xscores}, \code{corr.Y.xscores}, \code{corr.X.yscores},
#'       \code{corr.Y.yscores}.}
#'   }
#' @export
arun_cca <- function(X, Y, scale = TRUE, tol = 1e-8) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows (cells).")
  }

  # Optional standardization (center + unit variance) of each column.
  if (scale) {
    X <- scale(X)
    Y <- scale(Y)
  }

  K <- min(ncol(X), ncol(Y))

  # Center columns explicitly so the canonical variates / scores come out
  # mean-zero with the right variance (harmless if already scaled above).
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)

  # Covariance matrices.
  Sxx <- stats::cov(X)
  Syy <- stats::cov(Y)
  Sxy <- stats::cov(X, Y)

  # Symmetric inverse square root via eigendecomposition, with a floor on the
  # eigenvalues to stay well-conditioned.
  inv_sqrt <- function(S) {
    e <- eigen(S, symmetric = TRUE)
    vals <- e$values
    vals[vals < tol] <- tol
    e$vectors %*% diag(1 / sqrt(vals), length(vals)) %*% t(e$vectors)
  }

  Sxx_isqrt <- inv_sqrt(Sxx)
  Syy_isqrt <- inv_sqrt(Syy)

  # Whiten the cross-covariance; its singular values are the canonical
  # correlations and its singular vectors give the canonical directions.
  M  <- Sxx_isqrt %*% Sxy %*% Syy_isqrt
  sv <- svd(M, nu = K, nv = K)

  cors <- pmin(pmax(sv$d[seq_len(K)], 0), 1)

  # Back-transform singular vectors into canonical weight vectors in the
  # original feature spaces (this is what CCA::rcc returns as xcoef / ycoef).
  xcoef <- Sxx_isqrt %*% sv$u
  ycoef <- Syy_isqrt %*% sv$v
  rownames(xcoef) <- colnames(X)
  rownames(ycoef) <- colnames(Y)
  colnames(xcoef) <- paste0("CC", seq_len(K))
  colnames(ycoef) <- paste0("CC", seq_len(K))

  # Canonical variates (scores).
  xscores <- Xc %*% xcoef
  yscores <- Yc %*% ycoef
  rownames(xscores) <- rownames(X)
  rownames(yscores) <- rownames(Y)

  # Structure correlations (loadings) and cross-loadings, matching the names
  # CCA::cc uses inside $scores.
  list(
    cor   = cors,
    xcoef = xcoef,
    ycoef = ycoef,
    scores = list(
      xscores        = xscores,
      yscores        = yscores,
      corr.X.xscores = stats::cor(Xc, xscores),
      corr.Y.xscores = stats::cor(Yc, xscores),
      corr.X.yscores = stats::cor(Xc, yscores),
      corr.Y.yscores = stats::cor(Yc, yscores)
    )
  )
}
