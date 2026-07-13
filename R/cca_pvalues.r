#' Deflation-based permutation p-values for canonical correlations
#'
#' Computes per-component permutation p-values for CCA using sequential
#' deflation: after each accepted canonical variate pair, the shared variance
#' is projected out of both matrices before testing the next component. This
#' gives p-values that are conditional on earlier components, making them
#' more interpretable than marginal permutation p-values.
#'
#' Built entirely on \code{\link{run_cca}}; does not use \code{CCA::rcc},
#' \code{CCA::cc}, or any external CCA implementation.
#'
#' @param X Numeric matrix, cells x shape features. Should already be scaled
#'   if scaling is desired (e.g. pass the same \code{scale(X)} used for
#'   \code{run_cca}).
#' @param Y Numeric matrix, cells x NMF factors. Same scaling note as \code{X}.
#' @param K Integer. Number of canonical components to test. Defaults to
#'   \code{min(ncol(X), ncol(Y))}.
#' @param nperm Integer. Number of permutations. Default 1000.
#' @param statistic Character. Test statistic: \code{"abs_cor"} (default) or
#'   \code{"cor"}.
#' @param method Character. \code{"stepdown"} (default) enforces monotone
#'   p-values via cumulative max across components. \code{"unadjusted"} returns
#'   raw per-component p-values.
#' @param blocks Optional integer or factor vector of length \code{nrow(X)}.
#'   When supplied, rows are permuted within blocks rather than globally (useful
#'   for permuting within FOV or slide).
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{pvalues}{Named numeric vector of length K (names: CC1, CC2, ...).}
#'     \item{raw_p}{Raw (pre-monotone) p-values. Identical to \code{pvalues}
#'       when \code{method = "unadjusted"}.}
#'     \item{obs_cancor}{Observed deflated canonical correlations, length K.}
#'     \item{counts}{Exceedance counts, length K.}
#'     \item{method}{The method used.}
#'     \item{meta}{List of call parameters for provenance.}
#'   }
#'
#' @seealso \code{\link{run_cca}}, \code{\link{link_shape_and_factors}}
#' @export
cca_pvalues <- function(X, Y,
                        K          = NULL,
                        nperm      = 1000L,
                        statistic  = c("abs_cor", "cor"),
                        method     = c("stepdown", "unadjusted"),
                        blocks     = NULL,
                        seed       = NULL,
                        verbose    = TRUE) {

  statistic <- match.arg(statistic)
  method    <- match.arg(method)

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  n <- nrow(X)
  if (nrow(Y) != n) stop("X and Y must have the same number of rows.")

  p <- ncol(X); q <- ncol(Y)
  if (is.null(K)) K <- min(p, q)
  K <- min(K, p, q)

  nperm <- as.integer(nperm)
  if (nperm < 1L) stop("`nperm` must be >= 1.")

  # --- block permutation setup ------------------------------------------------

  if (!is.null(blocks)) {
    if (length(blocks) != n) stop("`blocks` must have length equal to nrow(X).")
    blocks       <- as.factor(blocks)
    block_idx    <- split(seq_len(n), blocks)
  } else {
    block_idx <- NULL
  }

  permute_rows <- function() {
    if (is.null(block_idx)) return(sample.int(n))
    idx <- integer(n)
    for (ids in block_idx) {
      idx[ids] <- if (length(ids) == 1L) ids else ids[sample.int(length(ids))]
    }
    idx
  }

  # --- deflation --------------------------------------------------------------
  # Runs CCA on (Xk, Yk), extracts the first canonical variate pair, projects
  # it out, and returns the first canonical correlation. Repeats K times.
  # Returns a numeric vector of length K (NA if a step fails numerically).

  deflated_cancors <- function(Xmat, Ymat) {
    Xk      <- Xmat
    Yk      <- Ymat
    cancors <- rep(NA_real_, K)

    for (k in seq_len(K)) {
      if (ncol(Xk) < 1L || ncol(Yk) < 1L) break

      fit <- tryCatch(
        run_cca(Xk, Yk, scale = FALSE, tol = .Machine$double.eps^0.5),
        error = function(e) NULL
      )
      if (is.null(fit) || length(fit$cor) == 0L) break

      cancors[k] <- fit$cor[1L]

      # first canonical variate pair
      u <- as.vector(Xk %*% fit$xcoef[, 1L])
      v <- as.vector(Yk %*% fit$ycoef[, 1L])

      uu <- sum(u * u)
      vv <- sum(v * v)

      # stop if variates are numerically zero
      eps <- .Machine$double.eps
      if (uu <= eps || vv <= eps) break

      # rank-1 deflation: remove the projection of u (v) from Xk (Yk)
      Xk <- Xk - tcrossprod(u, crossprod(u, Xk) / uu)
      Yk <- Yk - tcrossprod(v, crossprod(v, Yk) / vv)
    }

    cancors
  }

  # --- observed canonical correlations ----------------------------------------

  if (verbose) message("Computing observed deflated canonical correlations ...")
  obs_cancor <- deflated_cancors(X, Y)

  abs_flag   <- statistic == "abs_cor"
  obs_stat   <- if (abs_flag) abs(obs_cancor) else obs_cancor

  # --- permutation loop -------------------------------------------------------

  if (!is.null(seed)) set.seed(seed)

  ge_counts <- integer(K)

  if (verbose) message(sprintf("Running %d permutations ...", nperm))

  for (i in seq_len(nperm)) {
    if (verbose && i %% 100L == 0L) message(sprintf("  permutation %d / %d", i, nperm))

    perm_idx   <- permute_rows()
    perm_stat  <- deflated_cancors(X, Y[perm_idx, , drop = FALSE])
    perm_stat  <- if (abs_flag) abs(perm_stat) else perm_stat

    if (method == "stepdown") {
      # suffix maximum: each position gets the max of itself and all later ones,
      # so exceedance at position k means the permutation beat the observed
      # statistic conditional on earlier components
      perm_stat[is.na(perm_stat)] <- -Inf
      if (K >= 2L) {
        for (j in (K - 1L):1L) {
          perm_stat[j] <- max(perm_stat[j], perm_stat[j + 1L])
        }
      }
    }

    exceeded     <- !is.na(perm_stat) & (perm_stat >= obs_stat)
    ge_counts    <- ge_counts + as.integer(exceeded)
  }

  # --- p-values ---------------------------------------------------------------

  # Laplace smoothing: avoids p = 0
  raw_p  <- (as.numeric(ge_counts) + 1L) / (nperm + 1L)
  pvals  <- if (method == "stepdown") cummax(raw_p) else raw_p
  names(pvals) <- names(raw_p) <- paste0("CC", seq_len(K))

  list(
    pvalues    = pvals,
    raw_p      = raw_p,
    obs_cancor = obs_cancor,
    counts     = ge_counts,
    method     = method,
    meta       = list(
      n         = n,
      p         = p,
      q         = q,
      K         = K,
      nperm     = nperm,
      statistic = statistic,
      blocks    = !is.null(blocks),
      seed      = seed
    )
  )
}
