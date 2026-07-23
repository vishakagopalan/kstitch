# TODO: Once NMF results are stored as DimReducs inside the Seurat object
# (consistent with the TPCA bridge in PR #6), the explicit return of NMF_Matrix
# and NMF_Loading should be retired in favour of reading from the enriched
# object directly. For now results are returned as a named list and the caller
# is responsible for passing NMF_Matrix to link_shape_and_factors().

# ---- private helpers --------------------------------------------------------

.strip_reserved <- function(x, reserved) {
  if (length(x) == 0L) return(x)
  x[setdiff(names(x), reserved)]
}

.extract_obj_err <- function(obj, reduction_name) {
  out <- NA_real_
  try({ out <- obj[[reduction_name]]@misc$objErr }, silent = TRUE)
  out
}

.get_feature_loadings <- function(obj, reduction_name) {
  obj[[reduction_name]]@feature.loadings
}

.get_embeddings <- function(obj, reduction_name) {
  Seurat::Embeddings(obj, reduction = reduction_name)
}

.compute_kotliar_stability <- function(W_list, k) {
  if (length(W_list) < 2L) return(c(silhouette = NA_real_, entropy = NA_real_))

  W_list <- Filter(Negate(is.null), W_list)
  W_list <- lapply(W_list, function(W) {
    W <- as.matrix(W)
    storage.mode(W) <- "double"
    W[!is.finite(W)] <- 0
    W
  })

  common_genes <- Reduce(intersect, lapply(W_list, rownames))
  if (length(common_genes) < 2L) return(c(silhouette = NA_real_, entropy = NA_real_))

  W_list <- lapply(W_list, function(W) {
    W <- W[common_genes, , drop = FALSE]
    keep <- apply(W, 2, function(x) stats::sd(x) > 0 && all(is.finite(x)))
    W[, keep, drop = FALSE]
  })

  W_list <- Filter(function(W) is.matrix(W) && ncol(W) > 0L, W_list)
  if (length(W_list) < 2L) return(c(silhouette = NA_real_, entropy = NA_real_))

  W_list <- lapply(W_list, function(W) {
    norms <- sqrt(colSums(W^2))
    norms[!is.finite(norms) | norms == 0] <- 1
    sweep(W, 2, norms, "/")
  })

  pooled <- do.call(cbind, W_list)
  storage.mode(pooled) <- "double"
  if (ncol(pooled) < 2L) return(c(silhouette = NA_real_, entropy = NA_real_))

  sim <- stats::cor(pooled, use = "pairwise.complete.obs")
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  sim[sim < 0] <- 0

  dist_obj <- stats::as.dist(1 - sim)
  hc       <- stats::hclust(dist_obj, method = "average")
  k_use    <- min(as.integer(k), ncol(pooled))
  if (k_use < 2L) return(c(silhouette = NA_real_, entropy = NA_real_))

  cl  <- stats::cutree(hc, k = k_use)
  sil <- cluster::silhouette(cl, dist_obj)

  run_ids <- rep(seq_along(W_list), vapply(W_list, ncol, integer(1L)))
  entropy <- mean(vapply(unique(cl), function(c) {
    p <- table(run_ids[cl == c]) / sum(cl == c)
    -sum(p * log(p))
  }, numeric(1L))) / log(length(unique(run_ids)))

  c(silhouette = mean(sil[, "sil_width"], na.rm = TRUE), entropy = entropy)
}

.run_one_online_fit <- function(obj, k, mini_batch_size, seed, reduction, params, nCores) {
  params <- .strip_reserved(params, c("k", "lambda", "seed", "nCores", "verbose",
                                      "datasetVar", "layer", "assay", "reduction"))
  do.call(rliger::runOnlineINMF,
          c(list(obj, k = k, minibatchSize = mini_batch_size,
                 seed = seed, nCores = nCores, reduction = reduction),
            params))
}

# ---- per-group NMF ----------------------------------------------------------

.run_nmf_one_group <- function(obj,
                               assay_name,
                               num_top_genes_per_factor,
                               batch_var,
                               count_matrix_layer,
                               min_cells,
                               gene_pct_threshold,
                               use_normalized_factor_scores,
                               default_num_factors,
                               nCores,
                               normalize_params,
                               selectGenes_params,
                               scaleNotCenter_params,
                               runCINMF_params,
                               quantileNorm_params,
                               integration_method,
                               stability_n_runs,
                               stability_seed,
                               runOnlineINMF_params) {

  empty_result <- function() {
    list(
      Factor_Gene_List     = list(),
      NMF_Matrix           = NULL,
      NMF_Loading          = NULL,
      Fit_Error            = NA_real_,
      Fit_Summary          = data.frame(
        Method               = integration_method,
        k                    = NA_integer_,
        Stability_Score      = NA_real_,
        Entropy              = NA_real_,
        Reconstruction_Error = NA_real_,
        N_Stability_Runs     = 0L,
        stringsAsFactors     = FALSE
      )
    )
  }

  # gene filter
  data_mat          <- obj[[assay_name]]@layers[[count_matrix_layer]]
  rownames(data_mat) <- SeuratObject::Features(obj)
  colnames(data_mat) <- Seurat::Cells(obj)
  genes_to_use      <- names(which(Matrix::rowMeans(data_mat > 0) >= gene_pct_threshold))
  obj               <- subset(obj, features = genes_to_use)
  rm(data_mat)

  if (length(Seurat::Cells(obj)) < min_cells) return(empty_result())

  mini_batch_size <- if (length(Seurat::Cells(obj)) > 5000) 5000L else
    max(2L, floor(length(Seurat::Cells(obj)) / 2))

  # preprocess
  if (!"layer" %in% names(normalize_params))
    normalize_params[["layer"]] <- count_matrix_layer

  obj <- do.call(rliger::normalize, c(list(obj), normalize_params))

  for (param_list in c("selectGenes_params", "scaleNotCenter_params",
                       "runCINMF_params", "runOnlineINMF_params",
                       "quantileNorm_params")) {
    if (!is.null(batch_var))
      assign(param_list, c(get(param_list), list(datasetVar = batch_var)))
  }

  if (identical(selectGenes_params[["nGenes"]], "all"))
    selectGenes_params[["nGenes"]] <- length(genes_to_use)

  obj <- do.call(rliger::selectGenes,    c(list(obj), selectGenes_params))
  obj <- do.call(rliger::scaleNotCenter, c(list(obj), scaleNotCenter_params))

  reduction_name  <- "nmf"
  n_factors       <- default_num_factors
  fit_ok          <- FALSE

  # fit with k fallback
  while (!fit_ok) {
    fit_ok  <- TRUE
    obj_try <- tryCatch({
      if (integration_method == "online") {
        .run_one_online_fit(obj, n_factors, mini_batch_size,
                            seed   = runOnlineINMF_params[["seed"]] %||% stability_seed,
                            reduction = reduction_name,
                            params = runOnlineINMF_params,
                            nCores = nCores)
      } else {
        rp <- .strip_reserved(
          if (length(runCINMF_params) > 0L) runCINMF_params else runOnlineINMF_params,
          c("k", "lambda", "seed", "nCores", "verbose",
            "datasetVar", "layer", "assay", "reduction")
        )
        do.call(rliger::runCINMF,
                c(list(obj, k = n_factors, minibatchSize = mini_batch_size,
                       seed   = runCINMF_params[["seed"]] %||% stability_seed,
                       nCores = nCores, reduction = reduction_name),
                  rp))
      }
    }, error = function(e) { message(e$message); NULL })

    if (is.null(obj_try)) {
      fit_ok <- FALSE
      if (n_factors <= 1L) stop("NMF could not be computed for this group.")
      message(sprintf("Fit failed with k = %d; retrying with k = %d.",
                      n_factors, n_factors - 1L))
      n_factors <- n_factors - 1L
    } else {
      obj <- obj_try
    }
  }

  # stability runs
  stability_score   <- c(silhouette = NA_real_, entropy = NA_real_)
  replicate_errors  <- numeric(0)

  if (stability_n_runs >= 2L) {
    W_list <- vector("list", stability_n_runs)
    for (i in seq_len(stability_n_runs)) {
      rep_reduction <- paste0(reduction_name, "_stab_", i)
      rep_obj <- tryCatch(
        .run_one_online_fit(obj, n_factors, mini_batch_size,
                            seed = stability_seed + i - 1L,
                            reduction = rep_reduction,
                            params = runOnlineINMF_params,
                            nCores = nCores),
        error = function(e) NULL
      )
      if (!is.null(rep_obj)) {
        W_list[[i]]        <- .get_feature_loadings(rep_obj, rep_reduction)
        replicate_errors[i] <- .extract_obj_err(rep_obj, rep_reduction)
      }
    }
    W_list <- Filter(Negate(is.null), W_list)
    if (length(W_list) >= 2L)
      stability_score <- .compute_kotliar_stability(W_list, n_factors)
  }

  # quantile normalisation
  obj <- do.call(rliger::quantileNorm,
                 c(list(obj, reduction = reduction_name), quantileNorm_params))
  norm_reduction <- paste0(reduction_name, "Norm")

  # extract
  if (use_normalized_factor_scores) {
    nmf_mat     <- .get_embeddings(obj, norm_reduction)
    nmf_loading <- .get_feature_loadings(obj, norm_reduction)
    colnames(nmf_mat)     <- gsub(norm_reduction, "Factor", colnames(nmf_mat))
    colnames(nmf_loading) <- gsub(norm_reduction, "Factor", colnames(nmf_loading))
  } else {
    nmf_mat     <- .get_embeddings(obj, reduction_name)
    nmf_loading <- .get_feature_loadings(obj, reduction_name)
    colnames(nmf_mat)     <- gsub(reduction_name, "Factor", colnames(nmf_mat))
    colnames(nmf_loading) <- gsub(reduction_name, "Factor", colnames(nmf_loading))
  }

  n_top        <- min(num_top_genes_per_factor, nrow(nmf_loading))
  gene_list    <- lapply(seq_len(ncol(nmf_loading)), function(j)
    names(sort(nmf_loading[, j], decreasing = TRUE)[seq_len(n_top)]))
  names(gene_list) <- paste0("Factor_", seq_along(gene_list))

  fit_error <- .extract_obj_err(obj, reduction_name)
  if (is.na(fit_error) && length(replicate_errors) > 0L)
    fit_error <- mean(replicate_errors, na.rm = TRUE)

  list(
    Factor_Gene_List = gene_list,
    NMF_Matrix       = nmf_mat,
    NMF_Loading      = nmf_loading,
    Fit_Error        = fit_error,
    Fit_Summary      = data.frame(
      Method               = integration_method,
      k                    = n_factors,
      Stability_Score      = stability_score[["silhouette"]],
      Entropy              = stability_score[["entropy"]],
      Reconstruction_Error = fit_error,
      N_Stability_Runs     = if (stability_n_runs >= 2L) stability_n_runs else 0L,
      stringsAsFactors     = FALSE
    )
  )
}

# ---- exported function ------------------------------------------------------

#' Run LIGER NMF factorization per group
#'
#' Preprocesses a Seurat object, runs online iNMF or consensus iNMF via
#' \pkg{rliger}, and returns factor loadings, cell embeddings, top genes per
#' factor, and a fit-summary table (including a Kotliar-style stability score).
#' When \code{group.by} is supplied, factorization runs independently for each
#' unique value of that metadata column.
#'
#' @param obj A Seurat v5 object.
#' @param assay_name Name of the assay to use.
#' @param group.by Character scalar naming a column in \code{obj@meta.data} to
#'   stratify by. \code{NULL} (default) pools all cells.
#' @param num_top_genes_per_factor Number of top genes to record per factor.
#'   Default 100.
#' @param batch_var Optional metadata column for dataset/batch labels passed to
#'   rliger functions.
#' @param count_matrix_layer Layer name for raw counts. Default \code{"counts"}.
#' @param min_cells Minimum cells required to run factorization for a group.
#'   Default 100.
#' @param gene_pct_threshold Minimum fraction of cells expressing a gene.
#'   Default 0.01.
#' @param use_normalized_factor_scores If \code{TRUE} (default), return
#'   quantile-normalised embeddings and loadings.
#' @param default_num_factors Starting k. Default 30.
#' @param nCores Number of cores. Default 8.
#' @param normalize_params Extra arguments for \code{rliger::normalize()}.
#' @param selectGenes_params Extra arguments for \code{rliger::selectGenes()}.
#' @param scaleNotCenter_params Extra arguments for \code{rliger::scaleNotCenter()}.
#' @param runCINMF_params Extra arguments for \code{rliger::runCINMF()}.
#' @param quantileNorm_params Extra arguments for \code{rliger::quantileNorm()}.
#' @param integration_method \code{"consensus"} (default) or \code{"online"}.
#' @param stability_n_runs Number of repeated online fits for stability scoring.
#'   Default 10.
#' @param stability_seed Starting seed for stability runs. Default 1.
#' @param runOnlineINMF_params Extra arguments for \code{rliger::runOnlineINMF()}.
#' @param verbose Logical. Print group progress. Default \code{TRUE}.
#'
#' @return A named list, one entry per group (or \code{"all"} when
#'   \code{group.by = NULL}), each containing:
#'   \describe{
#'     \item{Factor_Gene_List}{Named list of top genes per factor.}
#'     \item{NMF_Matrix}{Cell x factor loading matrix.}
#'     \item{NMF_Loading}{Gene x factor loading matrix.}
#'     \item{Fit_Error}{Reconstruction error for the final fit.}
#'     \item{Fit_Summary}{One-row data frame: Method, k, Stability_Score,
#'       Entropy, Reconstruction_Error, N_Stability_Runs.}
#'   }
#'
#' @importFrom stats hclust as.dist cutree cor median sd
#' @importFrom cluster silhouette
#' @export
compute_nmf <- function(obj,
                        assay_name,
                        group.by                  = NULL,
                        num_top_genes_per_factor  = 100,
                        batch_var                 = NULL,
                        count_matrix_layer        = "counts",
                        min_cells                 = 100,
                        gene_pct_threshold        = 0.01,
                        use_normalized_factor_scores = TRUE,
                        default_num_factors       = 30,
                        nCores                    = 8,
                        normalize_params          = list(),
                        selectGenes_params        = list(),
                        scaleNotCenter_params     = list(),
                        runCINMF_params           = list(),
                        quantileNorm_params       = list(),
                        integration_method        = c("consensus", "online"),
                        stability_n_runs          = 10L,
                        stability_seed            = 1L,
                        runOnlineINMF_params      = list(),
                        verbose                   = TRUE) {

  integration_method <- match.arg(integration_method)

  if (!is.null(group.by) && !group.by %in% colnames(obj@meta.data)) {
    stop(sprintf("`group.by` column '%s' not found in obj@meta.data.", group.by))
  }

  all_cells <- rownames(obj@meta.data)

  if (is.null(group.by)) {
    groups <- list(all = all_cells)
  } else {
    meta_vec <- obj@meta.data[all_cells, group.by]
    groups   <- split(all_cells, meta_vec)
  }

  results <- list()

  for (grp in names(groups)) {
    if (verbose) message(sprintf("Running NMF for group '%s' ...", grp))

    grp_cells <- groups[[grp]]
    grp_obj   <- subset(obj, cells = grp_cells)

    results[[grp]] <- .run_nmf_one_group(
      obj                          = grp_obj,
      assay_name                   = assay_name,
      num_top_genes_per_factor     = num_top_genes_per_factor,
      batch_var                    = batch_var,
      count_matrix_layer           = count_matrix_layer,
      min_cells                    = min_cells,
      gene_pct_threshold           = gene_pct_threshold,
      use_normalized_factor_scores = use_normalized_factor_scores,
      default_num_factors          = default_num_factors,
      nCores                       = nCores,
      normalize_params             = normalize_params,
      selectGenes_params           = selectGenes_params,
      scaleNotCenter_params        = scaleNotCenter_params,
      runCINMF_params              = runCINMF_params,
      quantileNorm_params          = quantileNorm_params,
      integration_method           = integration_method,
      stability_n_runs             = stability_n_runs,
      stability_seed               = stability_seed,
      runOnlineINMF_params         = runOnlineINMF_params
    )
  }

  results
}
