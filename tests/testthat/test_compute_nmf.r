# Tests for compute_nmf focus on logic that does not require a real rliger
# call: group splitting, the early-exit path, features argument handling, and
# return structure validation. Full integration tests (actual NMF fit) require
# real data and are left for manual / CI testing with the full dependency stack.
#
# Uses make_seurat_stub() from helper-seurat-stub.R for all obj construction.
# .run_nmf_one_group is mocked throughout so no rliger calls are made.

empty_fit <- list(
  Factor_Gene_List = list(),
  NMF_Matrix       = NULL,
  NMF_Loading      = NULL,
  Fit_Error        = NA_real_,
  Fit_Summary      = data.frame(
    Method               = "online",
    k                    = NA_integer_,
    Stability_Score      = NA_real_,
    Entropy              = NA_real_,
    Reconstruction_Error = NA_real_,
    N_Stability_Runs     = 0L,
    stringsAsFactors     = FALSE
  )
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_grouped_stub <- function(groups = c("TypeA", "TypeB"), n_per_group = 50,
                              n_features = 200) {
  n_cells    <- n_per_group * length(groups)
  cell_names <- paste0("cell_", seq_len(n_cells))
  meta_cols  <- list(celltype = rep(groups, each = n_per_group))
  make_seurat_stub(cell_names, meta_cols = meta_cols, n_features = n_features)
}

# ---------------------------------------------------------------------------
# group.by validation
# ---------------------------------------------------------------------------

test_that("compute_nmf errors if group.by column is absent", {
  obj <- make_seurat_stub(paste0("cell_", 1:50))
  expect_error(
    compute_nmf(obj, assay_name = "RNA", group.by = "nonexistent"),
    "not found in obj@meta.data"
  )
})

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

test_that("compute_nmf returns a list keyed by group when group.by is supplied", {
  obj <- make_grouped_stub()

  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          res <- compute_nmf(obj, assay_name = "RNA", group.by = "celltype",
                             verbose = FALSE)
        }
      )
    }
  )

  expect_named(res, c("TypeA", "TypeB"), ignore.order = TRUE)
  expect_true(all(vapply(res, is.list, logical(1))))
})

test_that("compute_nmf returns a single 'all' group when group.by is NULL", {
  obj <- make_seurat_stub(paste0("cell_", 1:60))

  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          res <- compute_nmf(obj, assay_name = "RNA", group.by = NULL,
                             verbose = FALSE)
        }
      )
    }
  )

  expect_named(res, "all")
})

test_that("compute_nmf attaches group and is_groupwise fields — single group", {
  obj <- make_seurat_stub(paste0("cell_", 1:60))

  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          res <- compute_nmf(obj, assay_name = "RNA", verbose = FALSE)
        }
      )
    }
  )

  expect_equal(res$all$group, "all")
  expect_false(res$all$is_groupwise)
})

test_that("compute_nmf attaches group and is_groupwise fields — grouped", {
  obj <- make_grouped_stub()

  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          res <- compute_nmf(obj, assay_name = "RNA", group.by = "celltype",
                             verbose = FALSE)
        }
      )
    }
  )

  expect_equal(res$TypeA$group, "TypeA")
  expect_equal(res$TypeB$group, "TypeB")
  expect_true(res$TypeA$is_groupwise)
  expect_true(res$TypeB$is_groupwise)
})

# ---------------------------------------------------------------------------
# return_results = FALSE
# ---------------------------------------------------------------------------

test_that("compute_nmf return_results = FALSE writes RDS files and returns paths", {
  obj     <- make_grouped_stub()
  out_dir <- file.path(tempdir(), paste0("test_nmf_", kstitch:::.random_id()))

  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          res <- compute_nmf(obj, assay_name = "RNA", group.by = "celltype",
                             verbose = FALSE, return_results = FALSE,
                             output_dir = out_dir)
        }
      )
    }
  )

  expect_named(res, "output_paths")
  expect_named(res$output_paths, c("TypeA", "TypeB"), ignore.order = TRUE)
  expect_true(all(file.exists(res$output_paths)))

  loaded <- load_kstitch_results(res$output_paths[["TypeA"]], type = "nmf")
  expect_true(is.list(loaded))

  unlink(out_dir, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# .compute_kotliar_stability
# ---------------------------------------------------------------------------

test_that(".compute_kotliar_stability returns NA for fewer than 2 runs", {
  expect_true(all(is.na(kstitch:::.compute_kotliar_stability(list(), k = 3))))
  expect_true(all(is.na(kstitch:::.compute_kotliar_stability(
    list(matrix(1:6, 3, 2)), k = 2
  ))))
})

test_that(".compute_kotliar_stability returns named numeric for valid input", {
  make_W <- function(seed) {
    set.seed(seed)
    matrix(abs(rnorm(30 * 5)), 30, 5,
           dimnames = list(paste0("g", 1:30), paste0("F", 1:5)))
  }
  result <- kstitch:::.compute_kotliar_stability(lapply(1:4, make_W), k = 5)

  expect_named(result, c("silhouette", "entropy"))
  expect_true(is.numeric(result))
  expect_true(result[["silhouette"]] >= -1 && result[["silhouette"]] <= 1)
})

# ---------------------------------------------------------------------------
# features argument — .validate_features_arg
# ---------------------------------------------------------------------------

test_that(".validate_features_arg passes silently when features is NULL", {
  expect_silent(
    kstitch:::.validate_features_arg(NULL, c("A", "B"), paste0("g", 1:500))
  )
})

test_that(".validate_features_arg passes for a plain character vector", {
  genes <- paste0("g", 1:150)
  expect_silent(
    kstitch:::.validate_features_arg(genes, c("A", "B"), paste0("g", 1:500))
  )
})

test_that(".validate_features_arg passes for a correctly named list", {
  genes <- paste0("g", 1:150)
  expect_silent(
    kstitch:::.validate_features_arg(
      list(A = genes, B = genes),
      c("A", "B"),
      paste0("g", 1:500)
    )
  )
})

test_that(".validate_features_arg errors when named list is missing a group", {
  genes <- paste0("g", 1:150)
  expect_error(
    kstitch:::.validate_features_arg(
      list(A = genes),          # B is missing
      c("A", "B"),
      paste0("g", 1:500)
    ),
    "do not match group levels"
  )
})

test_that(".validate_features_arg errors when named list has extra group names", {
  genes <- paste0("g", 1:150)
  expect_error(
    kstitch:::.validate_features_arg(
      list(A = genes, B = genes, C = genes),   # C is extra
      c("A", "B"),
      paste0("g", 1:500)
    ),
    "do not match group levels"
  )
})

test_that(".validate_features_arg warns when fewer than 100 genes present — character vector", {
  # Only 50 of the supplied genes exist in the object
  supplied <- paste0("g", 1:200)
  present  <- paste0("g", 1:50)
  expect_warning(
    kstitch:::.validate_features_arg(supplied, c("A"), present),
    "very small"
  )
})

test_that(".validate_features_arg warns when fewer than 100 genes present — named list", {
  all_features <- paste0("g", 1:500)
  # Group A has enough genes; group B has only 30 present
  feat_list <- list(
    A = paste0("g", 1:150),
    B = paste0("x", 1:200)   # none of these are in all_features
  )
  expect_warning(
    kstitch:::.validate_features_arg(feat_list, c("A", "B"), all_features),
    "very small"
  )
})

test_that(".validate_features_arg does not warn when >= 100 genes present", {
  genes <- paste0("g", 1:150)
  expect_silent(
    kstitch:::.validate_features_arg(genes, c("A"), paste0("g", 1:500))
  )
})

# ---------------------------------------------------------------------------
# features argument — .resolve_features_for_group
# ---------------------------------------------------------------------------

test_that(".resolve_features_for_group returns NULL when features is NULL", {
  expect_null(kstitch:::.resolve_features_for_group(NULL, "A"))
})

test_that(".resolve_features_for_group returns the vector for a character features", {
  genes <- paste0("g", 1:100)
  expect_equal(kstitch:::.resolve_features_for_group(genes, "A"), genes)
  expect_equal(kstitch:::.resolve_features_for_group(genes, "B"), genes)
})

test_that(".resolve_features_for_group returns the group-specific vector for a named list", {
  feat_list <- list(A = paste0("g", 1:100), B = paste0("g", 101:200))
  expect_equal(kstitch:::.resolve_features_for_group(feat_list, "A"), feat_list$A)
  expect_equal(kstitch:::.resolve_features_for_group(feat_list, "B"), feat_list$B)
})

# ---------------------------------------------------------------------------
# features argument — compute_nmf integration (mocked)
# ---------------------------------------------------------------------------

test_that("compute_nmf errors on invalid features type", {
  obj <- make_seurat_stub(paste0("cell_", 1:60))
  expect_error(
    compute_nmf(obj, assay_name = "RNA", features = 42, verbose = FALSE),
    "must be NULL, a character vector, or a named list"
  )
})

test_that("compute_nmf errors when named features list misses a group", {
  obj   <- make_grouped_stub(n_features = 200)
  genes <- paste0("gene", 1:150)
  expect_error(
    compute_nmf(obj, assay_name = "RNA", group.by = "celltype",
                features = list(TypeA = genes),   # TypeB missing
                verbose = FALSE),
    "do not match group levels"
  )
})

test_that("compute_nmf passes features = NULL to .run_nmf_one_group", {
  obj            <- make_seurat_stub(paste0("cell_", 1:60))
  captured_feats <- NULL

  mock_run <- function(...) {
    args <- list(...)
    captured_feats <<- args[["features"]]
    empty_fit
  }

  with_mocked_bindings(
    .run_nmf_one_group = mock_run,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          compute_nmf(obj, assay_name = "RNA", features = NULL, verbose = FALSE)
        }
      )
    }
  )

  expect_null(captured_feats)
})

test_that("compute_nmf passes a character vector features to every group", {
  obj            <- make_grouped_stub(n_features = 200)
  genes          <- paste0("gene", 1:150)
  captured_feats <- list()

  mock_run <- function(...) {
    args <- list(...)
    captured_feats[[length(captured_feats) + 1L]] <<- args[["features"]]
    empty_fit
  }

  with_mocked_bindings(
    .run_nmf_one_group = mock_run,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          compute_nmf(obj, assay_name = "RNA", group.by = "celltype",
                      features = genes, verbose = FALSE)
        }
      )
    }
  )

  expect_length(captured_feats, 2L)
  expect_equal(captured_feats[[1]], genes)
  expect_equal(captured_feats[[2]], genes)
})

test_that("compute_nmf passes per-group features correctly", {
  obj   <- make_grouped_stub(n_features = 200)
  feats <- list(TypeA = paste0("gene", 1:150), TypeB = paste0("gene", 51:200))
  captured <- list()

  mock_run <- function(...) {
    args <- list(...)
    captured[[length(captured) + 1L]] <<- args[["features"]]
    empty_fit
  }

  with_mocked_bindings(
    .run_nmf_one_group = mock_run,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          res <- compute_nmf(obj, assay_name = "RNA", group.by = "celltype",
                             features = feats, verbose = FALSE)
        }
      )
    }
  )

  # group order matches names(groups) from split(), i.e. alphabetical
  expect_setequal(
    c(list(feats$TypeA), list(feats$TypeB)),
    captured
  )
})
