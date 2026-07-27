# Tests for compute_nmf. Focuses on logic that does not require a real rliger
# call: the features argument, early-exit path, return structure, and
# return_results = FALSE. Full integration tests (actual NMF fit) require real
# data and are left for manual / CI testing with the full dependency stack.
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
# features argument validation
# ---------------------------------------------------------------------------

test_that("compute_nmf errors when features is not NULL or character", {
  obj <- make_seurat_stub(paste0("cell_", 1:60))
  expect_error(
    compute_nmf(obj, assay_name = "RNA", features = 42, verbose = FALSE),
    "must be NULL or a character vector"
  )
})

test_that("compute_nmf warns when fewer than 100 supplied genes are present", {
  obj   <- make_seurat_stub(paste0("cell_", 1:60), n_features = 200)
  genes <- paste0("notreal_", 1:200)   # none present in stub
  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          expect_warning(
            compute_nmf(obj, assay_name = "RNA", features = genes, verbose = FALSE),
            "very small|unreliable"
          )
        }
      )
    }
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

test_that("compute_nmf passes a character vector to .run_nmf_one_group", {
  obj            <- make_seurat_stub(paste0("cell_", 1:60), n_features = 200)
  genes          <- paste0("gene", 1:150)
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
          suppressWarnings(
            compute_nmf(obj, assay_name = "RNA", features = genes, verbose = FALSE)
          )
        }
      )
    }
  )

  expect_equal(captured_feats, genes)
})

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

test_that("compute_nmf returns a flat list with expected fields", {
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

  expect_true(is.list(res))
  expect_true(all(c("Factor_Gene_List", "NMF_Matrix",
                    "Fit_Error", "Fit_Summary") %in% names(res)))
  # No group / is_groupwise fields on flat result
  expect_false("group" %in% names(res))
  expect_false("is_groupwise" %in% names(res))
})

# ---------------------------------------------------------------------------
# return_results = FALSE
# ---------------------------------------------------------------------------

test_that("compute_nmf return_results = FALSE writes a single RDS and returns its path", {
  obj     <- make_seurat_stub(paste0("cell_", 1:60))
  out_dir <- file.path(tempdir(), paste0("test_nmf_", kstitch:::.random_id()))

  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(x, ...) x,
        .package = "base",
        {
          res <- compute_nmf(obj, assay_name = "RNA", verbose = FALSE,
                             return_results = FALSE, output_dir = out_dir)
        }
      )
    }
  )

  expect_named(res, "output_path")
  expect_true(file.exists(res$output_path))
  expect_equal(basename(res$output_path), "all.rds")

  loaded <- load_kstitch_results(res$output_path, type = "nmf")
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
# .resolve_features_for_group (still used internally)
# ---------------------------------------------------------------------------

test_that(".resolve_features_for_group returns NULL when features is NULL", {
  expect_null(kstitch:::.resolve_features_for_group(NULL, "A"))
})

test_that(".resolve_features_for_group returns the vector for character features", {
  genes <- paste0("g", 1:100)
  expect_equal(kstitch:::.resolve_features_for_group(genes, "A"), genes)
})
