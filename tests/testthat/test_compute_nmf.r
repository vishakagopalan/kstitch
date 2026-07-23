# Tests for compute_nmf focus on logic that does not require a real rliger
# call: group splitting, the early-exit path, and return structure validation.
# Full integration tests (actual NMF fit) require real data and are left for
# manual / CI testing with the full dependency stack.
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

test_that("compute_nmf errors if group.by column is absent", {
  obj <- make_seurat_stub(paste0("cell_", 1:50))
  expect_error(
    compute_nmf(obj, assay_name = "RNA", group.by = "nonexistent"),
    "not found in obj@meta.data"
  )
})

test_that("compute_nmf returns a list keyed by group when group.by is supplied", {
  obj <- make_seurat_stub(
    cell_names = paste0("cell_", 1:100),
    meta_cols  = list(celltype = rep(c("TypeA", "TypeB"), each = 50))
  )

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



# ---------------------------------------------------------------------------
# .compute_kotliar_stability -- unchanged, no Seurat dependency
# ---------------------------------------------------------------------------

test_that(".compute_kotliar_stability returns NA for fewer than 2 runs", {
  expect_true(all(is.na(kstitch:::.compute_kotliar_stability(list(), k = 3))))
  expect_true(all(is.na(kstitch:::.compute_kotliar_stability(list(matrix(1:6, 3, 2)), k = 2))))
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
