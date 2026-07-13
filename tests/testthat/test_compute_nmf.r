# Tests for compute_nmf focus on logic that does not require a real rliger
# call: group splitting, the early-exit path, and return structure validation.
# Full integration tests (actual NMF fit) require real data and are left for
# manual / CI testing with the full dependency stack.

# Minimal stub obj that responds to [[ for meta.data without dispatching
# through SeuratObject S4 methods.
make_nmf_stub <- function(cell_names, group_vec = NULL) {
  meta <- if (!is.null(group_vec)) {
    data.frame(celltype = group_vec, row.names = cell_names)
  } else {
    data.frame(row.names = cell_names)
  }
  list(meta.data = meta)
}

test_that("compute_nmf errors if group.by column is absent", {

  obj <- make_nmf_stub(paste0("cell_", 1:50))

  expect_error(
    compute_nmf(obj, assay_name = "RNA", group.by = "nonexistent"),
    "not found in obj@meta.data"
  )
})

test_that("compute_nmf returns a list keyed by group when group.by is supplied", {

  # We mock .run_nmf_one_group so no rliger calls are made.
  cells     <- paste0("cell_", 1:100)
  grp_vec   <- rep(c("TypeA", "TypeB"), each = 50)
  obj       <- make_nmf_stub(cells, grp_vec)

  empty_fit <- list(
    Factor_Gene_List = list(),
    NMF_Matrix       = NULL,
    NMF_Loading      = NULL,
    Fit_Error        = NA_real_,
    Fit_Summary      = data.frame(
      Method = "online", k = NA_integer_,
      Stability_Score = NA_real_, Entropy = NA_real_,
      Reconstruction_Error = NA_real_, N_Stability_Runs = 0L,
      stringsAsFactors = FALSE
    )
  )

  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(obj, ...) obj,
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

  cells <- paste0("cell_", 1:60)
  obj   <- make_nmf_stub(cells)

  empty_fit <- list(
    Factor_Gene_List = list(), NMF_Matrix = NULL, NMF_Loading = NULL,
    Fit_Error = NA_real_,
    Fit_Summary = data.frame(Method = "online", k = NA_integer_,
                             Stability_Score = NA_real_, Entropy = NA_real_,
                             Reconstruction_Error = NA_real_,
                             N_Stability_Runs = 0L, stringsAsFactors = FALSE)
  )

  with_mocked_bindings(
    .run_nmf_one_group = function(...) empty_fit,
    .package = "kstitch",
    {
      with_mocked_bindings(
        subset = function(obj, ...) obj,
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

test_that(".compute_kotliar_stability returns NA for fewer than 2 runs", {

  expect_true(all(is.na(kstitch:::.compute_kotliar_stability(list(), k = 3))))
  expect_true(all(is.na(kstitch:::.compute_kotliar_stability(list(matrix(1:6, 3, 2)), k = 2))))
})

test_that(".compute_kotliar_stability returns named numeric for valid input", {

  set.seed(1)
  make_W <- function(seed) {
    set.seed(seed)
    W <- matrix(abs(rnorm(30 * 5)), 30, 5,
                dimnames = list(paste0("g", 1:30), paste0("F", 1:5)))
    W
  }
  W_list <- lapply(1:4, make_W)
  result <- kstitch:::.compute_kotliar_stability(W_list, k = 5)

  expect_named(result, c("silhouette", "entropy"))
  expect_true(is.numeric(result))
  expect_true(result[["silhouette"]] >= -1 && result[["silhouette"]] <= 1)
})
