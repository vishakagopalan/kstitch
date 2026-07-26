# Shared test helper: creates a minimal real Seurat object for testing.
# Uses SeuratObject::CreateSeuratObject() which is available without the
# full Seurat package.
make_seurat_stub <- function(cell_names, meta_cols = NULL, n_features = 10) {
  n      <- length(cell_names)
  counts <- matrix(
    rpois(n_features * n, lambda = 2), nrow = n_features, ncol = n,
    dimnames = list(paste0("gene", seq_len(n_features)), cell_names)
  )
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  if (!is.null(meta_cols)) {
    for (col_name in names(meta_cols)) {
      obj[[col_name]] <- meta_cols[[col_name]]
    }
  }
  obj
}
