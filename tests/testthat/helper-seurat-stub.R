# Shared test helper: creates a minimal real Seurat object for testing.
# Uses SeuratObject::CreateSeuratObject() which is available without the
# full Seurat package.

make_seurat_stub <- function(cell_names, meta_cols = NULL) {
  n     <- length(cell_names)
  # minimal count matrix: 10 genes x n cells
  counts <- matrix(
    rpois(10 * n, lambda = 2), nrow = 10, ncol = n,
    dimnames = list(paste0("gene", seq_len(10)), cell_names)
  )
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))

  if (!is.null(meta_cols)) {
    for (col_name in names(meta_cols)) {
      obj[[col_name]] <- meta_cols[[col_name]]
    }
  }

  obj
}
