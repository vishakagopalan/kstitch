# data-raw/prepare_vignette_data.R
#
# Builds the pre-computed kstitch vignette dataset from saved melanoma Xenium
# Prime objects. Run this script once to produce:
#
#   inst/extdata/kstitch_vignette_xenium.rds
#
# The saved object is a 3,171-cell keratinocyte Seurat v5 slice with TPCA,
# NMF, and CCA results stored inside via kstitch storage helpers. The vignette
# loads this object directly so users do not need Python, rliger, or large
# data files to run it.
#
# This script loads pre-computed results from RDS files saved outside the
# kstitch repo. Update the file paths below to point to your saved objects.

library(kstitch)
library(Seurat)
library(dplyr)

# ---- File paths (update these to point to your saved RDS files) -----------

input_files <- list(
  xenium_obj      = "~/Downloads/Temp_Keratinocyte_Object.rds",
  nmf_output      = "~/Downloads/Temp_Melanoma_NMF.rds",
  pga_embedding   = "~/Downloads/Temp_Melanoma_PGA_Out.rds",
  pga_info        = "~/Downloads/Temp_Melanoma_PGA_Info.rds",
  cca_cell        = "~/Downloads/Temp_Melanoma_CCA_Cells.rds",
  cca_nucleus     = "~/Downloads/Temp_Melanoma_CCA_Nuclei.rds"
)

# Verify all files exist
for (name in names(input_files)) {
  path <- input_files[[name]]
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s (%s)", name, path))
  }
}

cat("Loading input data...\n")
xenium_obj              <- readRDS(input_files$xenium_obj)
nmf_output              <- readRDS(input_files$nmf_output)
pga_embedding_list      <- readRDS(input_files$pga_embedding)
pga_info_list           <- readRDS(input_files$pga_info)
cca_cell_output         <- readRDS(input_files$cca_cell)
cca_nucleus_output      <- readRDS(input_files$cca_nucleus)

# ---- 1. Use the xenium object as-is (already keratinocyte subset) ----------

obj <- xenium_obj
cat(sprintf("Starting with: %d cells\n", ncol(obj)))

# ---- 2. Store TPCA results (Cell) -------------------------------------------

tpca_cell <- list(
  PGA_Embedding = pga_embedding_list$Cell$Keratinocytes,
  Info = list(
    variances           = pga_info_list$Cell$variances,
    v_matrix            = pga_info_list$Cell$v_matrix,
    frechet_mean        = pga_info_list$Cell$frechet_mean,
    pre_shape_embedding = NULL   # not shipped to keep file size small
  ),
  contour_type = "cell",
  output_dir   = NA_character_   # no temp dir — precomputed
)

obj <- store_tpca_results(obj, tpca_cell)
cat(sprintf("TPCA cell stored: %d cells with embeddings\n",
            sum(!is.na(Seurat::Embeddings(obj, "tpca_cell")[, 1]))))

# ---- 3. Store TPCA results (Nucleus) ----------------------------------------

tpca_nucleus <- list(
  PGA_Embedding = pga_embedding_list$Nucleus$Keratinocytes,
  Info = list(
    variances           = pga_info_list$Nucleus$variances,
    v_matrix            = pga_info_list$Nucleus$v_matrix,
    frechet_mean        = pga_info_list$Nucleus$frechet_mean,
    pre_shape_embedding = NULL
  ),
  contour_type = "nucleus",
  output_dir   = NA_character_
)

obj <- store_tpca_results(obj, tpca_nucleus)
cat(sprintf("TPCA nucleus stored: %d cells with embeddings\n",
            sum(!is.na(Seurat::Embeddings(obj, "tpca_nucleus")[, 1]))))

# ---- 4. Store NMF results ---------------------------------------------------

nmf_mat     <- nmf_output$Keratinocytes$NMF_Matrix
nmf_loading <- nmf_output$Keratinocytes$NMF_Loading

# Store as DimReduc
colnames(nmf_mat) <- paste0("NMF_", seq_len(ncol(nmf_mat)))
obj[["nmf"]] <- Seurat::CreateDimReducObject(
  embeddings = nmf_mat,
  key        = "NMF_",
  assay      = Seurat::DefaultAssay(obj)
)

# Store loading matrix in misc for downstream use
obj@misc$kstitch$nmf <- list(
  NMF_Loading      = nmf_loading,
  NMF_Factor_Names = colnames(nmf_loading)
)

cat(sprintf("NMF stored: %d cells x %d factors\n",
            nrow(nmf_mat), ncol(nmf_mat)))

# ---- 5. Store CCA results (Cell) --------------------------------------------

old_cca_cell <- cca_cell_output
k_cell       <- length(old_cca_cell$Keratinocytes$CC_Corr_Coefs)

csp_scores_cell  <- old_cca_cell$Keratinocytes$Shape_Canonical_Scores
cep_scores_cell  <- old_cca_cell$Keratinocytes$Factor_Canonical_Scores
csp_vectors_cell <- old_cca_cell$Keratinocytes$Shape_Canonical_Vectors
cep_vectors_cell <- old_cca_cell$Keratinocytes$Factor_Canonical_Vectors

colnames(csp_scores_cell)  <- paste0("CSP", seq_len(k_cell))
colnames(cep_scores_cell)  <- paste0("CEP", seq_len(k_cell))
colnames(csp_vectors_cell) <- paste0("CSP", seq_len(k_cell))
colnames(cep_vectors_cell) <- paste0("CEP", seq_len(k_cell))

misc_cca_cell <- old_cca_cell$Misc_CCA
misc_cca_cell$scores$corr.X.xscores <- old_cca_cell$Keratinocytes$Shape_CC_Self_Correlations
misc_cca_cell$scores$corr.Y.yscores <- old_cca_cell$Keratinocytes$Factor_CC_Self_Correlations
misc_cca_cell$scores$corr.X.yscores  <- old_cca_cell$Keratinocytes$Misc_CCA$scores$corr.X.yscores
misc_cca_cell$scores$corr.Y.xscores  <- old_cca_cell$Keratinocytes$Misc_CCA$scores$corr.Y.xscores

cca_result_cell <- list(
  CC_Corr_Coefs         = old_cca_cell$CC_Corr_Coefs,
  CSP_Scores            = csp_scores_cell,
  CEP_Scores            = cep_scores_cell,
  CSP_Vectors           = csp_vectors_cell,
  CEP_Vectors           = cep_vectors_cell,
  CSP_Self_Correlations = old_cca_cell$Keratinocytes$Shape_CC_Self_Correlations,
  CEP_Self_Correlations = old_cca_cell$Keratinocytes$Factor_CC_Self_Correlations,
  Misc_CCA              = misc_cca_cell
)

if (!is.null(old_cca_cell$Keratinocytes$Prioritized_Variables)) {
  cca_result_cell$Priority_Info <- old_cca_cell$Keratinocytes$Prioritized_Variables
}

cca_result_cell <- anchor_cca_signs(cca_result_cell)

obj@misc[["kstitch"]][["cell_cca"]] <- cca_result_cell[
  setdiff(names(cca_result_cell), c("CSP_Scores", "CEP_Scores"))
]

cat(sprintf("CCA (cell) stored: %d canonical components\n", k_cell))

# ---- 6. Store CCA results (Nucleus) -----------------------------------------

old_cca_nucleus <- cca_nucleus_output
k_nucleus       <- length(old_cca_nucleus$Keratinocytes$CC_Corr_Coefs)

csp_scores_nuc  <- old_cca_nucleus$Keratinocytes$Shape_Canonical_Scores
cep_scores_nuc  <- old_cca_nucleus$Keratinocytes$Factor_Canonical_Scores
csp_vectors_nuc <- old_cca_nucleus$Keratinocytes$Shape_Canonical_Vectors
cep_vectors_nuc <- old_cca_nucleus$Keratinocytes$Factor_Canonical_Vectors

colnames(csp_scores_nuc)  <- paste0("CSP", seq_len(k_nucleus))
colnames(cep_scores_nuc)  <- paste0("CEP", seq_len(k_nucleus))
colnames(csp_vectors_nuc) <- paste0("CSP", seq_len(k_nucleus))
colnames(cep_vectors_nuc) <- paste0("CEP", seq_len(k_nucleus))

misc_cca_nuc <- old_cca_nucleus$Keratinocytes$Misc_CCA
misc_cca_nuc$scores$corr.X.xscores <- old_cca_nucleus$Keratinocytes$Shape_CC_Self_Correlations
misc_cca_nuc$scores$corr.Y.yscores <- old_cca_nucleus$Keratinocytes$Factor_CC_Self_Correlations
misc_cca_nuc$scores$corr.X.yscores  <- old_cca_nucleus$Keratinocytes$Misc_CCA$scores$corr.X.yscores
misc_cca_nuc$scores$corr.Y.xscores  <- old_cca_nucleus$Keratinocytes$Misc_CCA$scores$corr.Y.xscores

cca_result_nuc <- list(
  CC_Corr_Coefs         = old_cca_nucleus$CC_Corr_Coefs,
  CSP_Scores            = csp_scores_nuc,
  CEP_Scores            = cep_scores_nuc,
  CSP_Vectors           = csp_vectors_nuc,
  CEP_Vectors           = cep_vectors_nuc,
  CSP_Self_Correlations = old_cca_nucleus$Keratinocytes$Shape_CC_Self_Correlations,
  CEP_Self_Correlations = old_cca_nucleus$Keratinocytes$Factor_CC_Self_Correlations,
  Misc_CCA              = misc_cca_nuc
)

if (!is.null(old_cca_nucleus$Keratinocytes$Prioritized_Variables)) {
  cca_result_nuc$Priority_Info <- old_cca_nucleus$Keratinocytes$Prioritized_Variables
}

cca_result_nuc <- anchor_cca_signs(cca_result_nuc)

obj@misc[["kstitch"]][["nucleus_cca"]] <- cca_result_nuc[
  setdiff(names(cca_result_nuc), c("CSP_Scores", "CEP_Scores"))
]

cat(sprintf("CCA (nucleus) stored: %d canonical components\n", k_nucleus))

# ---- 7. Trim and save -------------------------------------------------------

obj <- DietSeurat(obj, assays = Seurat::Assays(obj),
                  dimreducs = c("tpca_cell", "tpca_nucleus", "nmf"))

out_dir  <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "kstitch_vignette_xenium.rds")

saveRDS(obj, out_path)
cat(sprintf("\nSaved to: %s\n", normalizePath(out_path)))
cat(sprintf("File size: %.1f MB\n", file.size(out_path) / 1e6))
