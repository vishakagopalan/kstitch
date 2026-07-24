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
library(lme4)

# ---- File paths (update these to point to your saved RDS files) -----------

input_files <- list(
  xenium_obj                        = "~/Downloads/Temp_Keratinocyte_Object.rds",
  nmf_output                        = "~/Downloads/Temp_Melanoma_NMF.rds",
  tpca_embedding                    = "~/Downloads/Temp_Melanoma_PGA_Out.rds",
  tpca_info                         = "~/Downloads/Temp_Melanoma_PGA_Info.rds",
  nucleus_meta_data_df              = "~/Downloads/For_Vignette_Nucleus_TPCA/Shape_Metadata.csv.gz",
  cell_meta_data_df                 = "~/Downloads/For_Vignette_Cell_TPCA/Shape_Metadata.csv.gz",
  keratinocyte_genes                = "~/Downloads/Temp_Keratinocyte_Genes.rds",
  keratinocytes_with_single_nucleus = "~/Downloads/Temp_Keratinocyte_Single_Nuclei.rds"
)

# Verify all files exist
for (name in names(input_files)) {
  path <- input_files[[name]]
  if (!file.exists(path))
    stop(sprintf("File not found: %s (%s)", name, path))
}

cat("Loading input data...\n")
xenium_obj                       <- readRDS(input_files$xenium_obj)
nmf_output                       <- readRDS(input_files$nmf_output)
tpca_embedding_list              <- readRDS(input_files$tpca_embedding)
tpca_info_list                   <- readRDS(input_files$tpca_info)
keratinocytes_with_single_nuclei <- readRDS(input_files$keratinocytes_with_single_nucleus)
keratinocyte_genes               <- readRDS(input_files$keratinocyte_genes)
nucleus_meta_data_df             <- read.csv(input_files$nucleus_meta_data_df)
cell_meta_data_df                <- read.csv(input_files$cell_meta_data_df)

xenium_obj <- subset(xenium_obj, features = keratinocyte_genes)

# ---- 1. Use the xenium object as-is (already keratinocyte subset) ----------

obj <- xenium_obj
cat(sprintf("Starting with: %d cells\n", ncol(obj)))

# ---- 2. Store TPCA results (Cell) ------------------------------------------

tpca_cell <- load_kstitch_results("~/Downloads/For_Vignette_Cell_TPCA/", type = "tpca")
tpca_cell$all$is_groupwise <- FALSE
tpca_cell$all$contour_type <- "cell"

obj <- store_tpca_results(obj, tpca_cell, reduction_key_prefix = "CellShapePC_")
cat(sprintf("TPCA cell stored: %d cells with embeddings\n",
            sum(!is.na(Seurat::Embeddings(obj, "tpca_cell")[, 1]))))

# ---- 3. Store TPCA results (Nucleus) ---------------------------------------

tpca_nucleus <- load_kstitch_results("~/Downloads/For_Vignette_Nucleus_TPCA/", type = "tpca")
tpca_nucleus$all$is_groupwise <- FALSE
tpca_nucleus$all$contour_type <- "nucleus"

obj <- store_tpca_results(obj, tpca_nucleus, reduction_key_prefix = "NucleusShapePC_")
cat(sprintf("TPCA nucleus stored: %d cells with embeddings\n",
            sum(!is.na(Seurat::Embeddings(obj, "tpca_nucleus")[, 1]))))

# ---- 4. Store NMF results --------------------------------------------------

nmf_mat     <- nmf_output$Keratinocytes$NMF_Matrix
nmf_loading <- nmf_output$Keratinocytes$NMF_Loading

colnames(nmf_mat) <- paste0("NMF_", seq_len(ncol(nmf_mat)))
obj[["nmf_all"]] <- Seurat::CreateDimReducObject(
  embeddings = nmf_mat,
  key        = "NMF_",
  assay      = Seurat::DefaultAssay(obj)
)

obj@misc$kstitch$nmf[["all"]] <- list(
  NMF_Loading      = nmf_loading,
  NMF_Factor_Names = colnames(nmf_loading)
)

cat(sprintf("NMF stored: %d cells x %d factors\n", nrow(nmf_mat), ncol(nmf_mat)))

# ---- 5. Covariate regression helper ----------------------------------------

# Regresses scaled_log_depth and segmentation_method (random effect) from each
# column of mat using lmer. NA rows are excluded from model fitting; residuals
# are placed back into the correct positions.
.regress_lmer <- function(mat, scaled_log_depth, seg_method) {
  cell_ids  <- rownames(mat)
  resid_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat),
                      dimnames = dimnames(mat))
  for (j in seq_len(ncol(mat))) {
    df <- data.frame(
      y                = mat[, j],
      scaled_log_depth = scaled_log_depth[cell_ids],
      seg_method       = seg_method[cell_ids],
      stringsAsFactors = FALSE
    )
    complete_idx <- complete.cases(df)
    if (sum(complete_idx) < 10L) {
      warning(sprintf("Column %d: fewer than 10 complete cases, skipping.", j))
      next
    }
    fit <- lme4::lmer(
      y ~ scaled_log_depth + (1 | seg_method),
      data = df[complete_idx, ], REML = FALSE
    )
    resid_mat[complete_idx, j] <- residuals(fit)
  }
  resid_mat
}

# ---- 6. Prepare covariate vectors ------------------------------------------

meta <- obj@meta.data

# Shared cell IDs: present in both TPCA reduction and NMF
tpca_cell_cells   <- rownames(Seurat::Embeddings(obj, "tpca_cell"))
tpca_nuc_cells    <- rownames(Seurat::Embeddings(obj, "tpca_nucleus"))
nmf_cells         <- rownames(Seurat::Embeddings(obj, "nmf_all"))
shared_cell_cells <- intersect(tpca_cell_cells, nmf_cells)
shared_nuc_cells  <- intersect(tpca_nuc_cells,  nmf_cells)

# Log area from Shape_Metadata (column is "area"), aligned by cell_id
cell_log_area_vec <- setNames(
  log(cell_meta_data_df$area[match(rownames(meta), cell_meta_data_df$cell_id)]),
  rownames(meta)
)
nuc_log_area_vec <- setNames(
  log(nucleus_meta_data_df$area[match(rownames(meta), nucleus_meta_data_df$cell_id)]),
  rownames(meta)
)

# Log sequencing depth
log_depth_vec <- setNames(log(meta$nCount_Xenium), rownames(meta))

# Segmentation method
seg_method_vec <- setNames(meta$segmentation_method, rownames(meta))

# Scale within shared cell sets
scaled_cell_log_area  <- setNames(scale(cell_log_area_vec[shared_cell_cells])[, 1],  shared_cell_cells)
scaled_nuc_log_area   <- setNames(scale(nuc_log_area_vec[shared_nuc_cells])[, 1],    shared_nuc_cells)
scaled_log_depth_cell <- setNames(scale(log_depth_vec[shared_cell_cells])[, 1],      shared_cell_cells)
scaled_log_depth_nuc  <- setNames(scale(log_depth_vec[shared_nuc_cells])[, 1],       shared_nuc_cells)

# ---- 7. Build shape and expression matrices (Cell) -------------------------

shape_pcs_cell <- Seurat::Embeddings(obj, "tpca_cell")[shared_cell_cells, 1:10]
expr_mat_cell  <- Seurat::Embeddings(obj, "nmf_all")[shared_cell_cells, ]

# Unregressed: top 10 shape PCs + scaled log area concatenated
shape_mat_cell_unreg <- cbind(
  shape_pcs_cell,
  scaled_log_area = scaled_cell_log_area[shared_cell_cells]
)

# Regressed: regress log area + log depth + seg_method from shape and expression
cat("Regressing covariates (cell)...\n")
shape_mat_cell_reg <- .regress_lmer(
  mat              = shape_mat_cell_unreg,
  scaled_log_depth = scaled_log_depth_cell,
  seg_method       = seg_method_vec
)
expr_mat_cell_reg <- .regress_lmer(
  mat              = expr_mat_cell,
  scaled_log_depth = scaled_log_depth_cell,
  seg_method       = seg_method_vec
)

# ---- 8. Build shape and expression matrices (Nucleus) ----------------------

shape_pcs_nuc <- Seurat::Embeddings(obj, "tpca_nucleus")[shared_nuc_cells, 1:10]
expr_mat_nuc  <- Seurat::Embeddings(obj, "nmf_all")[shared_nuc_cells, ]

# Unregressed: top 10 shape PCs + scaled log area concatenated
shape_mat_nuc_unreg <- cbind(
  shape_pcs_nuc,
  scaled_log_area = scaled_nuc_log_area[shared_nuc_cells]
)

# Regressed
cat("Regressing covariates (nucleus)...\n")
shape_mat_nuc_reg <- .regress_lmer(
  mat              = shape_mat_nuc_unreg,
  scaled_log_depth = scaled_log_depth_nuc,
  seg_method       = seg_method_vec
)
expr_mat_nuc_reg <- .regress_lmer(
  mat              = expr_mat_nuc,
  scaled_log_depth = scaled_log_depth_nuc,
  seg_method       = seg_method_vec
)

# ---- 9. Run CCA — Cell, unregressed ----------------------------------------

cat("Running CCA without regression (cell)...\n")
cca_cell_unreg        <- link_shape_and_factors(
  obj       = obj,
  shape_mat = shape_mat_cell_unreg,
  expr_mat  = expr_mat_cell,
  group.by  = NULL
)
obj <- store_kstitch_results(obj, cca_cell_unreg,
                             reduction_key_csp  = "CSPCU_",
                             reduction_key_cep  = "CEPCU_")
obj@misc$kstitch$cell_cca_unregressed <- obj@misc$kstitch$all
obj@misc$kstitch$all <- NULL
cat(sprintf("CCA cell unregressed stored: %d components\n",
            length(cca_cell_unreg$CC_Corr_Coefs)))

# ---- 10. Run CCA — Cell, regressed -----------------------------------------

cat("Running CCA with regression (cell)...\n")
cca_cell_reg        <- link_shape_and_factors(
  obj       = obj,
  shape_mat = shape_mat_cell_reg,
  expr_mat  = expr_mat_cell_reg,
  group.by  = NULL
)
obj <- store_kstitch_results(obj, cca_cell_reg,
                             reduction_key_csp  = "CSPCR_",
                             reduction_key_cep  = "CEPCR_")
obj@misc$kstitch$cell_cca_regressed <- obj@misc$kstitch$all
obj@misc$kstitch$all <- NULL
cat(sprintf("CCA cell regressed stored: %d components\n",
            length(cca_cell_reg$CC_Corr_Coefs)))

# ---- 11. Run CCA — Nucleus, unregressed ------------------------------------

cat("Running CCA without regression (nucleus)...\n")
cca_nuc_unreg        <- link_shape_and_factors(
  obj       = obj,
  shape_mat = shape_mat_nuc_unreg,
  expr_mat  = expr_mat_nuc,
  group.by  = NULL
)
obj <- store_kstitch_results(obj, cca_nuc_unreg,
                             reduction_key_csp  = "CSPNU_",
                             reduction_key_cep  = "CEPNU_")
obj@misc$kstitch$nuc_cca_unregressed <- obj@misc$kstitch$all
obj@misc$kstitch$all <- NULL
cat(sprintf("CCA nucleus unregressed stored: %d components\n",
            length(cca_nuc_unreg$CC_Corr_Coefs)))

# ---- 12. Run CCA — Nucleus, regressed --------------------------------------

cat("Running CCA with regression (nucleus)...\n")
cca_nuc_reg        <- link_shape_and_factors(
  obj       = obj,
  shape_mat = shape_mat_nuc_reg,
  expr_mat  = expr_mat_nuc_reg,
  group.by  = NULL
)
obj <- store_kstitch_results(obj, cca_nuc_reg,
                             reduction_key_csp  = "CSPNR_",
                             reduction_key_cep  = "CEPNR_")

obj@misc$kstitch$nuc_cca_regressed <- obj@misc$kstitch$all
obj@misc$kstitch$all <- NULL

cat(sprintf("CCA nucleus regressed stored: %d components\n",
            length(cca_nuc_reg$CC_Corr_Coefs)))

# ---- 13. Trim and save -----------------------------------------------------

obj <- DietSeurat(obj,
                  assays    = Seurat::Assays(obj),
                  dimreducs = c("tpca_cell", "tpca_nucleus", "nmf_all",
                                "csp_cell_unreg", "cep_cell_unreg",
                                "csp_cell_reg",   "cep_cell_reg",
                                "csp_nuc_unreg",  "cep_nuc_unreg",
                                "csp_nuc_reg",    "cep_nuc_reg"))

out_dir  <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "kstitch_vignette_xenium.rds")

saveRDS(
  list(
    obj      = obj,
    metadata = list(
      tpca_cell    = tpca_cell$all$Metadata,
      tpca_nucleus = tpca_nucleus$all$Metadata
    )
  ),
  out_path
)

cat(sprintf("\nSaved to: %s\n", normalizePath(out_path)))
cat(sprintf("File size: %.1f MB\n", file.size(out_path) / 1e6))
