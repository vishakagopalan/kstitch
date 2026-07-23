library(Seurat)
library(kstitch)
library(tidyr)

xenium_data_path <- "../../data/Xenium_Prime_Human_Skin_FFPE_outs/"
precomputed_obj <- readRDS(system.file("extdata", "kstitch_vignette_xenium.rds",
                                              package = "kstitch"))
precomputed_nuc_tpca_mat <- Embeddings( precomputed_obj, "tpca_nucleus")
keratinocyte_cell_ids <- Cells(precomputed_obj)
keratinocytes_cells_with_single_nucleus <- rownames(na.omit(precomputed_nuc_tpca_mat))
keratinocyte_genes_to_use <- Features(precomputed_obj) #These are genes whose expression is >= 10 nTPM in RNA-seq of keratinocytes from Human Protein Atlas

melanoma_xenium_obj <- LoadXenium(xenium_data_path,molecule.coordinates=F)
keratinocyte_obj <- subset( melanoma_xenium_obj, cells = keratinocyte_cell_ids,
                            features=keratinocyte_genes_to_use) #You may receive warning messages about object validation. These can be ignored for this vignette.

nmf_output <- compute_nmf( keratinocyte_obj, assay_name = "Xenium", stability_seed = 235, default_num_factors = 10,
                           num_top_genes_per_factor=50,
                           integration_method="online",
                           selectGenes_params=list(nGenes="all",verbose=T), gene_pct_threshold = 0.01)

melanoma_xenium_segmentation_info <- ReadXenium(xenium_data_path, outs=c("segmentation_method"),
                                                type=c("segmentations","nucleus_segmentations","centroids"))

tpca_out <- export_seurat_contours(seg_list=melanoma_xenium_segmentation_info,
                                   contour_type = "nucleus",
                                   cell_ids = keratinocytes_cells_with_single_nucleus)

num_shape_pcs_to_use <- 10
tpca_mat <- tpca_out$TPCA_Embedding[,1:num_shape_pcs_to_use]

area_mat <- tpca_out$Metadata[, "area", drop = FALSE]
rownames(area_mat) <- tpca_out$Metadata[["cell_id"]]
area_mat <- as.matrix(area_mat)

common_ids <- intersect( rownames(area_mat), rownames(tpca_mat))
tpca_with_area_mat <- cbind( tpca_mat[common_ids,,drop=F], area_mat[common_ids,,drop=F] )

cca_out <- link_shape_and_factors(keratinocyte_obj, expr_mat=nmf_output$all$NMF_Matrix,
                                  shape_mat = tpca_with_area_mat, test_significance = F, verbose = T)

