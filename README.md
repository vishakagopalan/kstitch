# kstitch <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/vishakagopalan/kstitch/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vishakagopalan/kstitch/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/kstitch)](https://CRAN.R-project.org/package=kstitch)
<!-- badges: end -->

**kstitch** integrates cell and nucleus **morphology** with **gene expression**
in spatial transcriptomics data. It works with Seurat v5 objects from Xenium
and CosMx platforms, stitching together shape and molecular information via
Tangent PCA (TPCA) and canonical correlation analysis (CCA).

---

## Overview

Spatial transcriptomics captures where genes are expressed, but cell
morphology — how a cell looks — carries independent biological information.
kstitch links these two modalities by:

1. **TPCA** — computing principal geodesic components of cell/nucleus contour
   polygons in Kendall shape space (via a Python backend, `PGA.py`).
2. **NMF** — factorising gene expression counts into interpretable latent
   factors via `rliger`.
3. **CCA** — finding linear combinations of shape PCs and NMF factors that are
   maximally correlated, revealing axes of joint morpho-molecular variation.

Results live inside the Seurat object and are visualised with four standard
figure functions.

---

## Installation

```r
# Install from GitHub (requires remotes)
remotes::install_github("vishakagopalan/kstitch")
```

kstitch requires Python for the TPCA step. Python dependencies are provisioned
automatically via `reticulate` (using an ephemeral `uv` virtualenv) when the
package is loaded — no manual Python setup is needed.

```r
library(kstitch)
```

---

## Quick start

```r
library(kstitch)

# --- 1. Run TPCA on boundary parquet files ---
tpca_cell <- run_tpca(
  boundary_parquet_path = "cell_boundaries.parquet",
  output_dir            = "output/tpca_cell",
  contour_type          = "cell"
)
obj <- store_tpca_results(obj, tpca_cell)

# --- 2. Run NMF ---
nmf_results <- compute_nmf(obj, assay_name = "Xenium", group.by = "cell_type")

# --- 3. Link shape and expression via CCA ---
cca_results <- link_shape_and_factors(
  obj        = obj,
  expr_mat    = nmf_results$Keratinocytes$NMF_Matrix,
  shape_mat  = Seurat::Embeddings(obj, "tpca_cell"),
  group.by   = "cell_type"
)
obj <- store_kstitch_results(obj, cca_results)

# --- 4. Visualise ---
grp <- get_kstitch_results(obj, "Keratinocytes")

plot_shape_modes(tpca_cell)
plot_csp_loadings(grp, cc_idx = 1)
plot_cep_loadings(grp, cc_idx = 1)
plot_csp_boundary_montage(grp, tpca_cell, cc_idx = 1)
```

For a full walkthrough see `vignette("kstitch-xenium")`.

---

## Function reference

| Family | Functions |
|--------|-----------|
| **TPCA** | `run_tpca()`, `store_tpca_results()`, `get_tpca_results()` |
| **NMF** | `compute_nmf()` |
| **CCA** | `run_cca()`, `cca_pvalues()`, `link_shape_and_factors()` |
| **CCA helpers** | `anchor_cca_signs()`, `prioritize_cca_components()` |
| **Storage** | `store_kstitch_results()`, `get_kstitch_results()` |
| **Figures** | `plot_shape_modes()`, `plot_csp_loadings()`, `plot_cep_loadings()`, `plot_csp_boundary_montage()` |

---

## Citation

If you use kstitch in your research, please cite:

> Gopalan V. et al. (2025). kstitch: integrating cell morphology and gene
> expression in spatial transcriptomics. *In preparation.*

---

## License

MIT © Vishaka Gopalan
