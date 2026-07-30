# kstitch 

<!-- badges: start -->
[![R-CMD-check](https://github.com/vishakagopalan/kstitch/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vishakagopalan/kstitch/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**kstitch** is a framework for identifying interpretable covariation between cell or nuclear morphology and gene expression in imaging-based spatial transcriptomics data.

The package is designed to work with Seurat v5 objects and boundary data from platforms such as Xenium and CosMx. It represents contour shape using tangent principal component analysis (TPCA), models gene-expression programs using non-negative matrix factorization (NMF), and links the two modalities using canonical correlation analysis (CCA).

## Get started

See the **[complete Xenium workflow](https://vishakagopalan.github.io/kstitch/articles/kstitch-xenium.html)** for installation details, data download and extraction paths, NMF, cell and nuclear TPCA, covariate regression, CCA, and visualization.

After installing the package, the same vignette can be opened from R:

```r
vignette("kstitch-xenium", package = "kstitch")
```

The workflow is:

```text
Xenium output directory
        ↓
Load expression and boundary data
        ↓
NMF expression factors
        +
TPCA contour features and area
        ↓
Optional covariate regression
        ↓
CCA
        ↓
CSP/CEP loadings and boundary montages
```

---

## Overview

Imaging-based spatial transcriptomics provides paired measurements of gene expression and cell or nuclear boundaries. kstitch analyzes these data in three main steps:

1. **TPCA** — represents cell or nuclear boundary contours in Kendall shape space after removing translation, rotation, reflection, cyclic starting-point differences, and scale. TPCA then provides a low-dimensional Euclidean approximation of the dominant modes of contour-shape variation.
2. **NMF** — factorizes gene-expression count matrices into interpretable latent expression programs using `rliger`.
3. **CCA** — identifies linear combinations of area-augmented TPCA features and NMF factors whose cell-level projection scores are maximally correlated.

Before CCA, kstitch can regress technical covariates from the morphology and expression representations. Depending on the dataset, these may include total RNA count, distance from the field-of-view edge, sample or slide identity, and segmentation-related covariates.

The resulting canonical shape projections (CSPs), canonical expression projections (CEPs), and loading vectors can be stored in a Seurat object and visualized with the package plotting functions.

---

## Installation

First, install `remotes` if you don't have it:

```r
install.packages("remotes")
```

Then install kstitch from GitHub:

```r
remotes::install_github("vishakagopalan/kstitch")
```

kstitch depends on `rhdf5` and `HDF5Array`, which are available from Bioconductor and must be installed separately:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("rhdf5", "HDF5Array"))
```

kstitch uses a Python backend for TPCA. Python dependencies are managed through `reticulate`, which is installed automatically with kstitch. Required Python packages are declared via `py_require()` and will be installed automatically when Python is first initialized. In most cases no manual action is needed.

If you encounter issues with Python dependencies, you can install them manually:

```r
reticulate::virtualenv_install("r-reticulate",
  packages = c(
    "h5py>=3.16.0",
    "numpy>=2.5.1",
    "pandas>=3.0.3",
    "pyarrow>=25.0.0",
    "scipy>=1.18.0",
    "shapely>=2.1.2",
    "multiprocess>=0.70.19"
  )
)
```

**macOS users** may also need to install [Pandoc](https://pandoc.org/installing.html) and [XQuartz](https://www.xquartz.org/) before building vignettes.

To access the vignette:

```r
remotes::install_github(
  "vishakagopalan/kstitch",
  build_vignettes = TRUE,
  force = TRUE,
  upgrade = "never"
)
vignette("kstitch-xenium", package = "kstitch")
```

Alternatively, the vignette is always available at the [kstitch documentation site](https://vishakagopalan.github.io/kstitch/articles/kstitch-xenium.html) without any additional installation steps.

---

## Main outputs

For each analyzed cell population, kstitch can return:

- **Shape PC scores** describing size-independent contour variation.
- **Cell or nuclear area** stored separately from shape.
- **NMF factor scores and gene weights** describing expression programs.
- **Canonical shape projections (CSPs)** and their morphology coefficients.
- **Canonical expression projections (CEPs)** and their expression-factor coefficients.
- **Canonical correlations** describing the strength of association between paired CSPs and CEPs.

---

## Function reference

| Family | Functions |
|---|---|
| **TPCA** | `run_tpca()`, `run_tpca_from_seurat()`, `store_tpca_results()`, `get_tpca_results()`, `load_kstitch_results()` |
| **NMF** | `compute_nmf()`, `store_nmf_results()`, `get_nmf_results()` |
| **CCA** | `run_cca()`, `cca_pvalues()`, `link_shape_and_factors()`, `anchor_cca_signs()` |
| **Storage** | `store_kstitch_results()`, `get_kstitch_results()` |
| **Figures** | `plot_shape_modes()`, `plot_mu_history()`, `plot_frechet_convergence()`, `plot_csp_loadings()`, `plot_cep_loadings()`, `plot_csp_boundary_montage()` |

---

## Software authors

**Vishaka Gopalan** and **Shashwat Kumar** are co-authors of the kstitch software.

---

## Citation

If you use kstitch in your research, please cite:

> Kumar S, Shi Y, Vallius T, Day C-P, Absil P-A, Srivastava A, Hannenhalli S, Gopalan V. **kstitch links cellular morphology and gene expression in spatial transcriptomics.** bioRxiv (2026). https://doi.org/10.64898/2026.06.07.730714

---

## License

MIT © Vishaka Gopalan and Shashwat Kumar
