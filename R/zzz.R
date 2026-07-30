.onLoad <- function(libname, pkgname) {
  reticulate::py_require(c(
    "numpy", "scipy", "shapely", "pyarrow", "h5py", "multiprocess", "pandas"
  ))
}
