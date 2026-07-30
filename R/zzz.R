.onLoad <- function(libname, pkgname) {
  pkgs <- c("numpy", "scipy", "shapely", "pyarrow", "h5py", "multiprocess", "pandas")
  missing <- pkgs[!sapply(pkgs, reticulate::py_module_available)]
  if (length(missing) > 0) {
    reticulate::py_install(missing)
  }
}
