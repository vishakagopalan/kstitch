# ---- .serialize_group_results() --------------------------------------------

# Serialize a named list of per-group results to RDS files under output_dir.
# Returns a named character vector of file paths (same keys as `results`).
# `type` is a short string used only in the console message ("nmf" or "cca").
.serialize_group_results <- function(results, output_dir, type = "nmf") {
  paths <- stats::setNames(
    file.path(output_dir, paste0(names(results), ".rds")),
    names(results)
  )
  for (grp in names(results)) {
    saveRDS(results[[grp]], file = paths[[grp]])
  }
  message("Results not loaded into memory (return_results = FALSE).")
  message(sprintf("  %s results written to: %s", toupper(type), output_dir))
  message(sprintf("  Load individual groups with load_kstitch_results(path, type = \"%s\")", type))
  paths
}
