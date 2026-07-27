# ---- .serialize_result() ---------------------------------------------------

# Serialize a single result to an RDS file under output_dir.
# Returns the file path. `type` is used only in the console message.
.serialize_result <- function(result, output_dir, type = "nmf") {
  path <- file.path(output_dir, "all.rds")
  saveRDS(result, file = path)
  message(sprintf("Result not loaded into memory (return_results = FALSE)."))
  message(sprintf("  %s result written to: %s", toupper(type), path))
  message(sprintf(
    "  Load with load_kstitch_results(\"%s\", type = \"%s\")", path, type
  ))
  path
}
