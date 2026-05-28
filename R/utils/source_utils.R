# Load all shared utilities (call after R/setup.R).
utils_dir <- here::here("R", "utils")
for (f in list.files(utils_dir, pattern = "\\.R$", full.names = TRUE)) {
  if (!endsWith(f, "source_utils.R")) source(f, local = FALSE)
}
