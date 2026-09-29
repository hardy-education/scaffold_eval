# R/setup.R ----------------------------------------------------------------
#
# Single entry point for the analysis code. Sourcing this file loads packages,
# configuration, and every module in R/, in dependency order.
#
#   source(here::here("R", "setup.R"))
#
# It is idempotent: sourcing twice is harmless.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(glue)
})

source(here::here("config.R"))

# Modules, in dependency order. gtheory.R is the heart of the repository:
# every reliability number in the paper is produced by the engine it defines.
.modules <- c(
  "data.R",      # loading, validation, design descriptives
  "formulas.R",  # the random-effects structures of Eq. 2 and Eq. 6
  "fit.R",       # fit-or-load wrappers for brms and lme4
  "variance.R",  # variance components from a fit, as posterior draws
  "gtheory.R",   # G-coefficients, signal-to-noise, D-studies
  "ranks.R",     # latent model effects, posterior ranks, rank correlations
  "cost.R",      # design cost, reliability-aware allocation, task subsampling
  "harbor.R",    # the Harbor Index corroboration dataset
  "external.R",  # out-of-panel validation of the latent model effect
  "disco.R",     # nonparametric distance-components decomposition
  "plots.R"      # shared theme, palette, and the paper's figure builders
)

for (.m in .modules) source(here::here("R", .m))
rm(.m, .modules)

for (.p in PATHS[c("fits", "figures", "tables", "results")]) {
  if (!dir.exists(.p)) dir.create(.p, recursive = TRUE)
}
rm(.p)


#' Record the computing environment alongside the results.
#'
#' Package versions move, and variance-component estimates near a boundary can
#' be sensitive to optimiser and sampler versions. Writing the environment next
#' to the outputs makes a later discrepancy diagnosable rather than mysterious.
record_session <- function(path = file.path(PATHS$results, "session_info.txt")) {
  cmdstan <- tryCatch(
    paste0("CmdStan ", cmdstanr::cmdstan_version(), " at ", cmdstanr::cmdstan_path()),
    error = function(e) "CmdStan not installed"
  )
  writeLines(
    c(paste("Recorded:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      cmdstan, "",
      utils::capture.output(utils::sessionInfo())),
    path
  )
  message("wrote ", path)
  invisible(path)
}

invisible(TRUE)
