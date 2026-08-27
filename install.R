# install.R ----------------------------------------------------------------
#
# Installs the packages the analyses need.
#
#   Rscript install.R
#
# Bayesian estimation additionally needs CmdStan. If you only want to inspect
# the code, reproduce the design summary, or work from cached fits, you can
# skip the CmdStan step.

CRAN <- "https://cloud.r-project.org"

# Data handling, modelling, posterior summaries, and figures. Kept deliberately
# short: the analysis code depends on these and nothing else.
pkgs <- c(
  # data
  "here", "readr", "dplyr", "tidyr", "purrr", "tibble", "stringr",
  "forcats", "glue", "data.table",
  # modelling
  "lme4", "brms", "posterior", "bayestestR",
  # nonparametric decomposition
  "energy",
  # figures and tables
  "ggplot2", "ggrepel", "ggthemes", "scales", "kableExtra"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = CRAN)
} else {
  message("All CRAN packages already installed.")
}

# cmdstanr is not on CRAN.
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  message("Installing cmdstanr from the Stan repository...")
  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", CRAN))
}

if (requireNamespace("cmdstanr", quietly = TRUE)) {
  have_cmdstan <- tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE)

  if (have_cmdstan) {
    message("CmdStan found at ", cmdstanr::cmdstan_path(),
            " (version ", cmdstanr::cmdstan_version(), ").")
  } else {
    message(
      "CmdStan is not installed. Run cmdstanr::install_cmdstan() once ",
      "(a few minutes) before re-fitting any Bayesian model.\n",
      "This is only needed for REFIT=TRUE; cached fits load without it."
    )
  }
}

message("\nNext: Rscript tests/test_gtheory.R, then Rscript scripts/00_design_summary.R")
