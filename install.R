# Install required packages for scaffold_eval reliability analyses.
pkgs <- c(
  "data.table",
  "dtplyr",
  "here",
  "glue",
  "tictoc",
  "tidyverse",
  "forcats",
  "purrr",
  "lme4",
  "lmerTest",
  "janitor",
  "broom",
  "broom.mixed",
  "brms",
  "cmdstanr",
  "bayestestR",
  "posterior",
  "energy",
  "ggthemes",
  "ggrepel",
  "ggpubr",
  "kableExtra"
)

cran <- setdiff(pkgs, c("cmdstanr"))
install.packages(cran, repos = "https://cloud.r-project.org")

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages(
    "cmdstanr",
    repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
  )
}

message(
  "Packages installed. For Bayesian models, run cmdstanr::install_cmdstan() once if CmdStan is not on your PATH."
)
