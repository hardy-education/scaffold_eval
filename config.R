# config.R -----------------------------------------------------------------
#
# Central configuration for the leaderboard-reliability analyses.
#
# Everything a reader might reasonably want to change -- how long the samplers
# run, where outputs land, which designs the D-studies sweep over -- lives here
# rather than being scattered through the analysis scripts.
#
# Override any value from the shell before sourcing a script, e.g.
#
#   REFIT <- TRUE
#   source("scripts/01_fit_models.R")
#
# or from the command line:
#
#   REFIT=TRUE Rscript scripts/01_fit_models.R

#' Read a boolean flag from the environment.
#'
#' Accepts TRUE/true/T/1/yes and their negations, so `REFIT=1` and `REFIT=TRUE`
#' both work rather than one of them silently becoming NA.
env_flag <- function(name, default = FALSE) {
  raw <- trimws(tolower(Sys.getenv(name, "")))
  if (!nzchar(raw)) return(default)
  if (raw %in% c("1", "t", "true", "yes", "y")) return(TRUE)
  if (raw %in% c("0", "f", "false", "no", "n")) return(FALSE)
  stop("Cannot interpret ", name, "='", Sys.getenv(name), "' as TRUE or FALSE.",
       call. = FALSE)
}

# ---- Paths ----------------------------------------------------------------

PATHS <- list(
  data     = here::here("data", "hal_response_matrix.csv"),
  fits     = here::here("outputs", "fits"),
  figures  = here::here("outputs", "figures"),
  tables   = here::here("outputs", "tables"),
  results  = here::here("outputs", "results")
)

# ---- Estimation -----------------------------------------------------------

# REFIT = FALSE reuses cached fits in outputs/fits/ and errors if one is
# missing. REFIT = TRUE re-runs the sampler (hours; see README for timings).
REFIT <- env_flag("REFIT")

# Sampler settings, matched to Appendix C.1.3 of the paper.
SAMPLER <- list(
  # Pooled leaderboard fit (Eq. 6): 6 chains x 9,000 iter, 4,000 warmup,
  # thin 8 -> 3,750 retained draws. ~5.5 h on an Apple M1 Max.
  leaderboard = list(chains = 6, cores = 6, iter = 9000, warmup = 4000, thin = 8),
  # Per-benchmark fits (Eq. 2): ~6.5 min each.
  benchmark   = list(chains = 4, cores = 4, iter = 2000, warmup = 1000, thin = 3),
  # Leave-one-benchmark-out refits: 4 chains x 4,000 iter -> 1,600 draws.
  loo         = list(chains = 4, cores = 4, iter = 4000, warmup = 2000, thin = 5),
  adapt_delta = 0.95,
  threads     = 2,
  backend     = "cmdstanr",
  seed        = 20260101
)

# ---- Posterior summaries --------------------------------------------------

# Figures report medians with 0.68 HDI ("pseudo-standard-error") intervals,
# because several derived quantities have strongly skewed posteriors.
CI_LEVEL <- 0.68

# ---- D-study grids --------------------------------------------------------

DSTUDY <- list(
  # Total task budget swept in the leaderboard D-study (Figure 4).
  total_tasks     = 500,
  # Per-benchmark D-study task grid cap (Figures 2 and 3).
  bench_max_tasks = 300,
  # Benchmark counts to sweep.
  n_benchmarks    = 1:9,
  # Scaffold count assumed by the primary D-studies. The observed leaderboard
  # reports each system under a single scaffold, so na = 1 is the design that
  # current practice corresponds to.
  n_scaffolds     = 1
)

# Limit-of-detection reference band for signal-to-noise ratios
# (CDER 2024; Sheehan & Yost 2026). Used as a measurement-design diagnostic,
# not a pass/fail threshold.
LOD_BAND <- c(low = 2, mid = 2.5, high = 3)

# Reliability reference used when describing benchmarks as high- or low-signal.
RELIABILITY_TARGET <- 0.75

# ---- Presentation ---------------------------------------------------------

# Benchmarks featured in the two-panel figures of the main text.
FEATURED_BENCHMARKS <- c("taubench_airline", "swebench_verified_mini")

# Font family for figures. "serif" resolves to a Times-like face on every
# graphics device. Set a specific installed family (e.g. "Times New Roman") if
# you want to match a particular manuscript template.
PLOT_FONT <- "serif"  # R's device-independent Times-like serif
