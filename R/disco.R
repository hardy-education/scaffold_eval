# R/disco.R ----------------------------------------------------------------
#
# Nonparametric variance decomposition via distance components (DISCO).
#
# A distribution-free counterpart to the mixed-model decompositions, from the
# energy-statistics literature (Rizzo & Szekely 2010; Szekely & Rizzo 2017).
# DISCO splits total dispersion, measured through pairwise distances, into
# between- and within-group parts without assuming a link function, additive
# effects, or Gaussian random effects:
#
#   S_total = S_between + S_within,   Erho^2_DISCO = S_p / (S_p + S_within)
#
# It is included as a robustness check, not as a competing estimate. Two
# caveats from Appendix C.5 and G.3 matter when reading its output: it has no
# separate residual term, so unexplained variation is pushed into the named
# interactions, and it is sensitive to cluster imbalance. It attributes far
# more dispersion to task-model interaction than the parametric fits do, and is
# best read as an upper bound on interaction effects.
#
# Cost: roughly 7 hours pooled across eight benchmarks, 17 hours across nine.

# DISCO reports its facets by the column names it was given. Mapping them onto
# the component codes of the mixed-model registry is what lets the two
# decompositions be tabulated side by side in the estimator comparison.
DISCO_COMPONENT <- c(
  item              = "I",
  model             = "M",
  agent             = "A",
  benchmark         = "B",
  item_model        = "IM",
  item_agent        = "IA",
  model_agent       = "MA",
  bench_model       = "BM",
  bench_agent       = "BA",
  bench_model_agent = "BMA",
  item_model_agent  = "IMA"
)

#' Build the factor columns DISCO decomposes over.
#'
#' `energy::disco` takes a response and a data frame of factors, and decomposes
#' **sequentially** over the columns in the order given, so interaction facets
#' must be materialised as explicit factors and the column order is part of the
#' specification. Main effects come first, then two-way, then three-way, which
#' matches the order the components are entered in the mixed models.
#'
#' @param df Rollout data.
#' @param level "leaderboard" includes benchmark facets; "benchmark" assumes
#'   the data have already been subset to a single benchmark.
disco_factors <- function(df, level = c("leaderboard", "benchmark")) {
  level <- match.arg(level)

  item  <- factor(paste(df$benchmark, df$task_id, sep = ":"))
  model <- factor(df$model_name)
  agent <- factor(df$agent_name)
  bench <- factor(df$benchmark)

  out <- tibble(item = item, agent = agent, model = model)

  if (level == "leaderboard") out$benchmark <- bench

  out$item_agent  <- forcats::fct_drop(interaction(item, agent))
  out$item_model  <- forcats::fct_drop(interaction(item, model))
  out$model_agent <- forcats::fct_drop(interaction(model, agent))

  if (level == "leaderboard") {
    out$bench_agent       <- forcats::fct_drop(interaction(bench, agent))
    out$bench_model       <- forcats::fct_drop(interaction(bench, model))
    out$bench_model_agent <- forcats::fct_drop(interaction(bench, model, agent))
    out$item_model_agent  <- forcats::fct_drop(interaction(item, model, agent))
  }
  out
}


#' Run a DISCO decomposition.
#'
#' @param df Rollout data, already subset if `level = "benchmark"`.
#' @param level Which facet set to decompose over.
#' @param index Distance exponent alpha; 1 is the Euclidean default.
#' @return Tidy tibble with one row per facet, carrying the dispersion share
#'   (`pct`), the DISCO analogue of Erho^2, and the raw statistics.
disco_decompose <- function(df, level = c("leaderboard", "benchmark"), index = 1) {
  level <- match.arg(level)
  df <- filter(df, !is.na(score))

  resp    <- data.table::as.data.table(tibble(resp = as.numeric(df$score)))
  factors <- data.table::as.data.table(disco_factors(df, level))

  fit <- energy::disco(resp, factors = factors, distance = FALSE,
                       index = index, R = 0)

  tibble::as_tibble(fit$stats) |>
    select(-dplyr::any_of("p-value")) |>
    mutate(
      facet        = fit$factor.names,
      component    = unname(DISCO_COMPONENT[facet]),
      within_facet = Within,
      within_total = fit$within,
      total        = fit$total,
      ep2          = Trt / (Trt + within_facet),
      snr          = ep2 / (1 - ep2),
      pct          = Trt / (total - within_total),
      index        = index,
      .before = 1
    ) |>
    select(facet, component, ep2, snr, pct, everything())
}


#' DISCO decomposition within every benchmark.
disco_by_benchmark <- function(df, index = 1) {
  benches <- levels(droplevels(df$benchmark))
  purrr::map(benches, \(b) {
    message("disco: ", b)
    disco_decompose(filter(df, benchmark == b), level = "benchmark", index = index) |>
      mutate(benchmark = b, .before = 1)
  }) |>
    bind_rows()
}


#' Run DISCO, or read a cached result.
#'
#' DISCO runs for hours and its result is a small table, so it is cached as CSV
#' under outputs/results/ and re-read by default.
disco_or_cache <- function(df, name, level = "leaderboard", index = 1, refit = REFIT) {
  path <- file.path(PATHS$results, paste0(name, ".csv"))
  if (file.exists(path) && !refit) return(readr::read_csv(path, show_col_types = FALSE))

  out <- if (level == "leaderboard") {
    disco_decompose(df, level = "leaderboard", index = index)
  } else {
    disco_by_benchmark(df, index = index)
  }
  readr::write_csv(out, path)
  out
}
