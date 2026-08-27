# R/ranks.R ----------------------------------------------------------------
#
# Latent capability estimates, posterior rankings, and their agreement with
# published leaderboards.
#
# A published leaderboard sorts models by mean accuracy, which credits a model
# for every condition it happened to be paired with. The decomposition instead
# estimates the part of performance that persists after averaging over the
# benchmark and scaffold universes. These functions extract that estimate,
# carry its uncertainty through to ranks, and compare the two orderings.
#
# Everything is computed per posterior draw: models are ranked within a draw
# and then summarised, so a reported interval is an interval on the rank, not a
# rank of intervals.

# ---- Extracting random effects safely --------------------------------------
#
# brms labels the levels of an interaction grouping factor by pasting the
# component variables with "_". Because model names themselves contain
# underscores, splitting those labels with a regular expression is ambiguous
# and quietly mis-assigns effects. Instead we build the expected labels from
# the observed design and join on them, then assert that every level matched.

#' Posterior draws of a fit, as `rvar`s, converting only if needed.
#'
#' Callers that need several grouping terms convert once and pass the result
#' down: the pooled fit is a few hundred megabytes and converting it repeatedly
#' is the slowest thing in this file.
as_rvar_draws <- function(x) {
  if (inherits(x, "draws_rvars")) x else posterior::as_draws_rvars(x)
}


#' Draws of the random effects for one grouping term, joined back to its keys.
#'
#' Reads the `r_<group>` parameters out of the posterior directly rather than
#' going through `brms::ranef()`, which dispatches through rstan's S4 machinery
#' even for a cmdstanr-backed fit. Going to the draws keeps this on the same
#' path the variance components already use, and one fewer package has to be
#' installed and working.
#'
#' @param fit A `brmsfit`, or draws already converted by `as_rvar_draws()`.
#' @param group Grouping term as it appears in the formula, e.g.
#'   `"benchmark:model_name"`.
#' @param keys The variables composing that term, in formula order.
#' @param data The data the model was fitted to, used to reconstruct labels.
#' @param value Name for the effect column.
#' @return A tibble with one row per (draw, level) and the key columns restored.
ranef_draws <- function(fit, group, keys, data, value = "effect") {
  draws <- as_rvar_draws(fit)
  parameter <- paste0("r_", group)

  if (is.null(draws[[parameter]])) {
    stop("The fit has no random effect `", parameter, "`. Available: ",
         paste(grep("^r_", names(draws), value = TRUE), collapse = ", "),
         call. = FALSE)
  }

  # [draw, level, coefficient]; these models have intercepts only.
  arr <- posterior::draws_of(draws[[parameter]])
  levels_observed <- dimnames(arr)[[2]]

  lookup <- data |>
    distinct(across(all_of(keys))) |>
    mutate(.level = do.call(paste, c(across(all_of(keys)), sep = "_")))

  unmatched <- setdiff(levels_observed, lookup$.level)
  if (length(unmatched)) {
    stop(
      "Could not map ", length(unmatched), " level(s) of `", group,
      "` back to the data (e.g. '", unmatched[1], "'). ",
      "This usually means `keys` does not match the grouping term.",
      call. = FALSE
    )
  }

  # Rebuilding the matrix explicitly keeps it two-dimensional even when a
  # grouping term happens to have a single level.
  mat <- matrix(arr[, , "Intercept"], nrow = dim(arr)[1],
                dimnames = list(NULL, levels_observed))

  tibble::as_tibble(mat) |>
    mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(-.draw, names_to = ".level", values_to = value) |>
    inner_join(lookup, by = ".level") |>
    select(-.level)
}


#' Posterior draws of scaffold-marginalized latent capability.
#'
#' For benchmark b and model m the paper defines
#'
#'   theta_mb = u_m^(M) + u_bm^(BM)
#'
#' the model's globally persistent capability plus its benchmark-specific
#' departure, with scaffold effects marginalized away. This is the quantity
#' whose ranking is compared against the published leaderboard.
#'
#' @param fit The pooled leaderboard fit.
#' @param df The rollout data.
#' @return Tibble of `.draw`, `benchmark`, `model_name`, `theta`.
latent_capability_draws <- function(fit, df) {
  draws <- as_rvar_draws(fit)
  model_effect <- ranef_draws(draws, "model_name", "model_name", df, "u_model")
  bench_model  <- ranef_draws(draws, "benchmark:model_name",
                              c("benchmark", "model_name"), df, "u_bench_model")

  bench_model |>
    inner_join(model_effect, by = c(".draw", "model_name")) |>
    mutate(theta = u_model + u_bench_model) |>
    select(.draw, benchmark, model_name, u_model, u_bench_model, theta)
}


#' Posterior draws of latent capability for model-scaffold systems.
#'
#' Adds the benchmark-specific model-scaffold compatibility term, giving the
#' latent score of the deployable system rather than of the underlying model.
latent_system_draws <- function(fit, df) {
  draws <- as_rvar_draws(fit)
  latent_capability_draws(draws, df) |>
    inner_join(
      ranef_draws(draws, "benchmark:model_name:agent_name",
                  c("benchmark", "model_name", "agent_name"), df, "u_bma"),
      by = c(".draw", "benchmark", "model_name")
    ) |>
    mutate(theta = theta + u_bma) |>
    select(.draw, benchmark, model_name, agent_name, theta)
}


#' Latent capability from the independently fitted per-benchmark models.
#'
#' The comparison estimator of Appendix D.3: uses only the information inside
#' a single benchmark, with no cross-benchmark pooling.
latent_capability_per_benchmark <- function(fits, df) {
  purrr::imap(fits, \(fit, b) {
    ranef_draws(as_rvar_draws(fit), "model_name", "model_name",
                filter(df, benchmark == b), "theta") |>
      mutate(benchmark = b, .before = 1)
  }) |>
    bind_rows()
}


# ---- Posterior ranks -------------------------------------------------------

#' Summarise posterior rank distributions.
#'
#' Ranking happens inside each draw, so the reported interval reflects how much
#' the ordering itself moves under the posterior.
#'
#' Two conventions, each matching the thing it is compared against:
#'
#' - integer ranks count from the top, so **1 is best**, like `published_rank`;
#' - percentile ranks increase with capability, so **1 is best and 0 is worst**,
#'   like `published_pctrank` from `published_scores()`.
#'
#' Getting these the same way round matters: Figure 6 plots one against the
#' other with a diagonal reference line.
#'
#' @param draws Output of `latent_capability_draws()` (or the system variant).
#' @param keys Columns identifying the ranked object.
#' @param percentile Report percentile ranks in [0, 1] instead of integer ranks,
#'   which makes benchmarks with different model counts comparable on one axis.
posterior_ranks <- function(draws, keys = c("model_name"), percentile = FALSE) {
  rank_fun <- if (percentile) {
    \(x) dplyr::percent_rank(x)
  } else {
    \(x) rank(-x, ties.method = "average")
  }

  draws |>
    group_by(benchmark, .draw) |>
    mutate(rank = rank_fun(theta)) |>
    ungroup() |>
    group_by(across(all_of(c("benchmark", keys)))) |>
    summarize(
      theta_median = stats::median(theta),
      rank_median  = stats::median(rank),
      rank_mean    = mean(rank),
      rank_q2_5    = stats::quantile(rank, 0.025, names = FALSE),
      rank_q97_5   = stats::quantile(rank, 0.975, names = FALSE),
      rank_q5      = stats::quantile(rank, 0.05, names = FALSE),
      rank_q95     = stats::quantile(rank, 0.95, names = FALSE),
      .groups = "drop"
    )
}


#' Pairwise posterior ordering probabilities, Pr(theta_m > theta_m' | y, b).
#'
#' The basis for calling two published ranks indistinguishable: a pair whose
#' probability sits near 0.5 is a tie the data do not resolve.
#'
#' @param threshold Probabilities within `[1 - threshold, threshold]` are
#'   counted as unresolved.
pairwise_ordering <- function(draws, threshold = 0.95) {
  draws |>
    group_by(benchmark) |>
    group_modify(function(d, key) {
      wide <- d |>
        select(.draw, model_name, theta) |>
        tidyr::pivot_wider(names_from = model_name, values_from = theta)
      mat <- as.matrix(wide[, -1, drop = FALSE])
      models <- colnames(mat)

      pairs <- utils::combn(seq_along(models), 2)
      tibble(
        model_a = models[pairs[1, ]],
        model_b = models[pairs[2, ]],
        p_a_beats_b = apply(pairs, 2, \(ij) mean(mat[, ij[1]] > mat[, ij[2]], na.rm = TRUE))
      )
    }) |>
    ungroup() |>
    mutate(resolved = p_a_beats_b > threshold | p_a_beats_b < (1 - threshold))
}


#' Share of model pairs the posterior cannot order (paper Section 5.6).
indistinguishable_share <- function(draws, threshold = 0.95) {
  pairwise_ordering(draws, threshold) |>
    group_by(benchmark) |>
    summarize(
      n_pairs             = n(),
      n_resolved          = sum(resolved),
      prop_indistinguishable = 1 - mean(resolved),
      .groups = "drop"
    )
}


# ---- Agreement with published rankings -------------------------------------

#' Rank-agreement statistics between two orderings.
#'
#' Spearman's rho and Kendall's tau measure monotone and pairwise agreement.
#' The bias-corrected squared distance correlation (dCor^2_n) is reported
#' alongside them because it behaves sensibly under the extensive ties that
#' low-scoring benchmarks produce.
rank_agreement <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) {
    return(tibble(spearman = NA_real_, kendall = NA_real_, dcor2n = NA_real_))
  }
  tibble(
    spearman = stats::cor(x, y, method = "spearman"),
    kendall  = stats::cor(x, y, method = "kendall"),
    dcor2n   = energy::bcdcor(rank(x), rank(y))
  )
}


#' Agreement between published and latent rankings, per posterior draw.
#'
#' Computing the correlation inside each draw and then summarising -- rather
#' than correlating posterior medians -- keeps the uncertainty in the reported
#' agreement, which is the point of Table 2.
rank_agreement_draws <- function(draws, published, ci = c(0.05, 0.95)) {
  draws |>
    inner_join(select(published, benchmark, model_name, score),
               by = c("benchmark", "model_name")) |>
    group_by(benchmark, .draw) |>
    group_modify(\(d, key) rank_agreement(d$theta, d$score)) |>
    ungroup() |>
    tidyr::pivot_longer(c(spearman, kendall, dcor2n),
                        names_to = "metric", values_to = "value") |>
    group_by(benchmark, metric) |>
    summarize(
      median = stats::median(value, na.rm = TRUE),
      mean   = mean(value, na.rm = TRUE),
      low    = stats::quantile(value, ci[1], na.rm = TRUE, names = FALSE),
      high   = stats::quantile(value, ci[2], na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}


#' Rank shift table: published rank against posterior median rank.
#'
#' `delta` is positive when accounting for task, scaffold, and benchmark
#' variation moves a model down the leaderboard.
rank_shifts <- function(ranks, published, keys = "model_name") {
  ranks |>
    inner_join(published, by = c("benchmark", keys)) |>
    mutate(delta = rank_median - published_rank) |>
    arrange(benchmark, rank_median)
}
