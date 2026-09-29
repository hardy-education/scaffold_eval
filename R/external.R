# R/external.R -------------------------------------------------------------
#
# Out-of-panel validation of the latent model effect.
#
# The pooled fit estimates a shared model component, theta_m, meant to capture
# capability that persists across benchmarks rather than advantage specific to
# one. That is a claim about transportability, and it can be tested: if
# theta_m is real signal rather than an artefact of these nine benchmarks, it
# should predict performance on benchmarks the model never saw during
# estimation, and should do so better than the conventional leaderboard
# aggregate (the unweighted mean score across benchmarks).
#
# Four contemporaneous external benchmarks are used, each reporting per-model
# scores and overlapping the HAL model pool by at least six models.
#
# Contamination control. Two external benchmarks overlap a HAL benchmark in
# content: SWE-bench Verified is a superset of SWE-bench Verified Mini, and
# tau^2-bench Core contains tau-bench Airline. Correlating theta_m against
# those while the corresponding HAL benchmark was in the fit would credit the
# estimate for information it was handed. For those two comparisons both
# estimators are therefore rebuilt with the overlapping benchmark removed:
# theta_m from the matching leave-one-benchmark-out refit, and the mean-score
# baseline from the remaining eight benchmarks. `EXTERNAL_EXCLUSIONS` records
# the mapping, and `external_validation()` refuses to proceed if a required
# refit is missing rather than silently using the full-data fit.
#
# What this does and does not show. Agreement here is convergent predictive
# evidence: the estimate carries information that transports to evaluations
# outside the panel. It is not evidence of construct validity, nor that any
# single latent dimension captures agentic capability, nor that either
# ranking is correct.

#' External benchmarks, and the HAL benchmark each one overlaps.
#'
#' `NA` means no content overlap, so the full-data fit is used directly.
EXTERNAL_EXCLUSIONS <- c(
  bfcl               = NA,
  terminal_bench_2_0 = NA,
  swebench           = "swebench_verified_mini",
  tau_2_core         = "taubench_airline"
)

EXTERNAL_LABELS <- c(
  bfcl               = "BFCL v4",
  terminal_bench_2_0 = "Terminal-Bench 2.0",
  swebench           = "SWE-bench Verified",
  tau_2_core         = "tau^2-bench Core"
)


#' Load the external validation scores.
#'
#' Several external leaderboards report the same benchmark, so scores are
#' averaged per (benchmark, model) before use. The source columns are kept so
#' provenance stays inspectable.
load_external_data <- function(path = PATHS$external) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE) |>
    mutate(score = as.numeric(score)) |>
    filter(!is.na(score))
}


#' Mean external score per (benchmark, model).
external_scores <- function(edf) {
  edf |>
    group_by(benchmark, model_name) |>
    summarize(external_score = mean(score), n_sources = n(), .groups = "drop")
}


#' Posterior-median latent model effect from a pooled fit.
#'
#' theta_m is the model main effect: the component expected to persist across
#' benchmarks, tasks, and scaffolds. Note this is the leaderboard-level
#' `u_m^(M)` alone, without the benchmark-conditioned departure, because the
#' claim under test is exactly that the shared part transports.
latent_model_effect <- function(fit, df) {
  ranef_draws(as_rvar_draws(fit), "model_name", "model_name", df, "u_model") |>
    group_by(model_name) |>
    summarize(theta = stats::median(u_model), .groups = "drop")
}


#' The conventional leaderboard aggregate: unweighted mean across benchmarks.
#'
#' Benchmarks are weighted equally rather than by task count, which is what a
#' leaderboard reporting one row per model effectively does.
mean_score_aggregate <- function(df, exclude = NULL) {
  if (!is.null(exclude) && !is.na(exclude)) {
    df <- filter(df, benchmark != exclude)
  }
  df |>
    group_by(benchmark, model_name) |>
    summarize(score = mean(score), .groups = "drop") |>
    group_by(model_name) |>
    summarize(mean_score = mean(score), .groups = "drop")
}


#' Assemble the per-model comparison table for every external benchmark.
#'
#' For each external benchmark, pairs its score with the two competing
#' in-panel estimators, using the contamination-controlled versions where the
#' content overlaps.
#'
#' @param df HAL rollout data.
#' @param fit The pooled leaderboard fit.
#' @param edf External validation data.
#' @param loo_fits Named list of leave-one-benchmark-out fits, keyed by the
#'   benchmark removed. Required for the two overlapping comparisons.
#' @return One row per (external benchmark, model) with `external_score`,
#'   `mean_score`, and `theta`.
external_comparison_data <- function(df, fit, edf, loo_fits = list()) {
  ext <- external_scores(edf)

  # Converting a pooled fit to draws is the expensive step, so the full-data
  # effect is computed once and reused by every non-overlapping comparison.
  theta_full <- latent_model_effect(fit, df)
  base_full  <- mean_score_aggregate(df)

  purrr::map(names(EXTERNAL_EXCLUSIONS), function(b) {
    excluded <- EXTERNAL_EXCLUSIONS[[b]]
    overlaps <- !is.na(excluded)

    if (overlaps) {
      if (is.null(loo_fits[[excluded]])) {
        stop(
          "External benchmark '", b, "' overlaps HAL benchmark '", excluded,
          "', so it needs the leave-one-benchmark-out refit without it.\n",
          "Estimate it with either of:\n",
          "  REFIT=TRUE Rscript scripts/08_external_validation.R   (just the two needed)\n",
          "  RUN_LOBO_BAYES=TRUE REFIT=TRUE Rscript scripts/05_ablations.R  (all nine)\n",
          "Using the full-data fit here would let the estimate borrow ",
          "information from the benchmark it is being validated against.",
          call. = FALSE
        )
      }
      panel <- filter(df, benchmark != excluded)
      theta <- latent_model_effect(loo_fits[[excluded]], panel)
      base  <- mean_score_aggregate(panel)
    } else {
      theta <- theta_full
      base  <- base_full
    }

    ext |>
      filter(benchmark == b) |>
      inner_join(base, by = "model_name") |>
      inner_join(theta, by = "model_name") |>
      mutate(excluded_from_panel = excluded)
  }) |>
    bind_rows() |>
    mutate(benchmark_label = EXTERNAL_LABELS[as.character(benchmark)])
}


#' Kendall's tau of each in-panel estimator against each external benchmark.
#'
#' Reproduces the paper's external-agreement table. `stats::cor.test` reports
#' no interval for Kendall's tau, so intervals are bootstrapped over models --
#' the same resampling scheme `delta_tau_bootstrap()` uses for the difference,
#' which keeps the two tables mutually consistent.
#'
#' With six to twenty overlapping models these intervals are wide, and should
#' be read as such: the comparison between estimators (the paired difference)
#' is far better determined than either tau on its own.
external_agreement <- function(comparison, n_boot = 20000L,
                               conf_level = 0.95, seed = SAMPLER$seed) {
  set.seed(seed)
  alpha <- 1 - conf_level

  safe_tau <- function(x, y) {
    if (length(x) < 2L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
      return(NA_real_)
    }
    unname(stats::cor(x, y, method = "kendall"))
  }

  tau_ci <- function(x, y) {
    ok <- stats::complete.cases(x, y)
    x <- x[ok]; y <- y[ok]; n <- length(x)
    if (n < 3) return(tibble(tau = NA_real_, low = NA_real_, high = NA_real_))
    boot <- replicate(n_boot, {
      i <- sample.int(n, n, replace = TRUE)
      safe_tau(x[i], y[i])
    })
    boot <- boot[is.finite(boot)]
    tibble(
      tau  = safe_tau(x, y),
      low  = if (length(boot) > 1) stats::quantile(boot, alpha / 2, names = FALSE) else NA_real_,
      high = if (length(boot) > 1) stats::quantile(boot, 1 - alpha / 2, names = FALSE) else NA_real_
    )
  }

  comparison |>
    group_by(benchmark, benchmark_label, excluded_from_panel) |>
    group_modify(function(d, key) {
      bind_rows(
        tau_ci(d$mean_score, d$external_score) |> mutate(estimator = "Mean score"),
        tau_ci(d$theta,      d$external_score) |> mutate(estimator = "theta_hat")
      ) |>
        mutate(n_models = nrow(d))
    }) |>
    ungroup() |>
    select(benchmark, benchmark_label, excluded_from_panel, estimator,
           tau, low, high, n_models)
}


#' Paired bootstrap of the difference in Kendall's tau.
#'
#' The question is comparative -- does theta_hat agree with the external
#' benchmark *more than* the mean score does -- so the uncertainty that matters
#' is on the difference, computed on the same resampled models rather than
#' from two independent intervals.
#'
#' Also reports leave-one-model-out deltas, because with six to twenty
#' overlapping models a single model can carry the comparison. A sign reversal
#' under deletion is a warning the difference is not robust.
#'
#' @param comparison Output of `external_comparison_data()`.
#' @param n_boot Bootstrap replicates.
#' @param conf_level Interval mass.
#' @param seed RNG seed.
delta_tau_bootstrap <- function(comparison, n_boot = 20000L,
                                conf_level = 0.95, seed = SAMPLER$seed) {
  stopifnot(n_boot >= 2)
  set.seed(seed)

  safe_tau <- function(x, y) {
    if (length(x) < 2L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
      return(NA_real_)
    }
    unname(stats::cor(x, y, method = "kendall"))
  }
  deltas <- function(d) {
    tau_mean  <- safe_tau(d$mean_score, d$external_score)
    tau_theta <- safe_tau(d$theta,      d$external_score)
    c(tau_mean = tau_mean, tau_theta = tau_theta,
      delta = tau_theta - tau_mean)
  }

  alpha <- 1 - conf_level

  comparison |>
    group_by(benchmark, benchmark_label) |>
    group_modify(function(d, key) {
      d <- d[stats::complete.cases(d$mean_score, d$theta, d$external_score), ]
      n <- nrow(d)
      if (n < 3L) {
        stop("External benchmark '", key$benchmark, "' has ", n,
             " complete observations; at least 3 are needed.", call. = FALSE)
      }

      observed <- deltas(d)

      boot <- replicate(n_boot, {
        deltas(d[sample.int(n, n, replace = TRUE), , drop = FALSE])[["delta"]]
      })
      valid <- boot[is.finite(boot)]

      # Leave-one-model-out: how much does any single model move the answer?
      loo <- vapply(seq_len(n), \(i) deltas(d[-i, , drop = FALSE])[["delta"]],
                    numeric(1))
      loo <- loo[is.finite(loo)]

      tibble(
        n_models        = n,
        tau_mean_score  = observed[["tau_mean"]],
        tau_theta       = observed[["tau_theta"]],
        delta           = observed[["delta"]],
        boot_low        = if (length(valid) > 1) stats::quantile(valid, alpha / 2, names = FALSE) else NA_real_,
        boot_high       = if (length(valid) > 1) stats::quantile(valid, 1 - alpha / 2, names = FALSE) else NA_real_,
        boot_undefined  = (n_boot - length(valid)) / n_boot,
        p_delta_gt_0    = mean(valid > 0),
        loo_min_delta   = if (length(loo)) min(loo) else NA_real_,
        loo_max_delta   = if (length(loo)) max(loo) else NA_real_,
        loo_sign_flips  = sum(sign(loo) != sign(observed[["delta"]]))
      )
    }) |>
    ungroup()
}


#' Kendall's tau between every pair of in-panel rankings.
#'
#' Context for the external comparison: benchmarks within the panel agree with
#' each other only moderately, which is the same heterogeneity the pooled
#' decomposition attributes to benchmark-conditioned model effects.
benchmark_rank_agreement <- function(df, fit) {
  per_bench <- df |>
    group_by(benchmark, model_name) |>
    summarize(score = mean(score), .groups = "drop") |>
    tidyr::pivot_wider(names_from = benchmark, values_from = score)

  aggregates <- mean_score_aggregate(df) |>
    inner_join(latent_model_effect(fit, df), by = "model_name")

  wide <- per_bench |>
    inner_join(aggregates, by = "model_name") |>
    select(-model_name)

  m <- stats::cor(wide, method = "kendall", use = "pairwise.complete.obs")
  tibble::as_tibble(m, rownames = "benchmark")
}
