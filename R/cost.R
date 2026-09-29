# R/cost.R -----------------------------------------------------------------
#
# What an evaluation design costs, and what reliability it buys.
#
# The D-studies in R/gtheory.R project reliability as a function of a design.
# This file attaches a price to that design, which is what turns a measurement
# result into a budgeting one: given a fixed number of dollars, where should
# they go?
#
# Two units matter and are easy to confuse:
#
#   tasks   nb * ni          distinct items in the battery
#   trials  nb * ni * na     rollouts actually executed, one per
#                            (benchmark, task, model, scaffold) cell per model
#
# Reliability is a function of the design; cost is a function of the trials.
# A design that adds scaffolds rather than tasks holds the task count fixed
# while multiplying trials, which is precisely the trade the paper is about.
#
# The price is an *assumption*, not a measurement. The response matrices record
# outcomes, not spend, so `COST$per_trial` in config.R is a flat per-rollout
# rate back-calculated so that the full observed battery reproduces its
# reported total. Real cost varies by benchmark, model, and episode length, and
# none of this counts the engineering time to stand a benchmark up. Treat the
# ordering of designs as trustworthy and the dollar levels as a planning
# approximation; supply `COST$per_trial_by_benchmark` to do better.

#' Number of rollouts a design executes, per model ranked.
#'
#' @param n_tasks,n_benchmarks,n_scaffolds Design counts; vectors recycle.
evaluation_trials <- function(n_tasks, n_benchmarks = 1, n_scaffolds = 1) {
  n_tasks * n_benchmarks * n_scaffolds
}


#' Projected dollar cost of a design.
#'
#' @param n_models Models evaluated. Cost scales with the number of systems
#'   run; reliability does not, so this is a pure multiplier.
#' @param cost_per_trial Mean dollar cost of one rollout.
design_cost <- function(n_tasks, n_benchmarks = 1, n_scaffolds = 1,
                        n_models = COST$n_models,
                        cost_per_trial = COST$per_trial) {
  evaluation_trials(n_tasks, n_benchmarks, n_scaffolds) * n_models * cost_per_trial
}


#' Observed cost of the evaluation actually run.
#'
#' Counts the rollouts in the data rather than assuming a balanced design, so
#' this is the baseline the projected savings are measured against.
observed_cost <- function(df, cost_per_trial = COST$per_trial) {
  tibble(
    n_rollouts = nrow(df),
    cost_per_trial = cost_per_trial,
    total_cost = nrow(df) * cost_per_trial
  )
}


#' Reliability and cost across a grid of designs.
#'
#' The joining of the two is the point: every row is a design, its projected
#' model-ranking reliability, and what it would cost to run. Sorting by cost
#' within a reliability band answers "what is the cheapest way to get here".
#'
#' @param vc A `vcomp` from the pooled fit.
#' @param n_tasks,n_benchmarks,n_scaffolds Design grids to cross.
#' @param ci Credible mass for interval summaries.
cost_reliability_frontier <- function(vc,
                                      n_tasks = c(1, 2, 5, 10, 14, 20, 30, 50, 100, 200),
                                      n_benchmarks = 1:9,
                                      n_scaffolds = 1,
                                      n_models = COST$n_models,
                                      cost_per_trial = COST$per_trial,
                                      ci = CI_LEVEL) {
  dstudy(vc,
         n_tasks = n_tasks, n_benchmarks = n_benchmarks,
         n_scaffolds = n_scaffolds, statistics = "ep2", ci = ci) |>
    mutate(
      trials = evaluation_trials(n_tasks, n_benchmarks, n_scaffolds),
      cost   = design_cost(n_tasks, n_benchmarks, n_scaffolds,
                           n_models, cost_per_trial)
    )
}


#' Cheapest design reaching a target reliability.
#'
#' The allocation question stated directly: of all designs on the frontier that
#' clear `target`, which costs least?
#'
#' @param frontier Output of `cost_reliability_frontier()`.
#' @param target Reliability to reach, on the posterior median.
#' @param baseline_cost Cost to compare against, for the saving column.
cheapest_design_reaching <- function(frontier, target,
                                     baseline_cost = NULL) {
  hit <- frontier |>
    filter(median >= target) |>
    slice_min(cost, n = 1, with_ties = FALSE)

  if (!nrow(hit)) {
    warning("No design in the frontier reaches Erho^2 = ", target,
            "; the maximum projected is ", round(max(frontier$median), 3), ".")
    return(hit)
  }

  hit |>
    mutate(
      target = target,
      baseline_cost = if (is.null(baseline_cost)) NA_real_ else baseline_cost,
      saving = if (is.null(baseline_cost)) NA_real_ else 1 - cost / baseline_cost
    )
}


# ---- Repeated task subsampling ---------------------------------------------

#' Does a smaller task set preserve the published ranking?
#'
#' A design-free check on the same question the D-study answers analytically.
#' Sample k tasks per benchmark without replacement, use the *same* sampled
#' tasks for every model so the comparison is paired, average to the
#' benchmark-model level, then average across benchmarks so each benchmark
#' carries equal weight. Compare the resulting model ranking against the
#' complete-data ranking.
#'
#' Agreement with the full-data ranking establishes information retention. It
#' says nothing about whether that ranking is valid -- a subsample can
#' faithfully reproduce an ordering that was never well identified.
#'
#' @param df Rollout data.
#' @param k_values Tasks per benchmark to try.
#' @param n_rep Subsamples per k.
#' @param seed RNG seed, for reproducibility.
#' @return One row per (k, replicate) with rank correlations against the
#'   full-data ranking.
subsample_rankings <- function(df, k_values = c(5, 10, 15, 20, 25, 30),
                               n_rep = 500, seed = SAMPLER$seed) {
  set.seed(seed)

  # Reference: the ranking implied by all available tasks, benchmarks weighted
  # equally rather than by task count.
  full_rank <- df |>
    group_by(benchmark, model_name) |>
    summarize(score = mean(score), .groups = "drop") |>
    group_by(model_name) |>
    summarize(score = mean(score), .groups = "drop") |>
    mutate(full = rank(-score, ties.method = "average")) |>
    select(model_name, full_score = score, full)

  tasks_by_bench <- df |>
    distinct(benchmark, task_id) |>
    group_by(benchmark) |>
    summarize(tasks = list(task_id), .groups = "drop")

  one_draw <- function(k) {
    keep <- tasks_by_bench |>
      mutate(sampled = purrr::map(tasks, \(t) sample(t, min(k, length(t))))) |>
      select(benchmark, sampled) |>
      tidyr::unnest(sampled) |>
      rename(task_id = sampled)

    sub <- df |>
      inner_join(keep, by = c("benchmark", "task_id")) |>
      group_by(benchmark, model_name) |>
      summarize(score = mean(score), .groups = "drop") |>
      group_by(model_name) |>
      summarize(score = mean(score), .groups = "drop") |>
      inner_join(full_rank, by = "model_name")

    tibble(
      n_models  = nrow(sub),
      spearman  = stats::cor(sub$score, sub$full_score, method = "spearman"),
      kendall   = stats::cor(sub$score, sub$full_score, method = "kendall"),
      # Does the subsample still put the same model on top?
      top1      = sub$model_name[which.max(sub$score)] ==
                  sub$model_name[which.min(sub$full)],
      mean_abs_rank_change =
        mean(abs(rank(-sub$score, ties.method = "average") - sub$full))
    )
  }

  tidyr::expand_grid(k = k_values, replicate = seq_len(n_rep)) |>
    mutate(.res = purrr::map(k, one_draw)) |>
    tidyr::unnest(.res)
}


#' Price of one rollout on each benchmark.
#'
#' Returns the benchmark-specific price table when `COST$per_trial_by_benchmark`
#' is set, and the flat average rate otherwise. Every cost figure in the
#' repository routes through here, so swapping in real prices changes all of
#' them consistently.
benchmark_price <- function(benchmarks,
                            per_trial = COST$per_trial,
                            by_benchmark = COST$per_trial_by_benchmark) {
  if (is.null(by_benchmark)) {
    return(rep(per_trial, length(benchmarks)))
  }
  price <- by_benchmark[as.character(benchmarks)]
  missing <- is.na(price)
  if (any(missing)) {
    warning("No price for benchmark(s): ",
            paste(unique(benchmarks[missing]), collapse = ", "),
            "; falling back to the flat rate.")
    price[missing] <- per_trial
  }
  unname(price)
}


#' Summarise repeated task subsampling across replicates.
#'
#' Reproduces the paper's subsampling table: rank agreement by k, plus the
#' share of the full battery's cost that k tasks per benchmark represents.
#'
#' The cost column weights each retained rollout by its benchmark's price, so
#' it is dollar-weighted when a price table is supplied and
#' rollout-count-weighted under the flat default. Benchmarks differ in both
#' task count and per-rollout cost, so the two weightings do not agree; the
#' flat-rate figures are the ones this repository can reproduce from the
#' shipped data.
summarize_subsampling <- function(subsamples, df,
                                  cost_per_trial = COST$per_trial) {
  spend <- df |>
    mutate(.price = benchmark_price(benchmark, cost_per_trial)) |>
    group_by(benchmark) |>
    summarize(n_rollouts = n(), price = dplyr::first(.price), .groups = "drop") |>
    inner_join(bench_item_counts(df), by = "benchmark")

  total_cost <- sum(spend$n_rollouts * spend$price)

  cost_share <- function(k) {
    # Cost retained when each benchmark is capped at k tasks, holding the
    # observed model and scaffold coverage fixed.
    kept <- sum(spend$n_rollouts * pmin(k, spend$n_tasks) / spend$n_tasks *
                  spend$price)
    1 - kept / total_cost
  }

  subsamples |>
    group_by(k) |>
    summarize(
      mean_spearman   = mean(spearman),
      median_spearman = stats::median(spearman),
      min_spearman    = min(spearman),
      max_spearman    = max(spearman),
      mean_kendall    = mean(kendall),
      median_kendall  = stats::median(kendall),
      min_kendall     = min(kendall),
      max_kendall     = max(kendall),
      top1_agreement  = mean(top1),
      mean_abs_rank_change = mean(mean_abs_rank_change),
      .groups = "drop"
    ) |>
    mutate(cost_saving = purrr::map_dbl(k, cost_share))
}
