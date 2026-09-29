# tests/test_design.R ------------------------------------------------------
#
# Checks the data-facing machinery added for the corroboration, external
# validation, and cost analyses. None of it needs a fitted model, so this runs
# in seconds and is the fastest way to catch a broken data file.
#
# Run with:  Rscript tests/test_design.R

source(here::here("R", "setup.R"))

failures <- 0L
check <- function(label, actual, expected, tol = 1e-8) {
  ok <- if (is.character(expected) || is.character(actual)) {
    identical(as.character(actual), as.character(expected))
  } else {
    isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = tol))
  }
  cat(if (ok) "  ok   " else "  FAIL ", label, "\n", sep = "")
  if (!ok) {
    cat("         actual:   ", format(actual), "\n",
        "         expected: ", format(expected), "\n", sep = "")
    failures <<- failures + 1L
  }
  invisible(ok)
}

hal <- load_hal_data()


cat("\nHarbor Index loader\n")

harbor_all <- load_harbor_data(min_tasks = NULL)
harbor     <- load_harbor_data()

check("the adjudicated outcome is binary",
      all(harbor_all$score %in% c(0, 1)), TRUE)
check("the inclusion threshold drops benchmarks with too few tasks",
      all(harbor |> count(benchmark, task_id) |> count(benchmark) |>
            pull(n) >= HARBOR_MIN_TASKS), TRUE)
check("filtering removes rows but no models or scaffolds",
      c(n_distinct(harbor$model_name), n_distinct(harbor$agent_name)),
      c(n_distinct(harbor_all$model_name), n_distinct(harbor_all$agent_name)))
check("Harbor uses the same facet column names as HAL",
      all(c("benchmark", "task_id", "model_name", "agent_name", "score")
          %in% names(harbor)), TRUE)

# Harbor's value as a check is that every benchmark sees the same models and
# scaffolds, which is where HAL is uneven. If that stops being true the
# corroboration argument weakens, so it is asserted rather than assumed.
per_bench <- harbor |>
  group_by(benchmark) |>
  summarize(n_models = n_distinct(model_name),
            n_scaffolds = n_distinct(agent_name), .groups = "drop")
check("every Harbor benchmark sees the same number of models",
      n_distinct(per_bench$n_models), 1)
check("every Harbor benchmark sees the same number of scaffolds",
      n_distinct(per_bench$n_scaffolds), 1)

# Harbor is *not* fully crossed on model by scaffold, which bounds what it can
# corroborate. Pinned so the claim in the docs cannot drift from the data.
crossing <- harbor |>
  summarize(pairs = n_distinct(paste(model_name, agent_name)),
            possible = n_distinct(model_name) * n_distinct(agent_name))
check("Harbor model-scaffold coverage is partial, as documented",
      crossing$pairs < crossing$possible, TRUE)

# A decomposition needs the benchmark facet to vary.
check("more than one benchmark survives the threshold",
      n_distinct(harbor$benchmark) > 1, TRUE)


cat("\nExternal validation data\n")

edf <- load_external_data()
ext <- external_scores(edf)

check("every external benchmark has a declared overlap rule",
      sort(unique(as.character(edf$benchmark))),
      sort(names(EXTERNAL_EXCLUSIONS)))
check("external scores are averaged to one row per benchmark-model",
      nrow(ext), nrow(distinct(edf, benchmark, model_name)))
check("every external benchmark has a display label",
      all(!is.na(EXTERNAL_LABELS[as.character(unique(edf$benchmark))])), TRUE)

# The contamination rules must name benchmarks that actually exist, or the
# exclusion silently does nothing.
overlaps <- stats::na.omit(EXTERNAL_EXCLUSIONS)
check("benchmarks named in the exclusion rules exist in the HAL data",
      all(overlaps %in% levels(hal$benchmark)), TRUE)

# Each comparison needs enough overlapping models to correlate at all.
overlap_counts <- ext |>
  filter(model_name %in% unique(hal$model_name)) |>
  count(benchmark)
check("every external benchmark overlaps HAL by at least six models",
      all(overlap_counts$n >= 6), TRUE)

# Excluding the overlapping benchmark must actually change the baseline.
base_full <- mean_score_aggregate(hal)
base_less <- mean_score_aggregate(hal, exclude = "swebench_verified_mini")
check("excluding a benchmark changes the mean-score aggregate",
      isTRUE(all.equal(base_full$mean_score, base_less$mean_score)), FALSE)
# Dropping a benchmark can drop models evaluated only there. That is correct
# behaviour, and it is why the external comparison reports its own model count
# per benchmark rather than assuming a fixed panel.
check("excluding a benchmark cannot add models",
      nrow(base_less) <= nrow(base_full), TRUE)


cat("\nCost model\n")

check("trials multiply across facets", evaluation_trials(7, 3, 2), 42)
check("cost is linear in the per-trial price",
      design_cost(10, 2, 1, n_models = 3, cost_per_trial = 4) /
        design_cost(10, 2, 1, n_models = 3, cost_per_trial = 2), 2)
check("the flat price reproduces the reported battery total",
      observed_cost(hal)$total_cost, COST$observed_total, tol = 1e-6)

check("a benchmark price table overrides the flat rate",
      benchmark_price(c("gaia", "usaco"), by_benchmark = c(gaia = 5, usaco = 2)),
      c(5, 2))
check("a missing price falls back to the flat rate",
      suppressWarnings(benchmark_price("scicode", by_benchmark = c(gaia = 5))),
      COST$per_trial)


cat("\nTask subsampling\n")

set.seed(1)
sub <- subsample_rankings(hal, k_values = c(5, 30), n_rep = 12)

check("one row per (k, replicate)", nrow(sub), 2 * 12)
check("rank correlations are in range",
      all(sub$spearman >= -1 & sub$spearman <= 1), TRUE)
check("more tasks agree better with the full-data ranking",
      mean(sub$spearman[sub$k == 30]) > mean(sub$spearman[sub$k == 5]), TRUE)
check("rank displacement shrinks as tasks are added",
      mean(sub$mean_abs_rank_change[sub$k == 30]) <
        mean(sub$mean_abs_rank_change[sub$k == 5]), TRUE)

# Sampling every task must reproduce the full-data ranking exactly.
max_tasks <- max(bench_item_counts(hal)$n_tasks)
full <- subsample_rankings(hal, k_values = max_tasks, n_rep = 2)
check("sampling every task recovers the complete-data ranking",
      all(abs(full$spearman - 1) < 1e-9), TRUE)
check("and leaves no model displaced",
      all(full$mean_abs_rank_change == 0), TRUE)

summary_tbl <- summarize_subsampling(sub, hal)
check("cost saving falls as tasks are retained",
      summary_tbl$cost_saving[summary_tbl$k == 30] <
        summary_tbl$cost_saving[summary_tbl$k == 5], TRUE)
check("cost savings are proportions",
      all(summary_tbl$cost_saving >= 0 & summary_tbl$cost_saving <= 1), TRUE)


cat("\n", if (failures == 0L) "All checks passed.\n" else
    sprintf("%d check(s) FAILED.\n", failures), sep = "")
quit(status = if (failures == 0L) 0L else 1L)
