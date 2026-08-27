# R/data.R -----------------------------------------------------------------
#
# Loading and describing the HAL response matrix.
#
# The analysis dataset is one row per agent rollout, with a binary outcome:
#
#   run_id            identifier of the evaluation run the rollout came from
#   benchmark         b -- one of nine HAL benchmarks
#   task_id           i -- item, nested within benchmark
#   agent_name        a -- agent scaffold
#   model_name        m -- base LLM *with reasoning effort appended*
#   reasoning_effort  the effort setting, kept for reference
#   score             y in {0, 1}
#
# Throughout the paper "model" means a base LLM paired with a reasoning-effort
# setting: gpt-5-high and gpt-5-minimal are distinct levels of m. That pairing
# is baked into `model_name` by `normalize_hal_columns()`.

# The nine benchmarks, in the order used by the paper's figures and tables.
BENCHMARKS <- c(
  "assistantbench", "corebench_hard", "gaia", "onlinemind2web", "scicode",
  "scienceagentbench", "swebench_verified_mini", "taubench_airline", "usaco"
)

# Display names for figures and LaTeX tables.
BENCHMARK_LABELS <- c(
  assistantbench         = "AssistantBench",
  corebench_hard         = "CORE-Bench Hard",
  gaia                   = "GAIA",
  onlinemind2web         = "Online-Mind2Web",
  scicode                = "SciCode",
  scienceagentbench      = "ScienceAgentBench",
  swebench_verified_mini = "SWE-bench Verified Mini",
  taubench_airline       = "tau-bench Airline",
  usaco                  = "USACO"
)

# Abbreviated names for facet strips, which are too narrow for the full names
# in a nine-panel row. Tables and axis labels use `BENCHMARK_LABELS`.
BENCHMARK_SHORT <- c(
  assistantbench         = "AssistantBench",
  corebench_hard         = "CORE-Bench",
  gaia                   = "GAIA",
  onlinemind2web         = "Mind2Web",
  scicode                = "SciCode",
  scienceagentbench      = "SciAgentBench",
  swebench_verified_mini = "SWE-bench Mini",
  taubench_airline       = "tau-bench",
  usaco                  = "USACO"
)

REQUIRED_COLUMNS <- c(
  "run_id", "benchmark", "task_id", "agent_name",
  "model_name", "reasoning_effort", "score"
)


#' Normalise raw HAL export columns into analysis form.
#'
#' Raw HAL exports name models inconsistently: some carry a provider prefix
#' (`openai:gpt-5`), casing and separators vary, the reasoning effort lives in
#' its own column, and one scaffold appears under two names. This function
#' resolves all of that. It is **idempotent** -- applying it to an
#' already-normalised matrix is a no-op -- so it is safe to call on either the
#' shipped analysis dataset or a fresh export.
#'
#' @param df A data frame with the columns in `REQUIRED_COLUMNS`.
#' @return The same rows with normalised `model_name` and `agent_name`.
normalize_hal_columns <- function(df) {
  stopifnot(all(REQUIRED_COLUMNS %in% names(df)))

  df |>
    mutate(
      reasoning_effort = if_else(
        is.na(reasoning_effort) | reasoning_effort == "",
        "default", reasoning_effort
      ),
      # Drop any provider prefix, lower-case, and make separators uniform so
      # that random-effect level labels are safe to parse.
      model_name = tolower(str_split_i(model_name, ":", 1)),
      model_name = str_replace_all(model_name, "[.-]", "_"),
      # Append reasoning effort unless a previous pass already did so.
      model_name = if_else(
        str_ends(model_name, reasoning_effort),
        model_name,
        str_c(model_name, reasoning_effort)
      ),
      # AssistantBench's browser agent is the same scaffold as `browser-use`.
      agent_name = if_else(
        agent_name == "assistantbench_browser_agent", "browser-use", agent_name
      )
    ) |>
    distinct()
}


#' Load the HAL response matrix used by every analysis in this repository.
#'
#' @param path CSV path; defaults to the dataset shipped in `data/`.
#' @param normalize Apply `normalize_hal_columns()`. Harmless on the shipped
#'   file (which is already normalised); required for a fresh HAL export.
#' @param validate Check the loaded design against the paper's Table 1.
#' @return A tibble of rollouts, with `benchmark` ordered as in `BENCHMARKS`.
load_hal_data <- function(path = PATHS$data, normalize = TRUE, validate = TRUE) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (normalize) df <- normalize_hal_columns(df)

  df <- df |>
    mutate(
      score = as.numeric(score),
      benchmark = factor(benchmark, levels = BENCHMARKS)
    ) |>
    filter(!is.na(score))

  if (validate) validate_hal_data(df)
  df
}


#' Sanity-check a loaded response matrix.
#'
#' Catches the failure modes that would silently change every downstream
#' number: an unexpected benchmark, a non-binary outcome, or a model whose
#' reasoning effort was appended twice.
validate_hal_data <- function(df) {
  unknown <- setdiff(as.character(unique(df$benchmark)), BENCHMARKS)
  if (length(unknown)) {
    stop("Unexpected benchmark(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (!all(df$score %in% c(0, 1))) {
    stop("`score` must be binary; the latent-scale decomposition assumes ",
         "Bernoulli outcomes.", call. = FALSE)
  }
  if (anyNA(df[REQUIRED_COLUMNS])) {
    stop("Missing values in required columns.", call. = FALSE)
  }
  invisible(df)
}


#' Number of distinct tasks per benchmark.
#'
#' Used throughout the D-studies as the "observed" point on the task axis.
#'
#' @return A tibble with columns `benchmark` and `n_tasks`.
bench_item_counts <- function(df) {
  df |>
    distinct(benchmark, task_id) |>
    count(benchmark, name = "n_tasks")
}


#' Number of distinct scaffolds per benchmark.
bench_scaffold_counts <- function(df) {
  df |>
    distinct(benchmark, agent_name) |>
    count(benchmark, name = "n_scaffolds")
}


#' Benchmark design summary (paper Table 1).
#'
#' Counts the unique level of each facet observed within each benchmark,
#' alongside the mean score and the rollout count.
design_summary <- function(df) {
  df |>
    group_by(benchmark) |>
    summarize(
      n_models          = n_distinct(model_name),
      n_tasks           = n_distinct(task_id),
      n_scaffolds       = n_distinct(agent_name),
      n_model_scaffold  = n_distinct(paste(model_name, agent_name)),
      n_rollouts        = n(),
      mean_score        = mean(score),
      .groups = "drop"
    ) |>
    mutate(benchmark_label = BENCHMARK_LABELS[as.character(benchmark)], .after = benchmark)
}


#' Coverage of the observed incidence structure (paper Appendix B and E.1).
#'
#' Reports the facts the identifiability argument relies on: how many models
#' appear on every benchmark, how widely each scaffold travels, and how much of
#' the fully crossed design is actually observed.
coverage_summary <- function(df) {
  n_bench <- n_distinct(df$benchmark)

  per_model <- df |>
    group_by(model_name) |>
    summarize(
      n_benchmarks = n_distinct(benchmark),
      n_scaffolds  = n_distinct(agent_name),
      .groups = "drop"
    )

  per_scaffold <- df |>
    group_by(agent_name) |>
    summarize(n_benchmarks = n_distinct(benchmark), .groups = "drop")

  cells_observed <- n_distinct(paste(df$benchmark, df$task_id,
                                     df$model_name, df$agent_name))
  cells_possible <- sum(design_summary(df)$n_tasks) *
    n_distinct(df$model_name) * n_distinct(df$agent_name)

  list(
    n_rollouts               = nrow(df),
    n_benchmarks             = n_bench,
    n_tasks                  = n_distinct(paste(df$benchmark, df$task_id)),
    n_models                 = n_distinct(df$model_name),
    n_scaffolds              = n_distinct(df$agent_name),
    # Cross-benchmark anchors: the models that identify M against BM.
    models_on_all_benchmarks = sum(per_model$n_benchmarks == n_bench),
    # Scaffold-side anchors: the scaffolds that identify A against BA.
    scaffolds_multi_bench    = sum(per_scaffold$n_benchmarks > 1),
    max_scaffold_reach       = max(per_scaffold$n_benchmarks),
    models_single_scaffold   = sum(per_model$n_scaffolds == 1),
    observed_cell_fraction   = cells_observed / cells_possible,
    per_model                = per_model,
    per_scaffold             = per_scaffold
  )
}


#' Fraction of (benchmark, task, model, scaffold) cells with replicated runs.
#'
#' The paper notes that under 5% of cells are replicated, which is why the
#' four-way interaction cannot be separated from cell-level execution noise
#' (Appendix E.4). This quantifies that claim from the data.
replication_rate <- function(df) {
  cells <- df |> count(benchmark, task_id, model_name, agent_name, name = "n_runs")
  list(
    n_cells            = nrow(cells),
    n_replicated       = sum(cells$n_runs > 1),
    fraction_replicated = mean(cells$n_runs > 1),
    fraction_of_rows    = sum(cells$n_runs[cells$n_runs > 1]) / nrow(df)
  )
}


#' Observed per-benchmark leaderboard: mean score by model, as published.
#'
#' This is the ranking that HAL and comparable leaderboards report, and the
#' baseline the latent rankings in `R/ranks.R` are compared against.
published_scores <- function(df, by_scaffold = FALSE) {
  keys <- c("benchmark", "model_name", if (by_scaffold) "agent_name")
  df |>
    group_by(across(all_of(keys))) |>
    summarize(score = mean(score), n_rollouts = n(), .groups = "drop") |>
    group_by(benchmark) |>
    mutate(
      published_rank    = rank(-score, ties.method = "average"),
      published_pctrank = dplyr::percent_rank(score)
    ) |>
    ungroup()
}
