# R/harbor.R ---------------------------------------------------------------
#
# The Harbor Index corroboration dataset.
#
# A second, independently constructed meta-benchmark used to check whether the
# design conclusions drawn from HAL survive on different data. Harbor Index
# distils 82 tasks from a pool of roughly 6,000 candidates across 54
# benchmarks, filtered for difficulty and audited for correctness, and runs
# them through 9 models and 4 scaffolds.
#
# How it differs from HAL, precisely. Harbor runs the *same* 9 models and 4
# scaffolds on every benchmark, so benchmark-by-model and benchmark-by-scaffold
# coverage is complete; in HAL both vary from benchmark to benchmark and the
# pooled fit leans on a handful of shared anchors to connect them. That is the
# sparsity Harbor is a check against.
#
# Harbor is *not* fully crossed on model by scaffold: only 18 of 36 pairs are
# observed, and the coverage is lopsided -- terminus-2 runs all 9 models,
# claude-code 7, and codex and gemini-cli one each. Scaffold-indexed components
# are correspondingly less well identified here than model-indexed ones, and
# Harbor's value as corroboration lies in its benchmark coverage rather than in
# resolving scaffold effects.
#
# The decomposition is the same one the leaderboard fit uses
# (`FORMULA_LEADERBOARD`), so nothing in R/gtheory.R needs to know which
# dataset a `vcomp` came from.

# Harbor names its facets differently; these map onto the columns the shared
# formulas and modules expect.
HARBOR_COLUMNS <- c(
  model_name = "model",
  agent_name = "agent",
  task_id    = "task_id",
  benchmark  = "benchmark"
)

#' Minimum tasks a benchmark needs to enter the decomposition.
#'
#' Task and benchmark effects cannot be separated for a benchmark represented
#' by one or two items. Three is the smallest threshold that leaves both
#' estimable; the paper reports that raising it to four or five drops the
#' estimated model share further, which it attributes to the non-random way
#' tasks were selected to represent each benchmark rather than to a real
#' change in signal. `scripts/07_harbor_corroboration.R` re-runs the threshold
#' sweep so that sensitivity is visible rather than asserted.
HARBOR_MIN_TASKS <- 3


#' Load the Harbor Index trial-level data.
#'
#' @param path CSV path; defaults to the copy shipped in `data/`.
#' @param min_tasks Drop benchmarks with fewer than this many distinct tasks.
#'   `NULL` keeps everything, for the descriptive table.
#' @return A tibble using the same column names as the HAL data, so the shared
#'   formulas and fitting code apply unchanged.
load_harbor_data <- function(path = PATHS$harbor,
                             min_tasks = HARBOR_MIN_TASKS) {
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)

  required <- c("benchmark", "task_id", "model", "agent", "judged_score")
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop("Harbor data is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  out <- raw |>
    rename(model_name = model, agent_name = agent) |>
    # `judged_score` is the adjudicated outcome; see data description in paper
    #  for how it relates to the verifier reward and the judge's classification.
    mutate(score = as.numeric(judged_score)) |>
    filter(!is.na(score))

  if (!is.null(min_tasks)) {
    keep <- out |>
      group_by(benchmark) |>
      summarize(n_tasks = n_distinct(task_id), .groups = "drop") |>
      filter(n_tasks >= min_tasks)
    out <- filter(out, benchmark %in% keep$benchmark)
  }

  out |>
    mutate(benchmark = factor(benchmark, levels = sort(unique(benchmark)))) |>
    validate_harbor_data()
}


#' Sanity-check the Harbor data the same way the HAL loader does.
validate_harbor_data <- function(df) {
  if (!all(df$score %in% c(0, 1))) {
    stop("Harbor `score` must be binary after adjudication; found ",
         paste(utils::head(setdiff(unique(df$score), c(0, 1)), 3), collapse = ", "),
         call. = FALSE)
  }
  if (anyNA(df[c("benchmark", "task_id", "model_name", "agent_name")])) {
    stop("Missing values in Harbor facet columns.", call. = FALSE)
  }
  invisible(df)
}


#' Harbor design summary (paper Appendix table of descriptive statistics).
#'
#' Call with `min_tasks = NULL` data to reproduce the full table, which marks
#' which benchmarks clear the inclusion threshold.
harbor_design_summary <- function(df, min_tasks = HARBOR_MIN_TASKS) {
  df |>
    group_by(benchmark) |>
    summarize(
      n_models         = n_distinct(model_name),
      n_tasks          = n_distinct(task_id),
      n_scaffolds      = n_distinct(agent_name),
      n_model_scaffold = n_distinct(paste(model_name, agent_name)),
      n_rollouts       = n(),
      mean_score       = mean(score),
      .groups = "drop"
    ) |>
    mutate(included = n_tasks >= min_tasks) |>
    arrange(benchmark)
}


#' Domain each Harbor benchmark belongs to.
#'
#' Harbor groups its benchmarks into topical domains. These are carried in the
#' data rather than recomputed here; this returns the mapping for reporting.
harbor_domains <- function(df) {
  df |>
    distinct(benchmark, domain) |>
    arrange(domain, benchmark)
}


#' Fit the pooled decomposition to the Harbor Index data.
#'
#' Uses `FORMULA_LEADERBOARD` unchanged: Harbor has the same facet structure
#' as the HAL leaderboard, so the same variance components are estimated and
#' the same `vcomp` registry applies.
fit_harbor <- function(df, name = "harbor_bernoulli", ...) {
  fit_bayes(df, FORMULA_LEADERBOARD, name, sampler = SAMPLER$harbor, ...)
}
