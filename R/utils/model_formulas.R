# Random-effects structures used across analyses (see docs/reliability_models.qmd).

#' Pooled across all benchmarks: full crossed hierarchy.
FORMULA_POOLED_RE <- ~ 1 +
  (1 | benchmark) +
  (1 | benchmark:task_id) +
  (1 | model_name) +
  (1 | agent_name) +
  (1 | benchmark:task_id:model_name) +
  (1 | benchmark:task_id:agent_name) +
  (1 | model_name:agent_name) +
  (1 | benchmark:task_id:model_name:agent_name) +
  (1 | benchmark:model_name) +
  (1 | benchmark:agent_name) +
  (1 | benchmark:model_name:agent_name)

#' Per-benchmark fit: benchmark fixed implicitly by subsetting data.
FORMULA_PER_BENCH_RE <- ~ 1 +
  (1 | task_id) +
  (1 | model_name) +
  (1 | agent_name) +
  (1 | task_id:model_name) +
  (1 | task_id:agent_name) +
  (1 | model_name:agent_name)

#' Per-benchmark with model and agent collapsed into one facet.
FORMULA_MODEL_AGENT_RE <- ~ 1 +
  (1 | task_id) +
  (1 | model)

# Use update(FORMULA_*, score ~ .) to attach the response variable in analysis scripts.
