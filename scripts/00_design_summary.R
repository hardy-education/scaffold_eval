# scripts/00_design_summary.R ----------------------------------------------
#
# Describes the evaluation design before any model is fitted.
#
# Produces:
#   Table 1        benchmark design summary
#   Appendix B     coverage of the observed incidence structure
#   Appendix E.1   the overlaps the identifiability argument relies on
#
# Runs in seconds and requires no fitted models: this is the script to run
# first when checking that the dataset loaded as expected.

source(here::here("R", "setup.R"))

df <- load_hal_data()

# ---- Table 1: benchmark design summary -------------------------------------

summary_tbl <- design_summary(df) |>
  transmute(
    Benchmark   = benchmark_label,
    `Models`    = n_models,
    `Tasks`     = n_tasks,
    `Scaffolds` = n_scaffolds,
    `Model-Scaffold pairs` = n_model_scaffold,
    `Rollouts`  = n_rollouts,
    `Avg. score` = round(mean_score, 2)
  )

print(summary_tbl, n = Inf)
save_table(summary_tbl, "table1_design_summary", digits = 2,
           caption = "Benchmark design summary. Counts denote unique levels of each facet.")


# ---- Coverage and connectivity (Appendix B, E.1) ---------------------------

cov <- coverage_summary(df)
rep <- replication_rate(df)

cat("\nObserved design\n")
cat(sprintf("  rollouts                              %s\n", format(cov$n_rollouts, big.mark = ",")))
cat(sprintf("  benchmarks / tasks / models / scaffolds  %d / %s / %d / %d\n",
            cov$n_benchmarks, format(cov$n_tasks, big.mark = ","),
            cov$n_models, cov$n_scaffolds))
cat(sprintf("  fraction of the fully crossed design observed   %.1f%%\n",
            100 * cov$observed_cell_fraction))

cat("\nConnectivity anchors (what makes M separable from BM, and A from BA)\n")
cat(sprintf("  models evaluated on every benchmark   %d\n", cov$models_on_all_benchmarks))
cat(sprintf("  widest scaffold reach                 %d benchmarks\n", cov$max_scaffold_reach))
cat(sprintf("  scaffolds used on >1 benchmark        %d\n", cov$scaffolds_multi_bench))
cat(sprintf("  models seen under a single scaffold   %d\n", cov$models_single_scaffold))

cat("\nCell replication (why the four-way interaction is not separable)\n")
cat(sprintf("  (benchmark, task, model, scaffold) cells   %s\n", format(rep$n_cells, big.mark = ",")))
cat(sprintf("  cells with more than one rollout          %.1f%%\n", 100 * rep$fraction_replicated))

coverage_tbl <- bind_rows(
  tibble(quantity = "Rollouts",                        value = cov$n_rollouts),
  tibble(quantity = "Benchmarks",                      value = cov$n_benchmarks),
  tibble(quantity = "Tasks (nested in benchmark)",     value = cov$n_tasks),
  tibble(quantity = "Models (LLM x reasoning effort)", value = cov$n_models),
  tibble(quantity = "Scaffolds",                       value = cov$n_scaffolds),
  tibble(quantity = "Models on all benchmarks",        value = cov$models_on_all_benchmarks),
  tibble(quantity = "Widest scaffold reach",           value = cov$max_scaffold_reach),
  tibble(quantity = "Replicated cells (%)",            value = round(100 * rep$fraction_replicated, 1))
)
save_table(coverage_tbl, "coverage_summary", digits = 1)


# ---- Per-benchmark scaffold roster (Appendix B) ----------------------------

scaffold_tbl <- df |>
  distinct(benchmark, agent_name) |>
  arrange(benchmark, agent_name) |>
  group_by(benchmark) |>
  summarize(scaffolds = paste(agent_name, collapse = ", "), .groups = "drop") |>
  mutate(Benchmark = BENCHMARK_LABELS[as.character(benchmark)], .before = 1) |>
  select(-benchmark)

save_table(scaffold_tbl, "scaffold_roster", digits = 0)


# ---- Coverage figures ------------------------------------------------------

# How broadly is each model evaluated? Half the scaffolds appear on a single
# benchmark, and many models are only ever seen under one scaffold; this is the
# sparsity the partial-pooling model has to work around.
p_coverage <- cov$per_model |>
  tidyr::pivot_longer(c(n_benchmarks, n_scaffolds),
                      names_to = "facet", values_to = "count") |>
  mutate(facet = recode(facet,
                        n_benchmarks = "Benchmarks evaluated on",
                        n_scaffolds  = "Scaffolds evaluated under")) |>
  ggplot(aes(count)) +
  geom_bar(fill = COMPONENT_COLOURS[["M"]]) +
  facet_wrap(~facet, scales = "free_x") +
  labs(x = "Count per model", y = "Number of models") +
  paper_theme(legend = "none")

save_figure(p_coverage, "coverage_per_model", width = 8, height = 3.2)

# Score distributions: benchmarks differ enormously in difficulty, and several
# sit near the floor, which is why the decomposition is estimated on the latent
# logit scale rather than treating proportions as Gaussian.
p_scores <- df |>
  group_by(benchmark, agent_name, model_name) |>
  summarize(score = mean(score), .groups = "drop") |>
  ggplot(aes(score, fill = agent_name)) +
  geom_histogram(bins = 20) +
  facet_wrap(~benchmark, scales = "free_y", labeller = label_benchmarks()) +
  labs(x = "Mean score per (model, scaffold)", y = "Count", fill = "Scaffold") +
  paper_theme(legend = "bottom")

save_figure(p_scores, "score_distributions", width = 11, height = 6)

# Record the environment next to the results, so a later discrepancy in a
# variance component can be traced to a package or sampler version.
record_session()

message("\n00_design_summary.R complete.")
