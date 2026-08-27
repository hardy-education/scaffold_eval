# scripts/02_benchmark_reliability.R ---------------------------------------
#
# What one benchmark on its own can and cannot tell you.
#
# Produces:
#   Figure 1   inter-scaffold reliability, rho_AA' (Eq. 5)
#   Figure 2   D-study by object of measurement, model vs. model-scaffold
#   Figure 3   signal-to-noise against task count, with the detection band
#   Table      task-only reliability ceilings (Eq. 8) per benchmark
#   Appendix F per-benchmark variance decompositions
#
# Reads the cached per-benchmark fits from scripts/01_fit_models.R.
#
# The three results here make one argument: within a single benchmark, adding
# tasks buys less than it appears to, because the error that limits model
# ranking is indexed by scaffold rather than by task.

source(here::here("R", "setup.R"))

df       <- load_hal_data()
fits     <- fit_benchmarks(df)                      # cached; no estimation
vcs      <- purrr::map(fits, vcomp_from_brms, level = "benchmark")
n_tasks  <- bench_item_counts(df)
n_scaf   <- bench_scaffold_counts(df)

bench_levels <- names(vcs)


# ---- Figure 1: inter-scaffold reliability ----------------------------------
#
# Two scaffolds evaluate the same models on the same tasks. Would they agree
# about the ordering? Averaged over every task the benchmark contains.

inter_scaffold <- purrr::imap(vcs, \(vc, b) {
  ni <- n_tasks$n_tasks[n_tasks$benchmark == b]
  summarize_quantity(inter_scaffold_reliability(vc, n_tasks = ni)) |>
    mutate(benchmark = b, n_tasks = ni, .before = 1)
}) |>
  bind_rows() |>
  mutate(benchmark = factor(benchmark, levels = BENCHMARKS))

print(select(inter_scaffold, benchmark, n_tasks, median, mean, low, high, rhat))
save_table(select(inter_scaffold, benchmark, n_tasks, median, mean, low, high),
           "fig1_inter_scaffold_reliability",
           caption = "Inter-scaffold reliability (Eq. 5) by benchmark.")

save_figure(plot_inter_scaffold(inter_scaffold),
            "fig1_inter_scaffold_reliability", width = 9, height = 4.5)


# ---- Figure 2: D-study for two objects of measurement ----------------------
#
# The same posterior, read two ways. When the model is the object, scaffold
# compatibility is error and creates a ceiling no number of tasks can pass.
# When the deployable system is the object, that same variation is signal.

task_grid <- purrr::map(
  purrr::set_names(bench_levels),
  \(b) unique(round(exp(seq(log(1),
                            log(min(DSTUDY$bench_max_tasks,
                                    n_tasks$n_tasks[n_tasks$benchmark == b])),
                            length.out = 60))))
)

bench_dstudy <- dstudy_by_benchmark(
  vcs, n_tasks = task_grid,
  objects = c("model", "system"), statistics = c("ep2", "snr")
) |>
  mutate(benchmark = factor(benchmark, levels = BENCHMARKS))

readr::write_csv(bench_dstudy |> select(-any_of("draws")),
                 file.path(PATHS$results, "benchmark_dstudy.csv"))

# Eq. 8: the task-only ceiling, evaluated at the scaffold count each benchmark
# actually has. This is the dashed line in Figure 2.
ceilings <- purrr::imap(vcs, \(vc, b) {
  summarize_quantity(reliability_ceiling(vc, n_scaffolds = 1)) |>
    mutate(benchmark = b,
           n_scaffolds_observed = n_scaf$n_scaffolds[n_scaf$benchmark == b],
           .before = 1)
}) |>
  bind_rows() |>
  mutate(benchmark = factor(benchmark, levels = BENCHMARKS))

save_figure(
  plot_bench_reliability(filter(bench_dstudy, statistic == "ep2"), ceilings),
  "fig2_benchmark_reliability_by_object", width = 14.3, height = 3.55
)


# ---- Table: which benchmarks could ever reach the reliability target -------
#
# Section 5.2: only two of nine benchmarks are projected to exceed
# Erho^2 = 0.75 for model ranking, even with unlimited tasks of the same
# construction. Their limitation is not test length.

ceiling_tbl <- ceilings |>
  transmute(
    Benchmark = BENCHMARK_LABELS[as.character(benchmark)],
    `Scaffolds observed` = n_scaffolds_observed,
    `Ceiling (median)` = median,
    `68% HDI low` = low, `68% HDI high` = high,
    `Reaches 0.75?` = ifelse(median > RELIABILITY_TARGET, "yes", "no")
  )

print(ceiling_tbl, n = Inf)
message(sprintf(
  "\n%d of %d benchmarks exceed Erho^2 = %.2f in the task-only limit.",
  sum(ceilings$median > RELIABILITY_TARGET), nrow(ceilings), RELIABILITY_TARGET
))

save_table(ceiling_tbl, "table_task_only_ceilings",
           caption = paste("Task-only reliability ceilings (Eq. 8) for model",
                           "ranking, by benchmark."))


# ---- Figure 3: signal-to-noise ---------------------------------------------

save_figure(
  plot_bench_snr(filter(bench_dstudy, statistic == "snr", object == "model")),
  "fig3_benchmark_snr", width = 14.3, height = 3.55
)

# Where does each benchmark's S/N sit at its observed task count, relative to
# the detection band?
snr_observed <- purrr::imap(vcs, \(vc, b) {
  ni <- n_tasks$n_tasks[n_tasks$benchmark == b]
  summarize_quantity(snr(vc, design(n_tasks = ni))) |>
    mutate(benchmark = b, n_tasks = ni, .before = 1)
}) |>
  bind_rows() |>
  mutate(
    benchmark = factor(benchmark, levels = BENCHMARKS),
    reaches_detection = median >= LOD_BAND[["low"]]
  )

save_table(select(snr_observed, benchmark, n_tasks, median, low, high, reaches_detection),
           "table_snr_at_observed_task_counts",
           caption = "Signal-to-noise ratio at each benchmark's observed task count.")


# ---- Appendix F: per-benchmark variance decompositions ---------------------

bench_shares <- purrr::imap(vcs, \(vc, b) {
  variance_shares(vc) |> mutate(benchmark = b, .before = 1)
}) |>
  bind_rows() |>
  mutate(benchmark = factor(benchmark, levels = BENCHMARKS))

save_table(
  bench_shares |>
    select(benchmark, component, label, median, median_normalized, low, high) |>
    arrange(benchmark, desc(median)),
  "appendix_benchmark_variance_shares",
  caption = "Posterior variance shares within each benchmark (Eq. 2)."
)

p_bench_shares <- bench_shares |>
  mutate(component = factor(component, levels = names(COMPONENT_COLOURS))) |>
  ggplot(aes(x = benchmark, y = median_normalized, fill = component, alpha = component)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = COMPONENT_COLOURS,
                    labels = \(x) component_registry("benchmark")$label[
                      match(x, component_registry("benchmark")$component)]) +
  scale_alpha_manual(values = COMPONENT_ALPHAS, guide = "none") +
  scale_x_discrete(labels = BENCHMARK_LABELS) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of latent variance", fill = "Component") +
  paper_theme(legend = "right") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_figure(p_bench_shares, "appendix_benchmark_variance_shares",
            width = 10, height = 4.5)

message("\n02_benchmark_reliability.R complete.")
