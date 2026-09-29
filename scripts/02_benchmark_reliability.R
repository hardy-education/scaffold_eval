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


# ---- Every coefficient, per benchmark --------------------------------------
#
# The same posterior read seven ways. Each column answers a different question,
# and the spread across them is the paper's first finding: a benchmark can be
# excellent for one claim and near-useless for another.
#
#   Erho^2_M(i), Erho^2_A(i)  how task-dependent is each facet (descriptive;
#                             see `task_reliability()` on what it omits)
#   rho_AA'                   would two scaffolds rank the models alike
#   Erho^2_M(b)               ranking models, scaffold treated as error
#   Erho^2_A(b)               ranking scaffolds, model treated as error
#   Erho^2_MA(b)              ranking deployable model-scaffold systems
#   Erho^2_M(b)(inf)          the task-only ceiling for model ranking

coefficient_table <- purrr::imap(vcs, function(vc, b) {
  ni <- n_tasks$n_tasks[n_tasks$benchmark == b]
  med <- function(x) summarize_quantity(x)$median

  tibble(
    benchmark      = b,
    n_tasks        = ni,
    ep2_M_task     = med(task_reliability(vc, "model")),
    ep2_A_task     = med(task_reliability(vc, "scaffold")),
    rho_AA         = med(inter_scaffold_reliability(vc, n_tasks = ni)),
    ep2_M_bench    = med(ep2(vc, design(n_tasks = ni), object = "model")),
    ep2_A_bench    = med(ep2(vc, design(n_tasks = ni), object = "scaffold")),
    ep2_MA_bench   = med(ep2(vc, design(n_tasks = ni), object = "system")),
    ep2_M_ceiling  = med(reliability_ceiling(vc, object = "model"))
  )
}) |>
  bind_rows() |>
  mutate(benchmark = factor(benchmark, levels = BENCHMARKS)) |>
  arrange(benchmark)

print(coefficient_table, n = Inf)
save_table(coefficient_table, "table_benchmark_coefficients",
           caption = paste("Reliability and generalizability coefficients per",
                           "benchmark, by object of measurement."))

message(sprintf(
  "\nRanking deployable systems is reliable everywhere (Erho^2_MA in [%.3f, %.3f]);\nranking the underlying model is not (Erho^2_M in [%.3f, %.3f]).",
  min(coefficient_table$ep2_MA_bench), max(coefficient_table$ep2_MA_bench),
  min(coefficient_table$ep2_M_bench), max(coefficient_table$ep2_M_bench)
))
message(sprintf(
  "Inter-scaffold reliability ranges from %.3f (%s) to %.3f (%s).",
  min(coefficient_table$rho_AA),
  coefficient_table$benchmark[which.min(coefficient_table$rho_AA)],
  max(coefficient_table$rho_AA),
  coefficient_table$benchmark[which.max(coefficient_table$rho_AA)]
))


# ---- Model against scaffold, within each benchmark -------------------------
#
# The pooled fit answers this with more power, but the per-benchmark fits show
# whether the pattern is general or driven by one or two instruments. Two
# contrasts, because they point different ways: models and scaffolds differ
# comparably on average, while the scaffold has the larger effect on *which
# tasks* get solved.

bench_contrasts <- purrr::imap(vcs, function(vc, b) {
  model_vs_scaffold_contrast(vc, design(1, 1, 1), scope = "both") |>
    mutate(benchmark = b, .before = 1)
}) |>
  bind_rows() |>
  mutate(benchmark = factor(benchmark, levels = BENCHMARKS))

print(bench_contrasts |>
        select(benchmark, label, median, low, high, p_scaffold_exceeds_model),
      n = Inf)

save_table(bench_contrasts |>
             select(benchmark, scope, label, median, low, high,
                    p_model_exceeds_scaffold, p_scaffold_exceeds_model),
           "table_benchmark_model_vs_scaffold",
           caption = paste("Per-benchmark contrasts between model-related and",
                           "scaffold-related variance: main effects and",
                           "task-indexed interactions."))

save_figure(plot_benchmark_contrasts(bench_contrasts),
            "fig2_benchmark_model_vs_scaffold", width = 10, height = 4.5)

for (s_ in unique(bench_contrasts$scope)) {
  d <- filter(bench_contrasts, scope == s_)
  message(sprintf(
    "%-18s scaffold exceeds model on %d of %d benchmarks (P > 0.5).",
    dplyr::first(d$label), sum(d$p_scaffold_exceeds_model > 0.5), nrow(d)
  ))
}


# ---- Absolute reliability, on the observed scale ---------------------------
#
# Everything above concerns *rankings*. A threshold claim -- "this system
# clears 60%" -- is a claim about a level, and levels are reproducible only if
# the components that shift every model together are also small. That is the
# dependability coefficient Phi, and it is always at most Erho^2.
#
# Phi is quoted on the observed proportion scale rather than the latent
# log-odds scale, because that is the scale a threshold is stated on. That
# needs the Gaussian per-benchmark fits, which are opt-in:
#
#   RUN_ABSOLUTE=TRUE REFIT=TRUE Rscript scripts/02_benchmark_reliability.R

RUN_ABSOLUTE <- env_flag("RUN_ABSOLUTE")

if (RUN_ABSOLUTE) {
  message("\n=== Absolute reliability (Gaussian fits, observed scale) ===")

  gaussian_vcs <- purrr::map(purrr::set_names(bench_levels), function(b) {
    f <- fit_bayes(filter(df, benchmark == b), FORMULA_BENCHMARK,
                   paste0("benchmark_gaussian_", b),
                   family = gaussian(), sampler = SAMPLER$benchmark)
    vcomp_from_brms(f, level = "benchmark")
  })

  absolute <- purrr::imap(gaussian_vcs, function(vc, b) {
    ni <- n_tasks$n_tasks[n_tasks$benchmark == b]
    na <- n_scaf$n_scaffolds[n_scaf$benchmark == b]
    bind_cols(
      tibble(benchmark = b, n_tasks = ni, n_scaffolds = na),
      tibble(
        phi_M     = summarize_quantity(phi(vc, design(ni, 1, na)))$median,
        ep2_M     = summarize_quantity(ep2(vc, design(ni, 1, na)))$median,
        phi_M_one_task = summarize_quantity(phi(vc, design(1, 1, na)))$median
      )
    )
  }) |>
    bind_rows() |>
    mutate(benchmark = factor(benchmark, levels = BENCHMARKS)) |>
    arrange(benchmark)

  print(absolute, n = Inf)
  save_table(absolute, "table_absolute_reliability",
             caption = paste("Absolute (dependability) and relative reliability",
                             "for model scores, estimated on the observed scale."))

  message(sprintf(
    "Absolute reliability is at or below relative reliability on every benchmark (max Phi = %.3f, max Erho^2 = %.3f).",
    max(absolute$phi_M), max(absolute$ep2_M)
  ))
  message("A stable ranking does not imply reproducible score levels: ",
          "these data support ordering claims far better than threshold claims.")
} else {
  message("\nSkipping absolute reliability ",
          "(set RUN_ABSOLUTE=TRUE REFIT=TRUE; needs Gaussian per-benchmark fits).")
}


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
