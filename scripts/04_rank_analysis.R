# scripts/04_rank_analysis.R -----------------------------------------------
#
# What the reliability estimates mean for published leaderboards.
#
# Produces:
#   Figure 6   published rank against posterior rank, with credible intervals
#   Figure 7   slope chart of the reordering, for the featured benchmarks
#   Table 2    agreement between observed and latent rankings
#   Table D.5  full-model versus per-benchmark latent rankings
#   Section 5.6  share of model pairs the posterior cannot order
#
# The point of these outputs is not that the posterior ordering is a corrected
# "true" ranking. It is that published ranks are sensitive to a defensible
# choice -- which variation gets credited to the model -- and that many of the
# rank differences a leaderboard displays are not resolved by the data.

source(here::here("R", "setup.R"))

df   <- load_hal_data()
fit  <- fit_leaderboard(df)
pub  <- published_scores(df)

# theta_mb = u_m + u_bm : the model's persistent capability plus its
# benchmark-specific departure, with scaffold effects marginalized out.
theta <- latent_capability_draws(fit, df)


# ---- Posterior ranks -------------------------------------------------------

ranks_pct <- posterior_ranks(theta, percentile = TRUE)
ranks_int <- posterior_ranks(theta, percentile = FALSE)

shifts <- rank_shifts(ranks_int, pub)
readr::write_csv(shifts, file.path(PATHS$results, "rank_shifts.csv"))

biggest <- shifts |>
  mutate(abs_delta = abs(delta)) |>
  slice_max(abs_delta, n = 15) |>
  select(benchmark, model_name, score, published_rank, rank_median, delta)

message("Largest rank movements after separating task, scaffold and benchmark variation:")
print(biggest, n = Inf)
save_table(biggest, "table_largest_rank_shifts", digits = 2)


# ---- Figure 6: published rank against posterior rank -----------------------

fig6_data <- ranks_pct |>
  inner_join(pub, by = c("benchmark", "model_name")) |>
  filter(benchmark %in% FEATURED_BENCHMARKS)

save_figure(plot_rank_comparison(fig6_data), "fig6_rank_comparison",
            width = 9.5, height = 6.5)

# The same panel for every benchmark (Appendix Figure 11).
all_ranks <- ranks_pct |> inner_join(pub, by = c("benchmark", "model_name"))
save_figure(plot_rank_comparison(all_ranks), "appendix_rank_comparison_all",
            width = 14, height = 10)


# ---- Figure 7: rank-shift slope chart --------------------------------------

save_figure(
  plot_rank_shift(filter(shifts, benchmark %in% FEATURED_BENCHMARKS)),
  "fig7_rank_shift", width = 10.5, height = 5.5
)
save_figure(
  plot_rank_shift(shifts, ncol = 3),
  "appendix_rank_shift_all", width = 15, height = 12
)


# ---- Table 2: agreement between observed and latent rankings ---------------
#
# Correlations are computed inside each posterior draw and then summarised, so
# the reported spread is uncertainty in the agreement itself.

agreement <- rank_agreement_draws(theta, pub)

table2 <- agreement |>
  filter(metric %in% c("spearman", "dcor2n")) |>
  mutate(metric = recode(metric, spearman = "Spearman", dcor2n = "dCor^2_n"),
         benchmark = BENCHMARK_LABELS[as.character(benchmark)]) |>
  select(benchmark, metric, median) |>
  tidyr::pivot_wider(names_from = benchmark, values_from = median)

print(table2)
save_table(table2, "table2_rank_agreement",
           caption = "Median correlation of observed and latent ranks across posterior draws.")
save_table(agreement, "table2_rank_agreement_full", digits = 3)


# ---- Full model versus per-benchmark estimator (Appendix D.3) --------------
#
# The per-benchmark fits use only within-benchmark information, and shrink
# harder as a result. Lower agreement with raw ranks is an estimand difference,
# not evidence that either estimator is wrong.

bench_fits  <- fit_benchmarks(df)
theta_bench <- latent_capability_per_benchmark(bench_fits, df) |>
  mutate(benchmark = factor(benchmark, levels = BENCHMARKS))

estimator_comparison <- bind_rows(
  rank_agreement_draws(theta, pub)       |> mutate(estimator = "Full (pooled, Eq. 6)"),
  rank_agreement_draws(theta_bench, pub) |> mutate(estimator = "Per-benchmark (Eq. 2)")
) |>
  select(estimator, benchmark, metric, median) |>
  tidyr::pivot_wider(names_from = benchmark, values_from = median) |>
  arrange(metric, estimator)

print(estimator_comparison, n = Inf)
save_table(estimator_comparison, "table_estimator_rank_agreement",
           caption = paste("Association between observed rankings and latent-effect",
                           "rankings, pooled versus per-benchmark estimator."))


# ---- Section 5.6: how many ranks are actually resolved? --------------------
#
# A high rank correlation does not imply that individual differences are
# distinguishable. Online-Mind2Web is the paper's example: it can correlate
# well with the published order while supporting almost no pairwise
# distinction under the posterior.

pairs <- pairwise_ordering(theta, threshold = 0.95)
readr::write_csv(pairs, file.path(PATHS$results, "pairwise_ordering.csv"))

indist <- pairs |>
  group_by(benchmark) |>
  summarize(n_pairs = n(), n_resolved = sum(resolved),
            prop_indistinguishable = 1 - mean(resolved), .groups = "drop") |>
  mutate(benchmark_label = BENCHMARK_LABELS[as.character(benchmark)], .before = 1)

print(indist, n = Inf)
save_table(select(indist, benchmark_label, n_pairs, n_resolved, prop_indistinguishable),
           "table_indistinguishable_pairs",
           caption = paste("Share of model pairs whose ordering the posterior does",
                           "not resolve at the 95% level."))

# Correlation with the published ordering set against the share of pairs that
# ordering actually resolves. High agreement and low resolution can coexist.
resolution_vs_agreement <- agreement |>
  filter(metric == "spearman") |>
  select(benchmark, spearman = median) |>
  inner_join(select(indist, benchmark, prop_indistinguishable), by = "benchmark")

print(resolution_vs_agreement, n = Inf)
save_table(resolution_vs_agreement, "table_agreement_vs_resolution")

p_res <- ggplot(resolution_vs_agreement,
                aes(spearman, prop_indistinguishable, label = BENCHMARK_LABELS[as.character(benchmark)])) +
  geom_point(size = 2.5, colour = COMPONENT_COLOURS[["M"]]) +
  ggrepel::geom_text_repel(size = 3) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Spearman correlation with published ranking",
       y = "Model pairs the posterior cannot order") +
  paper_theme(legend = "none")

save_figure(p_res, "agreement_vs_resolution", width = 7, height = 5)

message("\n04_rank_analysis.R complete.")
