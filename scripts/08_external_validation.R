# scripts/08_external_validation.R -----------------------------------------
#
# Does the latent model effect transport outside the panel it was fitted on?
#
# The pooled fit splits a model's performance into a shared component that
# should persist across benchmarks (theta_m) and benchmark-conditioned
# departures that should not. If that split is real rather than an artefact of
# these nine benchmarks, theta_m should predict performance on benchmarks the
# fit never saw -- and should do so better than the conventional leaderboard
# aggregate, which credits a model for every condition it happened to be
# paired with.
#
# Produces:
#   Agreement of each in-panel estimator with four external benchmarks (Table 3)
#   Paired bootstrap of the difference in Kendall's tau, with leave-one-model-out
#   Kendall's tau matrix among in-panel benchmarks and the two aggregates
#
# Contamination control: two external benchmarks overlap a HAL benchmark in
# content, so both estimators are rebuilt without it. That needs the matching
# leave-one-benchmark-out refits:
#
#   RUN_LOBO_BAYES=TRUE REFIT=TRUE Rscript scripts/05_ablations.R
#
# The script stops with an explanation rather than silently falling back to
# the full-data fit, because that fallback would inflate exactly the number
# the section is about.

source(here::here("R", "setup.R"))

df  <- load_hal_data()
fit <- fit_leaderboard(df)
edf <- load_external_data()

message(sprintf("External validation data: %d benchmarks, %d models, %d rows.",
                n_distinct(edf$benchmark), n_distinct(edf$model_name), nrow(edf)))

# Provenance is part of the claim here: these are third-party leaderboards,
# not reruns, so it matters which source each score came from.
provenance <- edf |>
  group_by(benchmark, leaderboard) |>
  summarize(n_models = n_distinct(model_name),
            scaffold = dplyr::first(agent_scaffold_harness),
            metric = dplyr::first(metric), .groups = "drop") |>
  mutate(benchmark_label = EXTERNAL_LABELS[as.character(benchmark)], .before = 1)

print(provenance, n = Inf)
save_table(provenance, "appendix_external_provenance", digits = 0,
           caption = "Sources of the external validation scores.")


# ---- Load only the leave-one-out refits this analysis needs ----------------

needed <- unique(stats::na.omit(EXTERNAL_EXCLUSIONS))
message("\nContamination control requires leave-one-benchmark-out refits for: ",
        paste(needed, collapse = ", "))

loo_fits <- purrr::set_names(needed) |>
  purrr::map(\(b) fit_bayes(filter(df, benchmark != b), FORMULA_LEADERBOARD,
                            paste0("leaderboard_bernoulli_no_", b),
                            sampler = SAMPLER$loo))


# ---- Table 3: agreement with external benchmarks ---------------------------

comparison <- external_comparison_data(df, fit, edf, loo_fits)
readr::write_csv(comparison, file.path(PATHS$results, "external_comparison.csv"))

agreement <- external_agreement(comparison)
print(agreement, n = Inf)

table3 <- agreement |>
  select(benchmark_label, estimator, tau, low, high, n_models) |>
  mutate(across(c(tau, low, high), \(x) round(x, 3))) |>
  tidyr::pivot_wider(names_from = estimator,
                     values_from = c(tau, low, high)) |>
  mutate(difference = `tau_theta_hat` - `tau_Mean score`)

print(table3, n = Inf)
save_table(table3, "table3_external_agreement",
           caption = paste("Agreement with external agent benchmarks.",
                           "Kendall's tau for the unweighted HAL mean score and",
                           "the reliability-adjusted latent model effect."))

message(sprintf(
  "\nMean Kendall's tau across the four external benchmarks: %.3f (mean score) vs %.3f (theta_hat).",
  mean(agreement$tau[agreement$estimator == "Mean score"]),
  mean(agreement$tau[agreement$estimator == "theta_hat"])
))


# ---- Is the difference robust? ---------------------------------------------
#
# The claim is comparative, so the uncertainty that matters is on the paired
# difference. With six to twenty overlapping models, a single model can carry
# a rank correlation, so leave-one-model-out deltas are reported alongside.

delta <- delta_tau_bootstrap(comparison)
print(delta, n = Inf)

save_table(delta, "table_external_delta_tau", digits = 3,
           caption = paste("Paired bootstrap of the difference in Kendall's tau",
                           "between the latent model effect and the mean-score",
                           "aggregate, with leave-one-model-out sensitivity."))

for (i in seq_len(nrow(delta))) {
  r <- delta[i, ]
  message(sprintf(
    "%-20s delta tau = %+.3f  [95%% %.3f, %.3f]  n = %d  LOO range [%+.3f, %+.3f]  sign flips: %d",
    r$benchmark_label, r$delta, r$boot_low, r$boot_high, r$n_models,
    r$loo_min_delta, r$loo_max_delta, r$loo_sign_flips
  ))
}

if (any(delta$loo_sign_flips > 0)) {
  message("\nAt least one external benchmark changes the sign of the difference ",
          "when a single model is dropped. Treat that comparison as unresolved.")
}

save_figure(plot_delta_tau(delta), "fig_external_delta_tau",
            width = 7.5, height = 4)
save_figure(plot_external_scatter(comparison), "fig_external_scatter",
            width = 10, height = 6)


# ---- In-panel agreement, for context ---------------------------------------
#
# Benchmarks inside the panel agree with each other only moderately. That is
# the same heterogeneity the decomposition attributes to benchmark-conditioned
# model effects, seen directly in the raw rankings.

tau_matrix <- benchmark_rank_agreement(df, fit)
print(tau_matrix, n = Inf)
save_table(tau_matrix, "appendix_benchmark_tau_matrix",
           caption = paste("Kendall's tau among benchmark rankings, the",
                           "unweighted mean score, and the latent model effect."))

message("\n08_external_validation.R complete.")
