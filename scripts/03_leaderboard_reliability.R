# scripts/03_leaderboard_reliability.R -------------------------------------
#
# What the whole battery can tell you, and where its ceiling comes from.
#
# Produces:
#   Figure 4a  reliability against total task budget, 1 benchmark vs. 9
#   Figure 4b  signal-to-noise by battery size, against the detection band
#   Figure 5   the pooled variance decomposition
#   Section 5.1  posterior probability that scaffold variation exceeds model
#                variation, at benchmark scope and at task scope
#
# Reads the cached pooled fit from scripts/01_fit_models.R.
#
# The organising result: item count and evidence diversity are not
# interchangeable. Five hundred tasks drawn from one benchmark and 55 drawn
# from each of nine cost the same and do not buy the same thing.

source(here::here("R", "setup.R"))

df  <- load_hal_data()
fit <- fit_leaderboard(df)              # cached; no estimation
vc  <- vcomp_from_brms(fit, level = "leaderboard")

print(vc)


# ---- Figure 5: pooled variance decomposition -------------------------------
#
# Persistent model variation is smaller than the combined variation associated
# with scaffolds, benchmark-conditioned model performance, and task-conditioned
# interactions. Most of what an observed outcome contains is
# condition-specific rather than transferable model signal.

shares <- variance_shares(vc)
print(arrange(shares, desc(median)), n = Inf)

save_table(select(shares, component, label, facets, median, median_normalized, mean, low, high),
           "fig5_leaderboard_variance_shares",
           caption = "Posterior variance shares from the leaderboard decomposition (Eq. 6).")

save_figure(plot_variance_decomposition(shares),
            "fig5_leaderboard_variance_decomposition", width = 15.6, height = 2.66)

message(sprintf(
  "\nModel main effect carries %.1f%% of latent variance; the benchmark-model, model-scaffold\nand benchmark-model-scaffold terms together carry %.1f%%.",
  100 * shares$median[shares$component == "M"],
  100 * sum(shares$median[shares$component %in% c("BM", "MA", "BMA")])
))


# ---- Figure 4a: reliability against evaluation budget ----------------------
#
# Sweeping tasks-per-benchmark for each battery size, then plotting against
# total tasks, puts every design on a common budget axis.

leaderboard_dstudy <- dstudy(
  vc,
  n_tasks      = 1:DSTUDY$total_tasks,
  n_benchmarks = DSTUDY$n_benchmarks,
  n_scaffolds  = DSTUDY$n_scaffolds,
  statistics   = c("ep2", "snr"),
  max_total_tasks = DSTUDY$total_tasks
)

readr::write_csv(leaderboard_dstudy, file.path(PATHS$results, "leaderboard_dstudy.csv"))

save_figure(
  plot_leaderboard_reliability(filter(leaderboard_dstudy, statistic == "ep2")),
  "fig4a_leaderboard_reliability", width = 5.5, height = 4.4
)

save_figure(
  plot_leaderboard_snr(filter(leaderboard_dstudy, statistic == "snr")),
  "fig4b_leaderboard_snr", width = 6.2, height = 4.4
)


# ---- Ceilings: what breadth buys, and what it does not ---------------------
#
# Proposition 3.1 in numbers. Adding tasks removes item-indexed error; adding
# benchmarks attenuates benchmark-conditioned model differences; neither
# touches the model-scaffold coupling.

ceiling_grid <- tidyr::expand_grid(
  n_benchmarks = c(1, 3, 5, 7, 9, Inf),
  n_scaffolds  = c(1, 2, 3, Inf)
) |>
  mutate(.s = purrr::map2(n_benchmarks, n_scaffolds,
                          \(nb, na) summarize_quantity(reliability_ceiling(vc, nb, na)))) |>
  tidyr::unnest(.s)

print(select(ceiling_grid, n_benchmarks, n_scaffolds, median, low, high), n = Inf)
save_table(select(ceiling_grid, n_benchmarks, n_scaffolds, median, low, high),
           "table_leaderboard_ceilings",
           caption = paste("Task-only reliability ceilings for model ranking as",
                           "the number of benchmarks and scaffolds varies."))

headline <- function(nb, na) {
  round(summarize_quantity(reliability_ceiling(vc, nb, na))$median, 3)
}
message("\nTask-only ceilings for model ranking (Eq. 8):")
message(sprintf("  1 benchmark,  1 scaffold   Erho^2 -> %.2f", headline(1, 1)))
message(sprintf("  9 benchmarks, 1 scaffold   Erho^2 -> %.2f", headline(9, 1)))
message(sprintf("  unlimited benchmarks       Erho^2 -> %.2f", headline(Inf, 1)))
message(sprintf("  unlimited benchmarks and scaffolds  Erho^2 -> %.2f", headline(Inf, Inf)))


# ---- How many benchmarks reach the detection band? -------------------------

detection <- tidyr::expand_grid(n_benchmarks = 1:9, total_tasks = c(100, 250, 500)) |>
  mutate(
    n_tasks = total_tasks / n_benchmarks,
    .s = purrr::map2(n_tasks, n_benchmarks,
                     \(ni, nb) summarize_quantity(snr(vc, design(ni, nb, 1))))
  ) |>
  tidyr::unnest(.s) |>
  select(n_benchmarks, total_tasks, n_tasks, median, low, high) |>
  mutate(reaches_band = median >= LOD_BAND[["mid"]])

print(detection, n = Inf)
save_table(detection, "table_detection_by_battery_size",
           caption = paste("Signal-to-noise by battery size and total task",
                           "budget, against the S/N = 2.5 reference."))

smallest <- detection |>
  filter(total_tasks == 500, reaches_band) |>
  slice_min(n_benchmarks, n = 1, with_ties = FALSE)
if (nrow(smallest)) {
  message(sprintf(
    "\nAt a 500-task budget, %d benchmarks are needed to reach S/N = %.1f.",
    smallest$n_benchmarks, LOD_BAND[["mid"]]
  ))
}


# ---- Model variation against scaffold variation -----------------------------
#
# The paper makes a two-part claim here, and the parts point different ways:
#
#   main effects       sigma^2_M vs sigma^2_A. Do models differ more than
#                      scaffolds do on average? The posterior favours neither.
#   task interactions  sigma^2_IM[B] vs sigma^2_IA[B]. Does the scaffold change
#                      *which tasks* get solved more than the model does? Yes.
#
# Read together: a scaffold can change which tasks a system solves without
# making it uniformly stronger on the benchmark. Reporting only the first
# contrast would understate how much the harness matters; reporting only the
# second would overstate it.

contrasts <- model_vs_scaffold_contrast(vc, design(1, 1, 1), scope = "both")

for (i in seq_len(nrow(contrasts))) {
  message(sprintf(
    "P(scaffold > model), %-18s %.1f%%  [contrast median %.3f, %d%% HDI %.3f to %.3f]",
    contrasts$label[i], 100 * contrasts$p_scaffold_exceeds_model[i],
    contrasts$median[i], round(100 * CI_LEVEL),
    contrasts$low[i], contrasts$high[i]
  ))
}

save_table(select(contrasts, scope, label, left, right, median, mean, low, high,
                  p_model_exceeds_scaffold, p_scaffold_exceeds_model),
           "table_model_vs_scaffold_contrast",
           caption = paste("Posterior contrasts between model-related and",
                           "scaffold-related variance components: main effects",
                           "and task-indexed interactions."))

for (i in seq_len(nrow(contrasts))) {
  save_figure(plot_model_scaffold_contrast(contrasts[i, ]),
              paste0("model_vs_scaffold_contrast_", contrasts$scope[i]),
              width = 7, height = 4)
}

message("\n03_leaderboard_reliability.R complete.")
