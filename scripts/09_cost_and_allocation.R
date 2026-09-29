# scripts/09_cost_and_allocation.R -----------------------------------------
#
# What does reliability cost, and where should the budget go?
#
# The D-studies say which design interventions reduce which error. Attaching a
# price turns that into an allocation question: at a fixed budget, is it
# better to run more tasks on one benchmark or fewer tasks on several? The
# variance structure already answers the measurement half; this script adds
# the arithmetic.
#
# Produces:
#   Cost-reliability frontier over designs
#   Cheapest design reaching a target reliability, and the saving against the
#     observed battery
#   Repeated task subsampling as a design-free check on the same claim
#
# A caution carried through every output here: the response matrices record
# outcomes, not spend. The per-rollout price is an assumption set in `COST`
# in config.R, back-calculated so the full observed battery reproduces its
# reported total. Real cost varies by benchmark, model, and episode length, so
# treat dollar figures as planning projections whose ordering is more
# trustworthy than their levels. See the note in config.R for supplying a
# benchmark-specific price table.

source(here::here("R", "setup.R"))

df <- load_hal_data()
vc <- vcomp_from_brms(fit_leaderboard(df), level = "leaderboard")

baseline <- observed_cost(df)
message(sprintf(
  "Observed battery: %s rollouts at %s each = %s.",
  format(baseline$n_rollouts, big.mark = ","),
  scales::dollar(baseline$cost_per_trial, accuracy = 0.01),
  scales::dollar(baseline$total_cost)
))


# ---- Cost-reliability frontier ---------------------------------------------

frontier <- cost_reliability_frontier(
  vc,
  n_tasks = c(1, 2, 5, 10, 14, 20, 25, 30, 50, 75, 100, 150, 200),
  n_benchmarks = 1:9,
  n_scaffolds = DSTUDY$n_scaffolds
)

readr::write_csv(frontier, file.path(PATHS$results, "cost_reliability_frontier.csv"))

# At a matched budget, how much does breadth buy?
matched <- frontier |>
  filter(total_tasks %in% c(90, 126, 180, 270, 450)) |>
  select(n_benchmarks, n_tasks, total_tasks, trials, cost, median, low, high) |>
  arrange(total_tasks, n_benchmarks)

print(matched, n = Inf)
save_table(matched, "table_cost_at_matched_budgets", digits = 3,
           caption = paste("Projected model-ranking reliability and cost for",
                           "designs at matched total task budgets."))

save_figure(plot_cost_frontier(frontier), "fig_cost_reliability_frontier",
            width = 8, height = 5)


# ---- Cheapest design reaching a target -------------------------------------

# Targets are chosen to sit under the battery's own asymptote. The
# nine-benchmark ceiling is ~0.75 in the infinite-task limit, so asking the
# frontier for 0.75 or more returns nothing and says so -- which is itself the
# point of Proposition 3.1, not a gap in the grid.
targets <- c(0.6, 0.65, 0.70, 0.715)

cheapest <- purrr::map(targets, \(t) {
  cheapest_design_reaching(frontier, t, baseline_cost = baseline$total_cost)
}) |> bind_rows()

if (nrow(cheapest)) {
  print(select(cheapest, target, n_benchmarks, n_tasks, total_tasks,
               trials, cost, saving, median), n = Inf)

  for (i in seq_len(nrow(cheapest))) {
    r <- cheapest[i, ]
    message(sprintf(
      "Erho^2 >= %.2f: %d benchmarks x %d tasks = %d tasks, %s (%.0f%% below the observed battery).",
      r$target, r$n_benchmarks, r$n_tasks, r$total_tasks,
      scales::dollar(r$cost), 100 * r$saving
    ))
  }

  save_table(select(cheapest, target, n_benchmarks, n_tasks, total_tasks,
                    trials, cost, baseline_cost, saving, median, low, high),
             "table_cheapest_design_by_target", digits = 3,
             caption = paste("Cheapest projected design reaching each",
                             "model-ranking reliability target."))
}

# Depth against breadth, stated as sharply as the data allow: take the best
# any single benchmark can do at any budget in the grid, then find the
# cheapest broad design that matches it.
concentrated <- frontier |>
  filter(n_benchmarks == 1) |>
  slice_max(median, n = 1, with_ties = FALSE)

broad_equivalent <- frontier |>
  filter(n_benchmarks == max(n_benchmarks), median >= concentrated$median) |>
  slice_min(cost, n = 1, with_ties = FALSE)

message(sprintf(
  "\nDepth vs. breadth:\n  the best single-benchmark design in the grid (%d tasks) reaches Erho^2 = %.3f for %s.",
  concentrated$n_tasks, concentrated$median, scales::dollar(concentrated$cost)
))
if (nrow(broad_equivalent)) {
  message(sprintf(
    "  %d benchmarks match it with %s per benchmark, for %s -- %.0f%% of the cost.",
    broad_equivalent$n_benchmarks,
    if (broad_equivalent$n_tasks == 1) "1 task" else paste(broad_equivalent$n_tasks, "tasks"),
    scales::dollar(broad_equivalent$cost),
    100 * broad_equivalent$cost / concentrated$cost
  ))
  message("  Read this as a statement about where the error lives, not a ",
          "recommendation to run one task per benchmark: a design that thin ",
          "leaves nothing to diagnose individual results with.")
} else {
  message("  No design at full breadth in this grid matches it; widen the grid.")
}


# ---- Repeated task subsampling ---------------------------------------------
#
# The frontier is a projection from the fitted variance components. This is
# the same question asked without a model: sample k tasks per benchmark, use
# the same tasks for every model so the comparison is paired, and see how well
# the resulting ranking matches the full-data one.
#
# Agreement establishes that the information survives the cut. It says nothing
# about whether the full-data ranking was itself well resolved -- scripts/04
# addresses that, and finds many pairs unresolved.

N_SUBSAMPLES <- as.integer(Sys.getenv("N_SUBSAMPLES", "500"))
message(sprintf("\n=== Repeated task subsampling (%d replicates per k) ===",
                N_SUBSAMPLES))

subsamples <- subsample_rankings(df, k_values = c(5, 10, 15, 20, 25, 30),
                                 n_rep = N_SUBSAMPLES)
readr::write_csv(subsamples, file.path(PATHS$results, "task_subsampling_draws.csv"))

subsample_table <- summarize_subsampling(subsamples, df)
print(subsample_table, n = Inf)

save_table(subsample_table, "table_task_subsampling", digits = 3,
           caption = paste("Repeated task subsampling: rank agreement with the",
                           "complete-data ranking as a function of tasks",
                           "retained per benchmark."))

save_figure(plot_subsampling(subsamples), "fig_task_subsampling",
            width = 8, height = 4.5)

message(
  "\nThe cost column is rollout-count weighted under the flat default price. ",
  "Benchmarks differ in per-rollout cost as well as task count, so a ",
  "benchmark-specific price table would shift these figures; supply one via ",
  "COST$per_trial_by_benchmark in config.R."
)

message("\n09_cost_and_allocation.R complete.")
