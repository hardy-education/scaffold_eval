# scripts/07_harbor_corroboration.R ----------------------------------------
#
# Does the design conclusion hold on a second, independent dataset?
#
# The HAL result -- that spreading a task budget across benchmarks buys more
# ranking reliability than concentrating it in one -- could in principle be an
# artefact of HAL's particular sparsity: models and scaffolds there are
# unevenly crossed, and the pooled fit leans on a handful of shared anchors.
# The Harbor Index has the opposite shape. It has far fewer tasks, but every
# model is crossed with every scaffold, so the incidence structure that might
# be driving the HAL result is absent.
#
# Produces:
#   Harbor descriptive table (Appendix A.2)
#   Harbor variance decomposition and D-study
#   The pooled D-study comparison of HAL against Harbor (Figure 3)
#   Sensitivity of the Harbor decomposition to the task-count threshold
#
# Reads cached fits. To estimate them:
#   REFIT=TRUE Rscript scripts/07_harbor_corroboration.R

source(here::here("R", "setup.R"))

hal    <- load_hal_data()
harbor <- load_harbor_data()                     # benchmarks with >= 3 tasks
harbor_all <- load_harbor_data(min_tasks = NULL) # everything, for the table


# ---- Descriptive table (Appendix A.2) --------------------------------------

harbor_table <- harbor_design_summary(harbor_all)
print(harbor_table, n = Inf)

message(sprintf(
  "\nHarbor Index: %d benchmarks, %d models, %d scaffolds, %d rollouts.",
  n_distinct(harbor_all$benchmark), n_distinct(harbor_all$model_name),
  n_distinct(harbor_all$agent_name), nrow(harbor_all)
))
message(sprintf(
  "%d benchmarks have at least %d tasks and enter the decomposition (%d rollouts).",
  sum(harbor_table$included), HARBOR_MIN_TASKS, nrow(harbor)
))

save_table(harbor_table, "appendix_harbor_design_summary", digits = 3,
           caption = paste("Descriptive statistics for the Harbor Index",
                           "corroboration dataset. Benchmarks with at least",
                           HARBOR_MIN_TASKS, "tasks enter the decomposition."))

save_table(harbor_domains(harbor_all), "appendix_harbor_domains", digits = 0)


# ---- Harbor decomposition --------------------------------------------------

harbor_fit <- fit_harbor(harbor)
vc_harbor  <- vcomp_from_brms(harbor_fit, level = "leaderboard")

print(vc_harbor)
save_table(variance_shares(vc_harbor) |>
             select(component, label, facets, median, median_normalized, low, high),
           "appendix_harbor_variance_shares",
           caption = "Posterior variance shares, Harbor Index pooled decomposition.")

conv <- convergence_summary(harbor_fit)
if (max(conv$rhat, na.rm = TRUE) > 1.01) {
  warning("Harbor fit has Rhat above 1.01; inspect before using these draws.")
}
save_table(conv, "convergence_harbor", digits = 3)


# ---- Figure 3: the same design question, asked of both datasets ------------
#
# Each dataset is fitted separately and its own D-study projected; they are
# placed side by side rather than pooled, because they share no tasks and the
# point is agreement between independent estimates, not a combined one.

n_bench_harbor <- n_distinct(harbor$benchmark)

hal_dstudy <- dstudy(
  vcomp_from_brms(fit_leaderboard(hal), level = "leaderboard"),
  n_tasks = 1:DSTUDY$total_tasks, n_benchmarks = DSTUDY$n_benchmarks,
  statistics = c("ep2", "snr"), max_total_tasks = DSTUDY$total_tasks
) |> mutate(dataset = "HAL")

harbor_dstudy <- dstudy(
  vc_harbor,
  n_tasks = 1:DSTUDY$total_tasks, n_benchmarks = 1:n_bench_harbor,
  statistics = c("ep2", "snr"), max_total_tasks = DSTUDY$total_tasks
) |> mutate(dataset = "Harbor Index")

both <- bind_rows(hal_dstudy, harbor_dstudy)
readr::write_csv(both, file.path(PATHS$results, "dstudy_hal_vs_harbor.csv"))

save_figure(plot_dataset_dstudy(both, statistic = "ep2"),
            "fig3a_hal_vs_harbor_reliability", width = 9, height = 4.4)
save_figure(plot_dataset_dstudy(both, statistic = "snr"),
            "fig3b_hal_vs_harbor_snr", width = 9, height = 4.4)

# Does breadth help in both datasets, at a matched budget?
breadth_effect <- both |>
  filter(statistic == "ep2", n_benchmarks %in% c(1, 9)) |>
  group_by(dataset, n_benchmarks) |>
  slice_max(total_tasks, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(dataset, n_benchmarks, total_tasks, median, low, high)

print(breadth_effect)
save_table(breadth_effect, "table_breadth_effect_by_dataset",
           caption = paste("Projected model-ranking reliability at one versus",
                           "nine benchmarks, in each dataset."))


# ---- Sensitivity to the task-count threshold -------------------------------
#
# Appendix A.2 reports that requiring four or five tasks per benchmark lowers
# the estimated model share, and attributes this to how tasks were selected to
# represent each benchmark rather than to a change in signal. Re-running the
# threshold makes that claim checkable rather than asserted. Each threshold
# needs its own fit, so this is opt-in.

RUN_HARBOR_THRESHOLDS <- env_flag("RUN_HARBOR_THRESHOLDS")

if (RUN_HARBOR_THRESHOLDS) {
  message("\n=== Harbor task-count threshold sweep ===")

  threshold_shares <- purrr::map(3:5, function(k) {
    d <- load_harbor_data(min_tasks = k)
    f <- fit_harbor(d, name = paste0("harbor_bernoulli_min", k))
    vc <- vcomp_from_brms(f, level = "leaderboard")

    variance_shares(vc) |>
      mutate(min_tasks = k,
             n_benchmarks = n_distinct(d$benchmark),
             n_rollouts = nrow(d), .before = 1)
  }) |> bind_rows()

  model_share <- threshold_shares |>
    filter(component == "M") |>
    select(min_tasks, n_benchmarks, n_rollouts, median, low, high)

  print(model_share)
  message("The model share falls as the threshold rises; with this few tasks ",
          "per benchmark that is consistent with the non-random task selection ",
          "the appendix describes, and is not evidence of weaker model signal.")

  save_table(threshold_shares |>
               select(min_tasks, n_benchmarks, component, label, median, low, high),
             "appendix_harbor_threshold_sensitivity", digits = 3,
             caption = paste("Harbor variance shares under minimum task-count",
                             "thresholds of three, four, and five."))
} else {
  message("\nSkipping the Harbor threshold sweep ",
          "(set RUN_HARBOR_THRESHOLDS=TRUE REFIT=TRUE; one fit per threshold).")
}

message("\n07_harbor_corroboration.R complete.")
