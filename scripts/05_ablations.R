# scripts/05_ablations.R ---------------------------------------------------
#
# Is the decomposition driven by any single instrument, model, or condition?
#
# Three deletions of increasing structural severity (paper Appendix D.4-D.6):
#
#   leave-one-benchmark-out          removes an entire instrument, its tasks,
#                                    and all of its interactions. Full Bayesian
#                                    refits, because the reliability curves
#                                    themselves are the object of interest.
#   leave-one-model-out              removes one LLM and every cell it appears
#                                    in. Frequentist, 54 refits.
#   leave-one-benchmark-scaffold-out removes one (benchmark, scaffold) cell,
#                                    which for a benchmark-specific scaffold
#                                    deletes nearly all direct information
#                                    about that scaffold. Frequentist.
#
# The Bayesian arm needs REFIT=TRUE and roughly 57 hours; the frequentist arms
# run in tens of minutes and are enabled by default.
#
#   RUN_LOBO_BAYES=TRUE REFIT=TRUE Rscript scripts/05_ablations.R

source(here::here("R", "setup.R"))

df <- load_hal_data()

RUN_LOBO_BAYES <- env_flag("RUN_LOBO_BAYES")
N_CORES <- as.integer(Sys.getenv("N_CORES", max(1L, parallel::detectCores() - 2L)))


# ---- Shared helper ---------------------------------------------------------

#' Refit the leaderboard decomposition with lme4 on a subset, and return its
#' variance shares tagged by what was withheld.
ablate_lme4 <- function(data, removed, family = "binomial") {
  fit <- fit_lme4(data, FORMULA_LEADERBOARD, family = family, nagq = 0L)
  variance_shares(vcomp_from_lme4(fit, level = "leaderboard")) |>
    mutate(removed = removed, .before = 1)
}

#' Summarise stability of a component across refits (Tables 6 and 7).
stability_table <- function(shares) {
  shares |>
    group_by(component, label) |>
    summarize(
      mean = mean(estimate), median = stats::median(estimate),
      min = min(estimate), max = max(estimate),
      sd = stats::sd(estimate), n_refits = n(),
      se = sd / sqrt(n_refits),
      cv = sd / mean,
      .groups = "drop"
    ) |>
    arrange(desc(mean))
}


# ---- Leave-one-benchmark-out, Bayesian (Appendix D.4) ----------------------

if (RUN_LOBO_BAYES) {
  message("=== leave-one-benchmark-out, full Bayesian refits (~57 h) ===")

  lobo_vcs <- purrr::set_names(BENCHMARKS) |>
    purrr::map(\(b) {
      f <- fit_bayes(filter(df, benchmark != b), FORMULA_LEADERBOARD,
                     paste0("leaderboard_bernoulli_no_", b),
                     sampler = SAMPLER$loo)
      vcomp_from_brms(f, level = "leaderboard")
    })

  # Variance shares under each deletion.
  lobo_shares <- purrr::imap(lobo_vcs, \(vc, b) {
    variance_shares(vc) |> mutate(removed = b, .before = 1)
  }) |> bind_rows()

  save_table(lobo_shares |> select(removed, component, label, median, low, high),
             "ablation_lobo_variance_shares", digits = 3)

  # Reliability curves under each deletion. Withholding CORE-Bench Hard is the
  # case the paper singles out: it supplies both discrimination and
  # connectivity, so its removal widens uncertainty and weakens model signal.
  lobo_curves <- purrr::imap(lobo_vcs, \(vc, b) {
    dstudy(vc, n_tasks = 1:DSTUDY$total_tasks,
           n_benchmarks = c(1, 8), statistics = c("ep2", "snr"),
           max_total_tasks = DSTUDY$total_tasks) |>
      mutate(removed = b, .before = 1)
  }) |> bind_rows()

  readr::write_csv(lobo_curves, file.path(PATHS$results, "ablation_lobo_dstudy.csv"))

  save_figure(plot_loo_curves(lobo_curves, "ep2"),
              "ablation_lobo_reliability", width = 10, height = 4.5)
  save_figure(plot_loo_curves(lobo_curves, "snr"),
              "ablation_lobo_snr", width = 10, height = 4.5)

  # The model-vs-scaffold conclusion, recomputed under each deletion. The
  # paper reports the range of these probabilities alongside the headline.
  lobo_contrast <- purrr::imap(lobo_vcs, \(vc, b) {
    model_vs_scaffold_contrast(vc, design(1, 1, 1), scope = "both") |>
      mutate(removed = b, .before = 1)
  }) |> bind_rows()

  save_table(lobo_contrast |> select(removed, scope, label, median, low, high,
                                     p_scaffold_exceeds_model),
             "ablation_lobo_model_vs_scaffold", digits = 3)

  for (s in unique(lobo_contrast$scope)) {
    r <- range(lobo_contrast$p_scaffold_exceeds_model[lobo_contrast$scope == s])
    message(sprintf("P(scaffold > model), %s contrast, across LOBO refits: %.0f%% to %.0f%%",
                    s, 100 * r[1], 100 * r[2]))
  }

  # Which deletion costs the most reliability?
  ceilings <- purrr::imap_dfr(lobo_vcs, \(vc, b) {
    summarize_quantity(reliability_ceiling(vc, n_benchmarks = 8)) |>
      mutate(removed = b, .before = 1)
  }) |> arrange(median)
  message("\nTask-only ceiling (8 benchmarks) after each deletion, worst first:")
  print(select(ceilings, removed, median, low, high), n = Inf)
  save_table(select(ceilings, removed, median, low, high), "ablation_lobo_ceilings")
} else {
  message("Skipping Bayesian leave-one-benchmark-out refits ",
          "(set RUN_LOBO_BAYES=TRUE REFIT=TRUE to run; ~57 h).")
}


# ---- Leave-one-model-out (Appendix D.6) ------------------------------------
#
# A stronger perturbation than deleting one response: it removes a model main
# effect level and every benchmark-model, model-scaffold, task-model and
# benchmark-model-scaffold cell involving that LLM.

message("\n=== leave-one-model-out (", n_distinct(df$model_name), " refits) ===")

lom_shares <- parallel::mclapply(
  unique(df$model_name),
  \(m) ablate_lme4(filter(df, model_name != m), removed = m),
  mc.cores = N_CORES
) |> bind_rows()

lom_stability <- stability_table(lom_shares)
print(lom_stability, n = Inf)

save_table(lom_stability, "table7_leave_one_model_out",
           caption = paste("Leave-one-LLM-out proportions of variance.",
                           "SD is the standard deviation across refits;",
                           "CV is the coefficient of variation."))
readr::write_csv(lom_shares, file.path(PATHS$results, "ablation_leave_one_model.csv"))

save_figure(
  plot_ablation_stability(lom_shares |>
                            group_by(component) |>
                            mutate(share = estimate, share_mean = mean(estimate)) |>
                            ungroup()),
  "ablation_leave_one_model_stability", width = 7, height = 4.5
)


# ---- Leave-one-benchmark-scaffold-out (Appendix D.5) -----------------------
#
# Especially stringent for benchmark-specific scaffolds. Note that the
# components this destabilises most -- the scaffold main effect and the
# benchmark-scaffold interaction -- cancel from same-condition model contrasts
# and therefore do not enter relative model reliability at all.

message("\n=== leave-one-benchmark-scaffold-out ===")

bench_scaffold_cells <- df |> distinct(benchmark, agent_name)
message(nrow(bench_scaffold_cells), " (benchmark, scaffold) cells")

lobs_shares <- parallel::mclapply(
  seq_len(nrow(bench_scaffold_cells)),
  function(i) {
    cell <- bench_scaffold_cells[i, ]
    kept <- df |> filter(!(benchmark == cell$benchmark & agent_name == cell$agent_name))
    # A deletion that disconnects the design cannot support the decomposition.
    if (n_distinct(kept$benchmark) < 2 || n_distinct(kept$agent_name) < 2) return(NULL)
    tryCatch(
      ablate_lme4(kept, removed = paste(cell$benchmark, cell$agent_name, sep = " / ")),
      error = function(e) NULL
    )
  },
  mc.cores = N_CORES
) |> purrr::compact() |> bind_rows()

lobs_stability <- stability_table(lobs_shares)
print(lobs_stability, n = Inf)

save_table(lobs_stability, "table6_leave_one_benchmark_scaffold_out",
           caption = paste("Leave-one-benchmark-scaffold-out proportions of variance.",
                           "SE is the standard error across refits;",
                           "CV is the coefficient of variation."))
readr::write_csv(lobs_shares,
                 file.path(PATHS$results, "ablation_leave_one_bench_scaffold.csv"))

# Coefficients of variation are only interpretable next to the absolute scale:
# the task-model component is the most variable in relative terms and among
# the smallest in absolute terms.
message("\nComponents ranked by coefficient of variation (read alongside `mean`):")
print(lobs_stability |> arrange(desc(cv)) |> select(label, mean, min, max, cv), n = Inf)

message("\n05_ablations.R complete.")
