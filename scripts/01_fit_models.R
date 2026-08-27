# scripts/01_fit_models.R --------------------------------------------------
#
# Estimates the two decompositions the paper's results rest on:
#
#   Eq. 6   one pooled leaderboard model over all nine benchmarks
#   Eq. 2   one model within each benchmark
#
# This is the only expensive script. It refuses to run unless REFIT is set,
# because everything downstream reads the cached fits and a partial re-run
# would leave outputs/fits/ in a mixed state.
#
#   REFIT=TRUE Rscript scripts/01_fit_models.R
#
# Approximate wall-clock on an Apple M1 Max (Appendix C.1.3, C.6):
#
#   pooled leaderboard fit          ~5.5 h
#   nine benchmark-level fits       ~1 h total
#   convergence + LOO diagnostics   minutes (LOO with moment matching: ~38 h)

source(here::here("R", "setup.R"))

if (!REFIT) {
  stop(
    "01_fit_models.R re-estimates every model and takes hours.\n",
    "Run it as:  REFIT=TRUE Rscript scripts/01_fit_models.R\n",
    "The analysis scripts (02-06) read the cached fits and do not need this.",
    call. = FALSE
  )
}

df <- load_hal_data()

# ---- Pooled leaderboard decomposition (Eq. 6) ------------------------------

message("=== leaderboard fit (Eq. 6) ===")
leaderboard_fit <- fit_leaderboard(df)
print(summary(leaderboard_fit))

conv <- convergence_summary(leaderboard_fit)
print(conv, n = Inf)
save_table(conv, "convergence_leaderboard", digits = 3)

if (max(conv$rhat, na.rm = TRUE) > 1.01) {
  warning("Some Rhat exceeds 1.01; inspect the chains before using these draws.")
}


# ---- Benchmark-level decompositions (Eq. 2) --------------------------------

message("=== per-benchmark fits (Eq. 2) ===")
benchmark_fits <- fit_benchmarks(df)

conv_bench <- purrr::imap(benchmark_fits, \(f, b) {
  convergence_summary(f) |> mutate(benchmark = b, .before = 1)
}) |> bind_rows()

print(conv_bench, n = Inf)
save_table(conv_bench, "convergence_benchmark", digits = 3)


# ---- Sensitivity fits ------------------------------------------------------
#
# Two alternative specifications reported in the appendices. Both are optional:
# set FIT_SENSITIVITY to include them.

FIT_SENSITIVITY <- env_flag("FIT_SENSITIVITY")

if (FIT_SENSITIVITY) {
  # Appendix C.2: same structure with a Gaussian identity link, to show how
  # much of the decomposition depends on modelling the binary outcome as
  # binary rather than as a proportion.
  message("=== leaderboard fit, Gaussian link (Appendix C.2) ===")
  fit_bayes(df, FORMULA_LEADERBOARD, "leaderboard_gaussian",
            family = gaussian(), sampler = SAMPLER$leaderboard)

  # Appendix D.2: adding the four-way interaction, which acts as an
  # observation-level random effect in this design.
  message("=== leaderboard fit with four-way interaction (Appendix D.2) ===")
  fit_bayes(df, FORMULA_LEADERBOARD_FULL, "leaderboard_bernoulli_full",
            sampler = SAMPLER$leaderboard)
}


# ---- Predictive diagnostics (Appendix D.2) ---------------------------------

RUN_LOO <- env_flag("RUN_LOO")

if (RUN_LOO) {
  message("=== PSIS-LOO (slow; ~38 h with moment matching) ===")
  loo_obj <- psis_loo(leaderboard_fit, "leaderboard_bernoulli",
                      moment_match = TRUE)
  print(loo_obj)

  # High Pareto-k observations should be concentrated in singleton interaction
  # levels: cells whose only observation is the one being removed.
  infl <- influential_observations(loo_obj, df)
  message(sprintf(
    "%d observations with k > 0.7; %.1f%% of them are in singleton (b, i, m) groups.",
    nrow(infl), 100 * mean(infl$group_n == 1)
  ))
  save_table(select(infl, benchmark, task_id, model_name, agent_name,
                    score, pareto_k, group_n),
             "psis_loo_influential", digits = 3)
}

message("\n01_fit_models.R complete. Fits cached in ", PATHS$fits)
