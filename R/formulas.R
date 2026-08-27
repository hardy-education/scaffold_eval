# R/formulas.R -------------------------------------------------------------
#
# The two random-effects structures estimated in the paper, plus the reduced
# specifications used for method comparison.
#
# Notation, matching the paper: b = benchmark, i = task (item, nested within
# benchmark), m = model (LLM x reasoning effort), a = agent scaffold.
#
# Every term is a mean-zero Gaussian random intercept, so each corresponds to
# exactly one variance component in the decomposition. There are no fixed
# effects beyond the intercept: that is what makes the fit a G-study.

#' Leaderboard-level decomposition (paper Eq. 6 / Eq. 9).
#'
#' Tasks are nested within benchmarks; models and scaffolds are crossed with
#' benchmarks to the extent the observed incidence graph supports. Fitted
#' jointly to all nine benchmarks, so that models and scaffolds appearing in
#' more than one benchmark carry information across them.
FORMULA_LEADERBOARD <- score ~ 1 +
  (1 | benchmark) +                                    # B
  (1 | benchmark:task_id) +                            # I[B]
  (1 | model_name) +                                   # M   <- universe score
  (1 | agent_name) +                                   # A
  (1 | benchmark:model_name) +                         # BM
  (1 | benchmark:agent_name) +                         # BA
  (1 | model_name:agent_name) +                        # MA
  (1 | benchmark:task_id:model_name) +                 # IM[B]
  (1 | benchmark:task_id:agent_name) +                 # IA[B]
  (1 | benchmark:model_name:agent_name)                # BMA

#' Leaderboard decomposition including the four-way interaction.
#'
#' Under 5% of (b, i, m, a) cells are replicated, so this term acts as an
#' observation-level random effect competing with the fixed logistic residual
#' rather than identifying a genuine four-way interaction (Appendix D.2). It is
#' fitted only as a sensitivity check; the primary fit uses
#' `FORMULA_LEADERBOARD` and folds this variation into the terminal component.
FORMULA_LEADERBOARD_FULL <- update(
  FORMULA_LEADERBOARD,
  . ~ . + (1 | benchmark:task_id:model_name:agent_name)
)

#' Benchmark-level decomposition (paper Eq. 2).
#'
#' Fitted separately within each benchmark, so benchmark terms drop out. With
#' only two or three scaffolds per benchmark these fits cannot estimate
#' scaffold variation precisely -- which is exactly why the pooled fit exists.
FORMULA_BENCHMARK <- score ~ 1 +
  (1 | task_id) +                                      # I
  (1 | model_name) +                                   # M   <- universe score
  (1 | agent_name) +                                   # A
  (1 | task_id:model_name) +                           # IM
  (1 | task_id:agent_name) +                           # IA
  (1 | model_name:agent_name)                          # MA

#' Benchmark-level decomposition with model and scaffold collapsed.
#'
#' Treats the deployable (model, scaffold) system as a single facet. Used in
#' the method comparison of Appendix G.5 (Decomposition 2).
FORMULA_BENCHMARK_SYSTEM <- score ~ 1 +
  (1 | task_id) +
  (1 | model_system)                                   # M:A as one facet

#' Benchmark-level decomposition with no scaffold facet at all.
#'
#' Decomposition 1 of Appendix G.5: what a leaderboard implicitly assumes when
#' it reports one score per model and ignores the harness.
FORMULA_BENCHMARK_NO_SCAFFOLD <- score ~ 1 +
  (1 | task_id) +
  (1 | model_name)

#' Add the collapsed model-scaffold column required by `FORMULA_BENCHMARK_SYSTEM`.
add_model_system <- function(df) {
  mutate(df, model_system = paste(model_name, agent_name, sep = "@"))
}
