# R/variance.R -------------------------------------------------------------
#
# Turning a fitted model into variance components.
#
# Every reliability number in the paper is a ratio of variance components, so
# this file defines the one canonical representation those ratios are computed
# from: a `vcomp` object holding, for each facet or interaction, the variance
# on the latent log-odds scale.
#
# The same representation is produced from a Bayesian fit (components are
# `posterior::rvar` objects, carrying the full posterior) and from a
# frequentist lme4 fit (components are plain numbers). Downstream code in
# gtheory.R is written in ordinary arithmetic and therefore works with either,
# which is why the paper's Bayesian results and its lme4 sensitivity analyses
# share one implementation instead of two.

# ---- Component registry ---------------------------------------------------
#
# `facets` is the string of facets a component is indexed by, and is what the
# G-theory algebra in gtheory.R reasons over:
#
#   B = benchmark, I = task (item), M = model, A = agent scaffold
#
# Tasks are nested within benchmarks, so at the leaderboard level every
# item-indexed component carries a "B" as well: the number of distinct items
# sampled is n_benchmarks * n_tasks_per_benchmark.
#
# `terminal = TRUE` marks the component that absorbs unresolved cell-level
# variation. In a Bernoulli-logit fit without replicated (b, i, m, a) cells,
# that is the fixed logistic residual pi^2 / 3 plus any four-way interaction
# (paper Appendix E.4: sigma^2_{BIMA,e}).

#' Component registry for a given decomposition level.
#'
#' @param level "leaderboard" (Eq. 6) or "benchmark" (Eq. 2).
#' @return A tibble with one row per variance component.
component_registry <- function(level = c("leaderboard", "benchmark")) {
  level <- match.arg(level)

  if (level == "leaderboard") {
    tibble::tribble(
      ~component, ~facets,  ~label,                      ~brms_par,                                              ~lme4_grp,
      "B",        "B",      "Benchmark",                 "sd_benchmark__Intercept",                              "benchmark",
      "I",        "BI",     "Task",                      "sd_benchmark:task_id__Intercept",                      "benchmark_task_id",
      "M",        "M",      "Model",                     "sd_model_name__Intercept",                             "model_name",
      "A",        "A",      "Scaffold",                  "sd_agent_name__Intercept",                             "agent_name",
      "BM",       "BM",     "Benchmark-Model",           "sd_benchmark:model_name__Intercept",                   "benchmark_model_name",
      "BA",       "BA",     "Benchmark-Scaffold",        "sd_benchmark:agent_name__Intercept",                   "benchmark_agent_name",
      "MA",       "MA",     "Model-Scaffold",            "sd_model_name:agent_name__Intercept",                  "model_name_agent_name",
      "IM",       "BIM",    "Task-Model",                "sd_benchmark:task_id:model_name__Intercept",           "benchmark_task_id_model_name",
      "IA",       "BIA",    "Task-Scaffold",             "sd_benchmark:task_id:agent_name__Intercept",           "benchmark_task_id_agent_name",
      "BMA",      "BMA",    "Benchmark-Model-Scaffold",  "sd_benchmark:model_name:agent_name__Intercept",        "benchmark_model_name_agent_name",
      "IMA",      "BIMA",   "Task-Model-Scaffold",       "sd_benchmark:task_id:model_name:agent_name__Intercept","benchmark_task_id_model_name_agent_name",
      "e",        "BIMA",   "Residual",                  "sigma",                                                "Residual"
    ) |>
      mutate(terminal = component == "e")
  } else {
    tibble::tribble(
      ~component, ~facets, ~label,             ~brms_par,                              ~lme4_grp,
      "I",        "I",     "Task",             "sd_task_id__Intercept",                "task_id",
      "M",        "M",     "Model",            "sd_model_name__Intercept",             "model_name",
      "A",        "A",     "Scaffold",         "sd_agent_name__Intercept",             "agent_name",
      "MA",       "MA",    "Model-Scaffold",   "sd_model_name:agent_name__Intercept",  "model_name_agent_name",
      "IM",       "IM",    "Task-Model",       "sd_task_id:model_name__Intercept",     "task_id_model_name",
      "IA",       "IA",    "Task-Scaffold",    "sd_task_id:agent_name__Intercept",     "task_id_agent_name",
      "e",        "IMA",   "Residual",         "sigma",                                "Residual"
    ) |>
      mutate(terminal = component == "e")
  }
}

# Logistic residual variance: the latent-scale error floor implied by a
# Bernoulli-logit likelihood (paper Section 3.3).
LOGISTIC_RESIDUAL_VAR <- pi^2 / 3


# ---- Constructing vcomp objects -------------------------------------------

#' Assemble a `vcomp` object.
#'
#' @param values Named list of variances, keyed by component code.
#' @param level "leaderboard" or "benchmark".
#' @param source Free-text provenance label, e.g. "brms (bernoulli)".
new_vcomp <- function(values, level, source = NA_character_) {
  reg <- component_registry(level)
  present <- reg$component[reg$component %in% names(values)]
  missing <- setdiff(reg$component, present)

  # Components absent from a fit (for example the four-way interaction, which
  # the primary specification omits) contribute zero variance. Their role is
  # taken by the terminal component.
  for (m in missing) values[[m]] <- 0

  structure(
    list(
      values   = values[reg$component],
      registry = reg,
      level    = level,
      source   = source,
      fitted   = present
    ),
    class = "vcomp"
  )
}

#' @export
print.vcomp <- function(x, ...) {
  cat("<vcomp>", x$level, "decomposition", if (!is.na(x$source)) paste0("[", x$source, "]"), "\n")
  print(variance_shares(x), n = Inf)
  invisible(x)
}


#' Variance components from a Bayesian (brms) fit.
#'
#' Components are returned as `posterior::rvar` objects, so every downstream
#' reliability quantity is evaluated once per posterior draw. The paper's
#' summaries are therefore posterior distributions of ratios, not ratios of
#' posterior means -- which matters, because these ratios are nonlinear and
#' several have strongly skewed posteriors.
#'
#' For a Bernoulli-logit fit brms reports no `sigma`; the latent residual is
#' the fixed logistic variance pi^2 / 3, which this function supplies.
#'
#' @param fit A `brmsfit`.
#' @param level "leaderboard" or "benchmark".
vcomp_from_brms <- function(fit, level = c("leaderboard", "benchmark")) {
  level <- match.arg(level)
  reg <- component_registry(level)

  draws <- posterior::as_draws_rvars(fit)
  fit_family <- fit$family$family
  is_bernoulli <- identical(fit_family, "bernoulli")

  values <- list()
  for (k in seq_len(nrow(reg))) {
    par <- reg$brms_par[k]
    comp <- reg$component[k]

    if (comp == "e" && is_bernoulli) {
      # Fixed logistic residual on the latent scale.
      values[[comp]] <- posterior::rvar(
        rep(LOGISTIC_RESIDUAL_VAR, posterior::ndraws(draws))
      )
    } else if (!is.null(draws[[par]])) {
      values[[comp]] <- draws[[par]]^2
    }
  }

  new_vcomp(
    values, level,
    source = paste0("brms (", fit_family, ")")
  )
}


#' Variance components from a frequentist lme4 fit.
#'
#' Used for the estimation-method sensitivity analyses (Appendix C.3-C.4) and
#' for the leave-one-out ablations, where thousands of refits make full
#' Bayesian estimation impractical.
#'
#' @param fit An `lmerMod` or `glmerMod`.
#' @param level "leaderboard" or "benchmark".
#' @param add_logistic_residual Append pi^2 / 3 as the residual. Required for
#'   `glmer` fits, which report no residual variance.
vcomp_from_lme4 <- function(fit,
                            level = c("leaderboard", "benchmark"),
                            add_logistic_residual = inherits(fit, "glmerMod")) {
  level <- match.arg(level)
  reg <- component_registry(level)

  vc <- tibble::as_tibble(lme4::VarCorr(fit)) |>
    filter(is.na(var2)) |>
    mutate(grp = gsub(":", "_", grp)) |>
    select(grp, vcov)

  if (add_logistic_residual && !"Residual" %in% vc$grp) {
    vc <- bind_rows(vc, tibble(grp = "Residual", vcov = LOGISTIC_RESIDUAL_VAR))
  }

  values <- list()
  for (k in seq_len(nrow(reg))) {
    hit <- vc$vcov[vc$grp == reg$lme4_grp[k]]
    if (length(hit) == 1) values[[reg$component[k]]] <- hit
  }

  new_vcomp(
    values, level,
    source = paste0("lme4 (", if (inherits(fit, "glmerMod")) "binomial" else "gaussian", ")")
  )
}


# ---- Summaries ------------------------------------------------------------

#' Proportion of total latent variance carried by each component.
#'
#' This is the quantity plotted in the paper's variance-decomposition figure
#' and tabulated in Appendix G.4. Shares are computed within each posterior
#' draw before being summarised.
#'
#' @param vc A `vcomp`.
#' @param ci Credible-interval mass for the HDI (Bayesian input only).
variance_shares <- function(vc, ci = CI_LEVEL) {
  total <- Reduce(`+`, vc$values)
  out <- vc$registry |>
    mutate(share = purrr::map(component, \(k) summarize_quantity(vc$values[[k]] / total, ci = ci))) |>
    select(component, label, facets, share) |>
    tidyr::unnest(share)

  # Shares sum to one within each posterior draw, but their marginal medians do
  # not: the median of a sum is not the sum of medians. `median` is the right
  # number to quote for a single component; `median_normalized` rescales the
  # medians to sum to one so that a stacked bar chart is well formed. Never
  # quote the normalized value as an estimate.
  mutate(out, median_normalized = median / sum(median))
}


#' Summarise a scalar quantity that may be a posterior draw or a point estimate.
#'
#' Bayesian input returns posterior median, mean, HDI bounds at `ci`, 5th/95th
#' percentiles, MAP, and split-Rhat. Point-estimate input returns the value
#' with `NA` uncertainty, so that Bayesian and frequentist results can share
#' plotting and table code.
#'
#' The paper reports medians rather than means for derived quantities, because
#' ratios of variance components are skewed.
summarize_quantity <- function(x, ci = CI_LEVEL) {
  if (!inherits(x, "rvar")) {
    return(tibble(
      estimate = as.numeric(x), median = as.numeric(x), mean = as.numeric(x),
      sd = NA_real_, low = NA_real_, high = NA_real_,
      q5 = NA_real_, q95 = NA_real_, map = NA_real_, rhat = NA_real_
    ))
  }

  d <- as.numeric(posterior::draws_of(x))
  s <- posterior::summarise_draws(x, "mean", "median", "sd", "rhat")

  # Some designs make a quantity degenerate rather than uncertain: at
  # n_tasks = Inf the relative error of a system-object coefficient is exactly
  # zero, so every draw equals 1. Density-based summaries are undefined there,
  # and the honest answer is the constant itself.
  finite <- d[is.finite(d)]
  degenerate <- length(finite) < 2 || stats::sd(finite) < .Machine$double.eps^0.5

  if (degenerate) {
    v <- if (length(finite)) finite[1] else NA_real_
    return(tibble(
      estimate = v, median = v, mean = v, sd = 0,
      low = v, high = v, q5 = v, q95 = v, map = v, rhat = s$rhat
    ))
  }

  hdi <- bayestestR::hdi(finite, ci = ci)

  tibble(
    estimate = s$median,
    median   = s$median,
    mean     = s$mean,
    sd       = s$sd,
    low      = hdi$CI_low,
    high     = hdi$CI_high,
    q5       = stats::quantile(finite, 0.05, names = FALSE),
    q95      = stats::quantile(finite, 0.95, names = FALSE),
    map      = suppressWarnings(
      as.numeric(bayestestR::map_estimate(finite, method = "KernSmooth")$MAP_Estimate)
    ),
    rhat     = s$rhat
  )
}


#' Posterior probability that a quantity is positive.
#'
#' Reported in the paper as the probability that model-related variation
#' exceeds scaffold-related variation (Section 5.1).
prob_positive <- function(x) {
  if (!inherits(x, "rvar")) return(as.numeric(x > 0))
  mean(as.numeric(posterior::draws_of(x)) > 0)
}
