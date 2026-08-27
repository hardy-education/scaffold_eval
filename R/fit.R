# R/fit.R ------------------------------------------------------------------
#
# Fitting the variance-decomposition models, and caching them.
#
# The pooled Bayesian fit takes hours. Everything here therefore follows a
# fit-or-load discipline: a fit is written to outputs/fits/<name>.rds, and
# subsequent runs read it back unless REFIT is TRUE. Analysis scripts never
# refit implicitly, so re-running a figure script is cheap and cannot silently
# produce numbers from a different fit than the one it reports.

#' Path of the cache file for a named fit.
fit_path <- function(name) file.path(PATHS$fits, paste0(name, ".rds"))


#' Fit a Bayesian variance decomposition, or load the cached fit.
#'
#' Priors are the brms defaults described in Appendix C.1.3: Student-t(3, 0,
#' 2.5) truncated at zero on every group-level standard deviation, and the
#' corresponding weakly informative prior on the intercept. These regularize
#' components supported by few levels -- notably scaffold and
#' benchmark-scaffold effects -- without forcing them toward zero.
#'
#' @param data Rollout data.
#' @param formula One of the structures in R/formulas.R.
#' @param name Cache name, without extension.
#' @param family `brms::bernoulli("logit")` for the latent-scale decomposition
#'   used throughout the paper; `gaussian()` for the linear sensitivity
#'   analysis of Appendix C.2.
#' @param sampler A list from `SAMPLER`.
#' @param refit Force re-estimation.
fit_bayes <- function(data,
                      formula,
                      name,
                      family = brms::bernoulli("logit"),
                      sampler = SAMPLER$leaderboard,
                      refit = REFIT) {
  path <- fit_path(name)

  if (file.exists(path) && !refit) {
    message("loading cached fit: ", basename(path))
    return(readRDS(path))
  }
  if (!refit) {
    stop(
      "No cached fit at ", path, ".\n",
      "Set REFIT <- TRUE (or REFIT=TRUE in the environment) to estimate it. ",
      "See the timing table in README.md before you do.",
      call. = FALSE
    )
  }

  message("fitting: ", name, " (", nrow(data), " rollouts)")
  fit <- brms::brm(
    formula = formula,
    data    = data,
    family  = family,
    chains  = sampler$chains,
    cores   = sampler$cores,
    iter    = sampler$iter,
    warmup  = sampler$warmup,
    thin    = sampler$thin,
    seed    = SAMPLER$seed,
    backend = SAMPLER$backend,
    control = list(adapt_delta = SAMPLER$adapt_delta),
    threads = brms::threading(SAMPLER$threads),
    save_pars = brms::save_pars(all = TRUE),
    stan_model_args = list(stanc_options = list("O1"))
  )

  saveRDS(fit, path)
  fit
}


#' Fit the pooled leaderboard decomposition (Eq. 6).
fit_leaderboard <- function(df, name = "leaderboard_bernoulli", ...) {
  fit_bayes(df, FORMULA_LEADERBOARD, name,
            sampler = SAMPLER$leaderboard, ...)
}


#' Fit the benchmark-level decomposition (Eq. 2) within each benchmark.
#'
#' @return A named list of `brmsfit` objects, one per benchmark.
fit_benchmarks <- function(df, prefix = "benchmark_bernoulli", ...) {
  benches <- levels(droplevels(df$benchmark))
  purrr::set_names(benches) |>
    purrr::map(\(b) {
      fit_bayes(filter(df, benchmark == b), FORMULA_BENCHMARK,
                paste0(prefix, "_", b), sampler = SAMPLER$benchmark, ...)
    })
}


#' Frequentist variance decomposition via lme4.
#'
#' Used where a full Bayesian refit is impractical: the leave-one-model-out and
#' leave-one-benchmark-scaffold-out ablations require tens to hundreds of
#' refits (Appendix D.5-D.6), and the LME/GLME contrast of Appendix C.3-C.4
#' needs both link functions.
#'
#' @param family "binomial" for the logit-link decomposition, "gaussian" for
#'   the identity-link baseline that is standard in G-theory applications.
#' @param nagq Adaptive Gauss-Hermite points for `glmer`. `0` is far faster and
#'   is what the ablation loops use; the default Laplace approximation (`1`)
#'   gives more accurate variance components for single reported fits.
fit_lme4 <- function(data, formula, family = c("binomial", "gaussian"), nagq = 0L) {
  family <- match.arg(family)

  if (family == "gaussian") {
    lme4::lmer(
      formula, data = data, na.action = stats::na.exclude,
      control = lme4::lmerControl(
        optimizer = "bobyqa", calc.derivs = FALSE,
        optCtrl = list(maxfun = 2e5)
      )
    )
  } else {
    lme4::glmer(
      formula, data = data, family = stats::binomial(),
      na.action = stats::na.exclude, nAGQ = nagq,
      control = lme4::glmerControl(
        optimizer = "bobyqa", calc.derivs = FALSE,
        optCtrl = list(maxfun = 2e5)
      )
    )
  }
}


# ---- Convergence diagnostics ----------------------------------------------

#' Convergence summary for the variance parameters of a Bayesian fit.
#'
#' The paper reports Rhat and effective sample sizes for every variance
#' component; this reproduces that table for any fit.
convergence_summary <- function(fit) {
  posterior::as_draws_rvars(fit) |>
    (\(d) d[stringr::str_starts(names(d), "sd_|sigma")])() |>
    purrr::imap(\(x, k) {
      posterior::summarise_draws(x, "rhat", "ess_bulk", "ess_tail") |>
        mutate(parameter = k, .before = 1) |>
        select(-any_of("variable"))
    }) |>
    bind_rows()
}


#' Approximate leave-one-out cross-validation (paper Appendix D.2).
#'
#' Reports elpd and the Pareto-k diagnostic. Moment matching is expensive --
#' roughly a day and a half single-threaded on the pooled fit -- so it is off
#' by default; without it, observations flagged with k > 0.7 are the ones to
#' inspect, and the paper's finding is that they are almost all singleton
#' interaction levels rather than evidence of model failure.
psis_loo <- function(fit, name, moment_match = FALSE, refit = REFIT) {
  path <- fit_path(paste0(name, "_loo"))
  if (file.exists(path) && !refit) return(readRDS(path))
  if (!refit) stop("No cached LOO at ", path, "; set REFIT <- TRUE.", call. = FALSE)

  out <- brms::loo(fit, moment_match = moment_match, save_psis = TRUE)
  saveRDS(out, path)
  out
}


#' Observations whose predictive distribution is sensitive to their own removal.
#'
#' @return The rows of `df` with Pareto k above `threshold`, annotated with how
#'   many rollouts share their (benchmark, task, model) group. Singleton groups
#'   are the expected source of high k in this design.
influential_observations <- function(loo_obj, df, threshold = 0.7) {
  k <- loo_obj$psis_object$diagnostics$pareto_k
  group_size <- df |>
    add_count(benchmark, task_id, model_name, name = "group_n") |>
    pull(group_n)

  df[k > threshold, ] |>
    mutate(pareto_k = k[k > threshold], group_n = group_size[k > threshold]) |>
    arrange(desc(pareto_k))
}
