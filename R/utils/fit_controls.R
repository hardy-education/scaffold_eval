lmer_control <- function() {
  lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
}

lmer_control_fast <- function() {
  lme4::lmerControl(
    optimizer = "bobyqa",
    calc.derivs = FALSE,
    optCtrl = list(maxfun = 2e5)
  )
}

glmer_control <- function(nagq = 0L) {
  lme4::glmerControl(
    optimizer = "bobyqa",
    calc.derivs = isTRUE(nagq == 0L),
    optCtrl = list(maxfun = 2e5)
  )
}

brm_defaults <- function(
    iter = 2000,
    thin = 5,
    chains = 4,
    cores = 4) {
  list(
    chains = chains,
    cores = cores,
    thin = thin,
    iter = iter,
    backend = "cmdstanr",
    control = list(adapt_delta = 0.95),
    threads = brms::threading(2),
    stan_model_args = list(stanc_options = list("O1"))
  )
}

brm_defaults_fast <- function() {
  modifyList(brm_defaults(iter = 1000, thin = 2), list())
}
