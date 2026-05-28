# Generalizability (Ep^2) from Bayesian posterior SD draws.
# Numerators isolate a facet; denominators sum relevant variance components,
# scaling item-level terms by nitems (and benchmark-level terms by nbench when pooled).

plogis_ep2 <- function(x) stats::plogis(x)

# ---- Pooled (multi-benchmark) formulas ----

get_mod_ep2 <- function(drawlist, nitems = 100, nbench = 1) {
  num <- drawlist[["sd_model_name__Intercept"]]^2
  den <- num +
    drawlist[["sd_benchmark:model_name__Intercept"]]^2 +
    drawlist[["sd_benchmark:task_id:model_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_model_name:agent_name__Intercept"]]^2 +
    drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]^2 / nbench +
    drawlist[["sigma"]]^2 / (nitems * nbench)
  make_draws_stats(num / den)
}

get_bmod_ep2 <- function(drawlist, nitems = 100, nbench = 1) {
  num <- drawlist[["sd_benchmark:model_name__Intercept"]]^2
  den <- num +
    drawlist[["sd_benchmark:task_id:model_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]^2 / nbench +
    drawlist[["sigma"]]^2 / (nitems * nbench)
  make_draws_stats(num / den)
}

get_agnt_ep2 <- function(drawlist, nitems = 100, nbench = 1) {
  num <- drawlist[["sd_agent_name__Intercept"]]^2
  den <- num +
    drawlist[["sd_benchmark:agent_name__Intercept"]]^2 +
    drawlist[["sd_benchmark:task_id:agent_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_model_name:agent_name__Intercept"]]^2 +
    drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]^2 / nbench +
    drawlist[["sigma"]]^2 / (nitems * nbench)
  make_draws_stats(num / den)
}

get_bagnt_ep2 <- function(drawlist, nitems = 100, nbench = 1) {
  num <- drawlist[["sd_benchmark:agent_name__Intercept"]]^2
  den <- num +
    drawlist[["sd_benchmark:task_id:agent_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]^2 / nbench +
    drawlist[["sigma"]]^2 / (nitems * nbench)
  make_draws_stats(num / den)
}

get_agmod_ep2 <- function(drawlist, nitems = 100, nbench = 1) {
  num <- drawlist[["sd_model_name:agent_name__Intercept"]]^2
  den <- num +
    drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]^2 / nbench +
    drawlist[["sigma"]]^2 / (nitems * nbench)
  make_draws_stats(num / den)
}

get_bagmod_ep2 <- function(drawlist, nitems = 100, nbench = 1) {
  num <- drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]^2
  den <- num / nbench +
    drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2 / (nitems * nbench) +
    drawlist[["sigma"]]^2 / (nitems * nbench)
  make_draws_stats(num / den)
}

# Logit-scale (generalized / Bernoulli) variants
get_mod_ep2g <- function(drawlist, nitems = 100, nbench = 1) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_model_name__Intercept"]])
  den <- num + p(drawlist[["sd_benchmark:model_name__Intercept"]]) +
    p(drawlist[["sd_benchmark:task_id:model_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_model_name:agent_name__Intercept"]]) +
    p(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]) / nbench +
    p(drawlist[["sigma"]]) / (nitems * nbench)
  make_draws_stats(num / den)
}

get_bmod_ep2g <- function(drawlist, nitems = 100, nbench = 1) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_benchmark:model_name__Intercept"]])
  den <- num + p(drawlist[["sd_benchmark:task_id:model_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]) / nbench +
    p(drawlist[["sigma"]]) / (nitems * nbench)
  make_draws_stats(num / den)
}

get_agnt_ep2g <- function(drawlist, nitems = 100, nbench = 1) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_agent_name__Intercept"]])
  den <- num + p(drawlist[["sd_benchmark:agent_name__Intercept"]]) +
    p(drawlist[["sd_benchmark:task_id:agent_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_model_name:agent_name__Intercept"]]) +
    p(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]) / nbench +
    p(drawlist[["sigma"]]) / (nitems * nbench)
  make_draws_stats(num / den)
}

get_agmod_ep2g <- function(drawlist, nitems = 100, nbench = 1) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_model_name:agent_name__Intercept"]])
  den <- num +
    p(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]) / nbench +
    p(drawlist[["sigma"]]) / (nitems * nbench)
  make_draws_stats(num / den)
}

get_bagnt_ep2g <- function(drawlist, nitems = 100, nbench = 1) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_benchmark:agent_name__Intercept"]])
  den <- num + p(drawlist[["sd_benchmark:task_id:agent_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]) / nbench +
    p(drawlist[["sigma"]]) / (nitems * nbench)
  make_draws_stats(num / den)
}

get_bagmod_ep2g <- function(drawlist, nitems = 100, nbench = 1) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_benchmark:model_name:agent_name__Intercept"]])
  den <- p(drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]) / nbench +
    p(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]) / (nitems * nbench) +
    p(drawlist[["sigma"]]) / (nitems * nbench)
  make_draws_stats(num / den)
}

# ---- Per-benchmark formulas (no benchmark random effects in fit) ----

get_bench_mod_ep2 <- function(drawlist, nitems = 100) {
  num <- drawlist[["sd_model_name__Intercept"]]^2
  den <- num +
    drawlist[["sd_task_id:model_name__Intercept"]]^2 / nitems +
    drawlist[["sd_model_name:agent_name__Intercept"]]^2 +
    drawlist[["sigma"]]^2 / nitems
  make_draws_stats(num / den)
}

get_bench_agnt_ep2 <- function(drawlist, nitems = 100) {
  num <- drawlist[["sd_agent_name__Intercept"]]^2
  den <- num +
    drawlist[["sd_task_id:agent_name__Intercept"]]^2 / nitems +
    drawlist[["sd_model_name:agent_name__Intercept"]]^2 +
    drawlist[["sigma"]]^2 / nitems
  make_draws_stats(num / den)
}

get_bench_agmod_ep2 <- function(drawlist, nitems = 100) {
  num <- drawlist[["sd_model_name:agent_name__Intercept"]]^2
  den <- num + drawlist[["sigma"]]^2 / nitems
  make_draws_stats(num / den)
}

get_bench_mod_ep2g <- function(drawlist, nitems = 100) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_model_name__Intercept"]])
  den <- num + p(drawlist[["sd_task_id:model_name__Intercept"]]) / nitems +
    p(drawlist[["sd_model_name:agent_name__Intercept"]]) +
    p(drawlist[["sigma"]]) / nitems
  make_draws_stats(num / den)
}

get_bench_agnt_ep2g <- function(drawlist, nitems = 100) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_agent_name__Intercept"]])
  den <- num + p(drawlist[["sd_task_id:agent_name__Intercept"]]) / nitems +
    p(drawlist[["sd_model_name:agent_name__Intercept"]]) +
    p(drawlist[["sigma"]]) / nitems
  make_draws_stats(num / den)
}

get_bench_agmod_ep2g <- function(drawlist, nitems = 100) {
  p <- function(x) plogis_ep2(x^2)
  num <- p(drawlist[["sd_model_name:agent_name__Intercept"]])
  den <- num + p(drawlist[["sigma"]]) / nitems
  make_draws_stats(num / den)
}

# Collapsed model–agent facet (single `model` random effect)
get_collapsed_mod_ep2 <- function(drawlist, nitems = 100) {
  num <- drawlist[["sd_model__Intercept"]]^2
  den <- num + drawlist[["sigma"]]^2 / nitems
  make_draws_stats(num / den)
}
