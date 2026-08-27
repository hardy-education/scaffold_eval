# scripts/06_method_comparison.R -------------------------------------------
#
# Does the decomposition depend on the estimator?
#
# Compares four ways of splitting the same variance (paper Appendix C.2-C.5,
# G.3):
#
#   LME    Gaussian identity link, REML. The estimator most G-theory
#          applications use, and the one this paper argues against for pass/fail
#          outcomes: its residual absorbs the Bernoulli variance floor p(1-p),
#          so named facets are systematically understated.
#   GLME   Bernoulli logit link, Laplace approximation. Respects the binary
#          outcome but reports a point estimate near the posterior mode, which
#          misses the right skew in scaffold variance.
#   Bayes  Bernoulli logit, the paper's primary estimator. Partial pooling for
#          sparse cells, and uncertainty propagated into every derived ratio.
#   DISCO  Nonparametric distance components. No residual term, so unexplained
#          variation is pushed into the named interactions; best read as an
#          upper bound on interaction effects.
#
# The frequentist fits run here directly and are the slow part: the pooled
# lme4 fits over 29,923 rows with ten crossed random effects take roughly
# 30-60 minutes between them. The Bayesian columns come from the cached fits;
# DISCO is cached as CSV because it takes hours.

source(here::here("R", "setup.R"))

df <- load_hal_data()
RUN_DISCO <- env_flag("RUN_DISCO")


# ---- Leaderboard level (Table 10) ------------------------------------------

message("=== leaderboard decomposition across estimators ===")

message("fitting LME (Gaussian identity link)...")
lme_fit  <- fit_lme4(df, FORMULA_LEADERBOARD, family = "gaussian")

message("fitting GLME (Bernoulli logit link)...")
glme_fit <- fit_lme4(df, FORMULA_LEADERBOARD, family = "binomial", nagq = 0L)

vc_lme   <- vcomp_from_lme4(lme_fit,  level = "leaderboard")
vc_glme  <- vcomp_from_lme4(glme_fit, level = "leaderboard")
vc_bayes <- vcomp_from_brms(fit_leaderboard(df), level = "leaderboard")

shares_by_method <- bind_rows(
  variance_shares(vc_lme)   |> mutate(method = "LME"),
  variance_shares(vc_glme)  |> mutate(method = "GLME"),
  variance_shares(vc_bayes) |> mutate(method = "Bayes")
) |>
  select(method, component, label, estimate)

if (RUN_DISCO) {
  message("running DISCO on the pooled data (~17 h)...")
  disco_pooled <- disco_or_cache(df, "disco_leaderboard", level = "leaderboard")
  shares_by_method <- bind_rows(
    shares_by_method,
    disco_pooled |> transmute(method = "DISCO", component = component,
                              label = facet, estimate = pct)
  )
}

leaderboard_comparison <- shares_by_method |>
  tidyr::pivot_wider(names_from = method, values_from = estimate) |>
  arrange(desc(coalesce(Bayes, 0)))

print(leaderboard_comparison, n = Inf)
save_table(leaderboard_comparison, "table10_leaderboard_method_comparison",
           caption = "Proportion of variation explained, full leaderboard, by estimator.")

p_methods <- shares_by_method |>
  mutate(component = factor(component, levels = names(COMPONENT_COLOURS))) |>
  ggplot(aes(x = method, y = estimate, fill = component, alpha = component)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = COMPONENT_COLOURS, na.value = "grey70") +
  scale_alpha_manual(values = COMPONENT_ALPHAS, guide = "none", na.value = 1) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = "Proportion of variance explained", fill = "Component") +
  paper_theme()

save_figure(p_methods, "method_comparison_leaderboard", width = 9, height = 4)


# ---- Where the estimators disagree, and why --------------------------------

resid_share <- shares_by_method |>
  filter(component == "e") |>
  select(method, residual_share = estimate)

print(resid_share)
message(
  "\nThe Gaussian fit assigns a markedly larger share to the residual: an ",
  "identity-link\nmodel of a binary outcome has to absorb the Bernoulli ",
  "variance floor p(1-p) there,\nwhereas the logit-link models separate it ",
  "out through the distributional assumption."
)


# ---- Benchmark level (Table 9) ---------------------------------------------

message("\n=== per-benchmark decomposition across estimators ===")

bench_levels <- levels(droplevels(df$benchmark))
bayes_bench  <- purrr::map(fit_benchmarks(df), vcomp_from_brms, level = "benchmark")

bench_comparison <- purrr::map(bench_levels, function(b) {
  d <- filter(df, benchmark == b)
  bind_rows(
    variance_shares(vcomp_from_lme4(
      fit_lme4(d, FORMULA_BENCHMARK, family = "gaussian"), "benchmark")) |>
      mutate(method = "LME"),
    variance_shares(vcomp_from_lme4(
      fit_lme4(d, FORMULA_BENCHMARK, family = "binomial"), "benchmark")) |>
      mutate(method = "GLME"),
    variance_shares(bayes_bench[[b]]) |> mutate(method = "Bayes")
  ) |>
    mutate(benchmark = b, .before = 1)
}) |> bind_rows()

if (RUN_DISCO) {
  disco_bench <- disco_or_cache(df, "disco_benchmarks", level = "benchmark")
  bench_comparison <- bind_rows(
    bench_comparison,
    disco_bench |> transmute(benchmark, component = component,
                             label = facet, estimate = pct, method = "DISCO")
  )
}

bench_table <- bench_comparison |>
  select(benchmark, component, label, method, estimate) |>
  tidyr::pivot_wider(names_from = method, values_from = estimate) |>
  arrange(benchmark, component)

print(bench_table, n = 40)
save_table(bench_table, "table9_benchmark_method_comparison",
           caption = "Proportion of variation explained per benchmark, by estimator.")


# ---- Do the estimators change the design conclusion? -----------------------
#
# The variance shares differ; the question that matters is whether the
# reliability projection does. This recomputes the headline ceiling under each
# estimator.

ceiling_by_method <- tibble(
  method = c("LME", "GLME", "Bayes"),
  ceiling_1_bench = c(
    reliability_ceiling(vc_lme),
    reliability_ceiling(vc_glme),
    summarize_quantity(reliability_ceiling(vc_bayes))$median
  ),
  ceiling_9_bench = c(
    reliability_ceiling(vc_lme,   n_benchmarks = 9),
    reliability_ceiling(vc_glme,  n_benchmarks = 9),
    summarize_quantity(reliability_ceiling(vc_bayes, n_benchmarks = 9))$median
  )
)

print(ceiling_by_method)
save_table(ceiling_by_method, "table_ceiling_by_method",
           caption = paste("Task-only model-ranking ceiling under each",
                           "variance-component estimator."))

message("\n06_method_comparison.R complete.")
