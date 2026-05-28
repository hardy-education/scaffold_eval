# Per-benchmark decomposition with model and agent collapsed into one facet (`model`).
# Demonstrates ranking stability when treating scaffold+LLM as a single unit.

RUN_ESTIMATION <- FALSE
FIT_MODELS <- FALSE
RUN_DISCO <- FALSE
MAKE_PLOTS <- FALSE

source(here::here("R", "setup.R"))
source(here::here("R", "utils", "source_utils.R"))

df <- load_hal_data(collapse_model_agent = TRUE)
benches <- sort(unique(df$benchmark))
item_counts <- bench_item_counts(df)
fml_ma <- update(FORMULA_MODEL_AGENT_RE, score ~ .)

if (!RUN_ESTIMATION) {
  message("03_per_benchmark_model_agent: skipping estimation (RUN_ESTIMATION=FALSE).")
  invisible(NULL)
} else {

# ---- Frequentist ----
lmsums2 <- list()
glmsums2 <- list()
for (b in benches) {
  d_b <- df |> filter(benchmark == b)
  lmsums2[[b]] <- varcorr_summary(lmer(
    fml_ma, data = d_b, na.action = na.exclude, control = lmer_control_fast()
  ))
  glmsums2[[b]] <- varcorr_summary(glmer(
    fml_ma, data = d_b, family = binomial, nAGQ = 0,
    na.action = na.exclude, control = glmer_control()
  ), add_glmer_residual = TRUE)
}
lmsums2 <- bind_rows(lmsums2, .id = "benchmark")
glmsums2 <- bind_rows(glmsums2, .id = "benchmark")

# ---- Bayesian ----
has_ma_cache <- any(file.exists(here("data", glue("brmfit2_{benches}.rds"))))
if (!FIT_MODELS && !has_ma_cache) {
  message("Skipping Bayesian model-agent fits (no brmfit2_{benchmark}.rds cache).")
  gblms2 <- blmdraws2 <- NULL
} else {

blmdraws2 <- list()
gblmsums2 <- list()
for (b in benches) {
  d_b <- df |> filter(benchmark == b)
  blm <- fit_or_load_brm(
    d_b, fml_ma,
    file = here("data", glue("brmfit2_{b}.rds")),
    fit_models = FIT_MODELS,
    save_pars = save_pars(all = TRUE)
  )
  blmdraws2[[b]] <- as_draws_rvars(blm)

  gblm <- fit_or_load_brm(
    d_b, fml_ma,
    file = here("data", glue("gbrmfit2_{b}.rds")),
    fit_models = FIT_MODELS,
    family = bernoulli("logit"),
    save_pars = save_pars(all = TRUE),
    iter = 1000,
    thin = 2
  )
  gdraw <- as_draws_rvars(gblm)
  gdraw$sigma <- posterior::rvar(rnorm(1000, mean = pi^2 / 3, sd = 0))
  gblmsums2[[b]] <- summarize_brm_variance_draws(gdraw) |> mutate(benchmark = b)
}
gblms2 <- bind_rows(gblmsums2)
}

get_bmod2_ep2 <- get_collapsed_mod_ep2

# ---- Disco ----
disco_path <- here("data", "results", "disco_all_benches2.csv")
if (RUN_DISCO) {
  disco_bench_res2 <- map(benches, \(b) {
    disco_decompose(df, benchmark = b, collapse_model_agent = TRUE) |>
      mutate(benchmark = b)
  }) |> bind_rows()
  write_csv(disco_bench_res2, disco_path)
} else if (file.exists(disco_path)) {
  disco_bench_res2 <- read_csv(disco_path, show_col_types = FALSE)
} else {
  warning("No disco results at ", disco_path)
  disco_bench_res2 <- NULL
}

if (!MAKE_PLOTS || is.null(gblms2)) {
  invisible(list(lmsums2 = lmsums2, gblms2 = gblms2))
} else {

plot_ni <- 350
gbaydf2 <- map(benches, \(b) {
  nis <- min(item_counts$n[item_counts$benchmark == b], plot_ni)
  tibble(
    nitems = 1:nis,
    m_ep2 = get_bmod2_ep2(blmdraws2[[b]], nitems = nitems) |> pull(median),
    m_ep2_low = get_bmod2_ep2(blmdraws2[[b]], nitems = nitems) |> pull(HDI_low),
    m_ep2_hi = get_bmod2_ep2(blmdraws2[[b]], nitems = nitems) |> pull(HDI_high),
    benchmark = b
  )
}) |> bind_rows()

p_ci <- gbaydf2 |>
  inner_join(item_counts, by = "benchmark") |>
  filter(nitems <= n) |>
  pivot_longer(m_ep2:m_ep2_hi, names_to = c("facet", "stat", "val"), names_sep = "_", values_to = "bayes") |>
  mutate(val = if_else(is.na(val), "est", val)) |>
  pivot_wider(names_from = val, values_from = bayes) |>
  filter(facet == "m") |>
  ggplot(aes(x = nitems, y = est, color = facet, fill = facet)) +
  geom_ribbon(aes(ymin = low, ymax = hi), alpha = 0.2, linewidth = 0) +
  geom_line() +
  scale_x_log10() +
  facet_wrap(~benchmark, scales = "free_x", nrow = 1) +
  ylim(0, 1) +
  labs(
    x = "Number of Agentic Tasks",
    y = expression("Estimated Reliability: E"*hat(rho)^2)
  ) +
  paper_theme() +
  theme(legend.position = "none")

ggsave(here("figures", "bayes_CI_mod-age_per_bench.pdf"), p_ci, width = 7, height = 5, device = cairo_pdf)

p_bar <- gblms2 |>
  group_by(benchmark) |>
  mutate(across(mean:MAP_Estimate, \(x) x / sum(x)), facet = rename_brm_facet(var)) |>
  ggplot(aes(x = benchmark, y = median, fill = facet)) +
  geom_col() +
  scale_fill_manual(values = component_cols[c("item", "model", "sigma")], labels = component_labs) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Benchmark", y = "Limit of Proportion of Variance Explained", fill = "Component") +
  paper_theme() +
  theme(legend.position = "bottom")

ggsave(here("figures", "mod_agent_reliab_bar.pdf"), p_bar, width = 15.1, height = 4.8, device = cairo_pdf)

invisible(list(lmsums2 = lmsums2, gblms2 = gblms2))
}
} # RUN_ESTIMATION
