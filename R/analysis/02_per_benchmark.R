# Per-benchmark variance decomposition (separate mixed models per benchmark).
# See docs/reliability_models.qmd and R/utils/model_formulas.R.

RUN_ESTIMATION <- FALSE
FIT_MODELS <- FALSE
RUN_DISCO <- FALSE
MAKE_PLOTS <- FALSE

source(here::here("R", "setup.R"))
source(here::here("R", "utils", "source_utils.R"))

df <- load_hal_data()
benches <- sort(unique(df$benchmark))
item_counts <- bench_item_counts(df)
fml_bench <- update(FORMULA_PER_BENCH_RE, score ~ .)

if (!RUN_ESTIMATION) {
  message("02_per_benchmark: skipping estimation (RUN_ESTIMATION=FALSE).")
  invisible(NULL)
} else {

# ---- Frequentist LME / GLME ----
lmmods <- list()
lmsums <- list()
glmmods <- list()
glmsums <- list()

for (b in benches) {
  d_b <- df |> filter(benchmark == b)
  lmmods[[b]] <- lmer(
    fml_bench,
    data = d_b,
    na.action = na.exclude,
    control = lmer_control_fast()
  )
  lmsums[[b]] <- varcorr_summary(lmmods[[b]])

  glmmods[[b]] <- glmer(
    fml_bench,
    data = d_b,
    family = binomial,
    nAGQ = 0,
    na.action = na.exclude,
    control = glmer_control()
  )
  glmsums[[b]] <- varcorr_summary(glmmods[[b]], add_glmer_residual = TRUE)
}

lmsums <- bind_rows(lmsums, .id = "benchmark")
glmsums <- bind_rows(glmsums, .id = "benchmark")

# ---- Bayesian (linear + Bernoulli), cached per benchmark ----
has_brm_cache <- any(file.exists(here("data", glue("brmfit_{benches}.rds"))))
if (!FIT_MODELS && !has_brm_cache) {
  message("Skipping Bayesian per-benchmark fits (no brmfit_{benchmark}.rds cache).")
  gblms <- NULL
} else {

blmmods <- list()
blmdraws <- list()
blmsums <- list()
gblmmods <- list()
gblmdraws <- list()
gblmsums <- list()

for (b in benches) {
  d_b <- df |> filter(benchmark == b)
  brm_file <- here("data", glue("brmfit_{b}.rds"))
  gbrm_file <- here("data", glue("gbrmfit_{b}.rds"))

  blmmods[[b]] <- fit_or_load_brm(
    d_b,
    fml_bench,
    file = brm_file,
    fit_models = FIT_MODELS,
    save_pars = save_pars(all = TRUE)
  )
  blmdraws[[b]] <- as_draws_rvars(blmmods[[b]])
  blmsums[[b]] <- summarize_brm_variance_draws(blmdraws[[b]]) |>
    mutate(benchmark = b)

  gblmmods[[b]] <- fit_or_load_brm(
    d_b,
    fml_bench,
    file = gbrm_file,
    fit_models = FIT_MODELS,
    family = bernoulli("logit"),
    save_pars = save_pars(all = TRUE),
    iter = 1000,
    thin = 2
  )
  gblmdraws[[b]] <- as_draws_rvars(gblmmods[[b]])
  gblmdraws[[b]]$sigma <- posterior::rvar(rnorm(1000, mean = pi^2 / 3, sd = 0))
  gblmsums[[b]] <- summarize_brm_variance_draws(gblmdraws[[b]]) |>
    mutate(benchmark = b)
}

gblms <- bind_rows(gblmsums)
}

# Per-benchmark Ep2 helpers (aliases for readability in plots)
get_bmod_ep2 <- get_bench_mod_ep2
get_bagnt_ep2 <- get_bench_agnt_ep2
get_bagmod_ep2 <- get_bench_agmod_ep2
get_rmod_ep2 <- get_bench_mod_ep2g
get_ragnt_ep2 <- get_bench_agnt_ep2g
get_ragmod_ep2 <- get_bench_agmod_ep2g

# ---- Nonparametric disco ----
disco_path <- here("data", "results", "disco_all_benches.csv")
if (RUN_DISCO) {
  disco_bench_res <- run_disco_by_benchmark(df, benches)
  write_csv(disco_bench_res, disco_path)
} else if (file.exists(disco_path)) {
  disco_bench_res <- read_csv(disco_path, show_col_types = FALSE)
} else {
  warning("No disco results at ", disco_path, "; set RUN_DISCO=TRUE or add cached CSV.")
  disco_bench_res <- NULL
}

# ---- Combined comparison table ----
if (!is.null(disco_bench_res)) {
  combbenchdf <- disco_bench_res |>
    mutate(
      facet = rename_disco_facet(facet),
      method = "nonp"
    ) |>
    select(method, benchmark, facet, pct) |>
    bind_rows(
      lmsums |>
        mutate(
          facet = finalize_facet(rename_lme_facet(grp)),
          method = "lme",
          pct = pct / 100
        ) |>
        select(facet, benchmark, method, pct)
    ) |>
    bind_rows(
      glmsums |>
        mutate(
          facet = finalize_facet(rename_lme_facet(grp)),
          method = "glme",
          pct = pct / 100
        ) |>
        select(facet, benchmark, method, pct)
    ) |>
    bind_rows(
      if (!is.null(gblms)) {
        gblms |>
          mutate(
            facet = finalize_facet(rename_brm_facet(var)),
            pct = mean,
            method = "bayes"
          ) |>
          select(facet, benchmark, method, pct)
      }
    )
}

if (!MAKE_PLOTS || is.null(gblms)) {
  message("02_per_benchmark.R: models summarized; plotting skipped.")
  invisible(list(lmsums = lmsums, glmsums = glmsums, gblms = gblms))
} else {

plot_ni <- 350

gbaydf <- map(benches, \(b) {
  nis <- min(item_counts$n[item_counts$benchmark == b], plot_ni)
  tibble(
    nitems = 1:nis,
    m_ep2 = get_rmod_ep2(gblmdraws[[b]], nitems = nitems) |> pull(median),
    m_ep2_low = get_rmod_ep2(gblmdraws[[b]], nitems = nitems) |> pull(HDI_low),
    m_ep2_hi = get_rmod_ep2(gblmdraws[[b]], nitems = nitems) |> pull(HDI_high),
    am_ep2 = get_ragmod_ep2(gblmdraws[[b]], nitems = nitems) |> pull(median),
    am_ep2_low = get_ragmod_ep2(gblmdraws[[b]], nitems = nitems) |> pull(HDI_low),
    am_ep2_hi = get_ragmod_ep2(gblmdraws[[b]], nitems = nitems) |> pull(HDI_high),
    benchmark = b
  )
}) |> bind_rows()

mvamp <- gbaydf |>
  inner_join(item_counts, by = "benchmark") |>
  filter(nitems <= n) |>
  pivot_longer(m_ep2:am_ep2_hi, names_to = c("facet", "stat", "val"), names_sep = "_", values_to = "bayes") |>
  mutate(val = if_else(is.na(val), "est", val)) |>
  pivot_wider(names_from = val, values_from = bayes) |>
  filter(facet %in% c("m", "am")) |>
  mutate(facet = if_else(facet == "m", "model", "agnt-mod")) |>
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
  paper_theme()

ggsave(here("figures", "bayes_CI_per_bench_mvamp.pdf"), mvamp, width = 14.3, height = 3.55, device = cairo_pdf)

p <- gbaydf |>
  inner_join(item_counts, by = "benchmark") |>
  filter(nitems == n) |>
  pivot_longer(m_ep2:m_ep2_hi, names_to = c("facet", "stat", "val"), names_sep = "_", values_to = "bayes") |>
  mutate(val = if_else(is.na(val), "est", val)) |>
  pivot_wider(names_from = val, values_from = bayes) |>
  filter(facet == "m") |>
  ggplot(aes(x = benchmark, y = est, color = benchmark, fill = benchmark)) +
  geom_pointrange(aes(ymin = low, ymax = hi), linewidth = 1.4, size = 1.4) +
  ylim(0, 1) +
  labs(y = expression("Estimated Reliability: E"*hat(rho)^2)) +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")

ggsave(here("figures", "bayes_CI_per_bench.pdf"), p, width = 7, height = 5, device = cairo_pdf)

gblms |>
  group_by(benchmark) |>
  mutate(across(mean:MAP_Estimate, \(x) x / sum(x)), facet = rename_brm_facet(var)) |>
  ggplot(aes(x = benchmark, y = mean, fill = facet)) +
  geom_col() +
  scale_fill_manual(values = component_cols, labels = component_labs) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Benchmark", y = "Proportion of Variance Explained", fill = "Component") +
  paper_theme(14) +
  theme(legend.position = "bottom")

invisible(list(lmsums = lmsums, glmsums = glmsums, gblms = gblms, combbenchdf = combbenchdf))
}
} # RUN_ESTIMATION
