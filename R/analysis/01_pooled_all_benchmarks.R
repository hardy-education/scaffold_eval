# Pooled variance decomposition across all benchmarks (full crossed random-effects hierarchy).
# Produces paper figures: variance bars, reliability curves, BLUP ranking diagnostics.

# Set RUN_ESTIMATION <- TRUE only when you intend to refit (pooled LME/GLME ~30+ min).
RUN_ESTIMATION <- FALSE
FIT_MODELS <- FALSE
RUN_DISCO <- FALSE
MAKE_PLOTS <- FALSE

source(here::here("R", "setup.R"))
source(here::here("R", "utils", "source_utils.R"))

df <- load_hal_data()
benches <- sort(unique(df$benchmark))
fml_pooled <- update(FORMULA_POOLED_RE, score ~ .)

if (!RUN_ESTIMATION) {
  message(
    "01_pooled: skipping estimation (RUN_ESTIMATION=FALSE). ",
    "Load cached data/*.rds and set MAKE_PLOTS=TRUE to regenerate figures only."
  )
  invisible(NULL)
} else {

# ---- Frequentist LME ----
tic("lme pooled")
m <- lmer(fml_pooled, data = df, na.action = na.exclude, control = lmer_control())
mvar <- varcorr_summary(m)
toc()

# ---- GLME (accurate vs fast Laplace) ----
tic("glme pooled")
gm <- glmer(
  fml_pooled, data = df, family = binomial, na.action = na.exclude,
  control = glmer_control()
)
gmvar <- varcorr_summary(gm, add_glmer_residual = TRUE)
toc()

tic("glme pooled fast")
gm2 <- glmer(
  fml_pooled, data = df, family = binomial, nAGQ = 0,
  na.action = na.exclude, control = glmer_control(nagq = 0L)
)
gmvar2 <- varcorr_summary(gm2, add_glmer_residual = TRUE)
toc()

# ---- Bayesian linear & Bernoulli ----
has_pooled_brm <- file.exists(here("data", "brmfit_.rds")) ||
  file.exists(here("data", "gbrmfit_.rds"))

if (!FIT_MODELS && !has_pooled_brm) {
  message("Skipping pooled Bayesian fits (no brmfit_.rds / gbrmfit_.rds). Frequentist results still available.")
  bms <- gbms <- NULL
  gbmd <- NULL
} else {

bm <- fit_or_load_brm(
  df, fml_pooled,
  file = here("data", "brmfit_.rds"),
  fit_models = FIT_MODELS,
  save_pars = save_pars(all = TRUE)
)
bmd <- as_draws_rvars(bm)
bms <- summarize_brm_variance_draws(bmd)

gbm <- fit_or_load_brm(
  df, fml_pooled,
  file = here("data", "gbrmfit_.rds"),
  fit_models = FIT_MODELS,
  family = bernoulli(link = "logit"),
  save_pars = save_pars(all = TRUE)
)
gbmd <- as_draws_rvars(gbm)
gbmd$sigma <- posterior::rvar(rnorm(800, mean = pi^2 / 3, sd = 0))
gbms <- summarize_brm_variance_draws(gbmd)
}

# ---- Nonparametric disco (optional; slow) ----
disco_path <- here("data", "results", "disco_all_model.csv")
if (RUN_DISCO) {
  tic("disco pooled")
  disco_results <- disco_decompose(df, benchmark = NULL, alpha = 1)
  write_csv(disco_results, disco_path)
  toc()
} else if (file.exists(disco_path)) {
  disco_results <- read_csv(disco_path, show_col_types = FALSE)
} else {
  warning("No disco cache at ", disco_path)
  disco_results <- NULL
}

# ---- Combined method comparison ----
combdf <- NULL
if (!is.null(disco_results)) {
  combdf <- disco_results |>
    mutate(facet = rename_disco_facet(facet), method = "nonp") |>
    select(method, facet, pct) |>
    bind_rows(
      mvar |> mutate(facet = finalize_facet(rename_lme_facet(grp)), method = "lme", pct = pct / 100) |>
        select(facet, method, pct)
    ) |>
    bind_rows(
      gmvar |> mutate(facet = finalize_facet(rename_lme_facet(grp)), method = "glme", pct = pct / 100) |>
        select(facet, method, pct)
    ) |>
    bind_rows(
      gmvar2 |> mutate(facet = finalize_facet(rename_lme_facet(grp)), method = "glme2", pct = pct / 100) |>
        select(facet, method, pct)
    ) |>
    bind_rows(if (!is.null(bms)) {
      bms |> mutate(facet = finalize_facet(rename_brm_facet(var)), pct = mean, method = "bayes") |>
        select(facet, method, pct)
    }) |>
    bind_rows(if (!is.null(gbms)) {
      gbms |> mutate(facet = finalize_facet(rename_brm_facet(var)), pct = mean, method = "gbayes") |>
        select(facet, method, pct)
    })
}

if (!MAKE_PLOTS || is.null(gbms)) {
  message("01_pooled: estimation complete; plots skipped.")
  invisible(list(mvar = mvar, gbms = gbms, combdf = combdf))
} else {

# ---- Paper: Bernoulli variance bar ----
bp <- gbms |>
  mutate(
    across(mean:MAP_Estimate, \(x) x / sum(x)),
    facet = finalize_facet(rename_brm_facet(var)),
    method = "Bern. (Bayes)",
    facet = factor(facet, levels = names(component_labs), ordered = TRUE)
  ) |>
  inner_join(tibble(facet = names(component_labs), labels = unname(component_labs)), by = "facet") |>
  ggplot(aes(x = method, y = median, fill = facet, alpha = facet)) +
  geom_col() +
  geom_label_repel(aes(label = labels), position = position_stack(vjust = 0.5), direction = "y", size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = component_cols, labels = component_labs) +
  scale_alpha_manual(values = component_alphas, labels = component_labs) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Proportion of Variance Explained", fill = "Component", alpha = "Component") +
  paper_theme(14) +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(here("figures", "bayes_variance_decomp_full.pdf"), bp, width = 15.6, height = 1.66, device = cairo_pdf)

# ---- Paper: reliability vs number of items (generalized Bayes) ----
plot_ni <- 50
n_bench <- 7
gbayescidf <- tibble(
  nitems = 1:plot_ni,
  m_ep2 = get_bmod_ep2g(gbmd, nitems = nitems, nbench = n_bench) |> pull(mean),
  m_ep2_low = get_bmod_ep2g(gbmd, nitems = nitems, nbench = n_bench) |> pull(q5),
  m_ep2_hi = get_bmod_ep2g(gbmd, nitems = nitems, nbench = n_bench) |> pull(q95)
)

gbayescidf1 <- tibble(
  nitems = 1:350,
  m_ep2 = get_bmod_ep2g(gbmd, nitems = nitems, nbench = 1) |> pull(mean),
  m_ep2_low = get_bmod_ep2g(gbmd, nitems = nitems, nbench = 1) |> pull(q5),
  m_ep2_hi = get_bmod_ep2g(gbmd, nitems = nitems, nbench = 1) |> pull(q95)
)

gbayescidf2 <- tibble(
  nitems = 1:175,
  m_ep2 = get_bmod_ep2g(gbmd, nitems = nitems, nbench = 2) |> pull(mean),
  m_ep2_low = get_bmod_ep2g(gbmd, nitems = nitems, nbench = 2) |> pull(q5),
  m_ep2_hi = get_bmod_ep2g(gbmd, nitems = nitems, nbench = 2) |> pull(q95)
)

p_rel <- gbayescidf |>
  pivot_longer(m_ep2:m_ep2_hi, names_to = c(".value", "band"), names_pattern = "m_ep2_(.*)") |>
  mutate(nitems = nitems * 7, Benches_Sampled = "7 benches") |>
  bind_rows(
    gbayescidf1 |>
      pivot_longer(m_ep2:m_ep2_hi, names_to = c(".value", "band"), names_pattern = "m_ep2_(.*)") |>
      filter(nitems > 6) |>
      mutate(Benches_Sampled = "1 bench")
  ) |>
  bind_rows(
    gbayescidf2 |>
      pivot_longer(m_ep2:m_ep2_hi, names_to = c(".value", "band"), names_pattern = "m_ep2_(.*)") |>
      mutate(nitems = nitems * 2) |>
      filter(nitems > 6) |>
      mutate(Benches_Sampled = "2 benches")
  ) |>
  pivot_wider(names_from = band, values_from = value) |>
  ggplot(aes(x = nitems, y = mean, color = Benches_Sampled, fill = Benches_Sampled)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2, linewidth = 0) +
  geom_line(linewidth = 1.5) +
  scale_x_log10() +
  ylim(0, 1) +
  labs(
    x = "Total Number of Agentic Tasks",
    y = expression("Estimated Reliability: E"*hat(rho)^2)
  ) +
  paper_theme() +
  theme(legend.position = "bottom")

ggsave(here("figures", "bayes_reliability_by_benchmark_count.pdf"), p_rel, width = 8, height = 5, device = cairo_pdf)

# ---- BLUP ranking diagnostics (generalized model) ----
bench_pat <- paste0("(", paste(benches, collapse = "|"), ")")

blupdf <- df |>
  group_by(benchmark, model_name) |>
  summarize(score = mean(score, na.rm = TRUE), .groups = "drop") |>
  group_by(benchmark) |>
  mutate(oldrank = percent_rank(score)) |>
  ungroup() |>
  inner_join(
    ranef(gbm, groups = "benchmark:model_name")[[1]] |>
      as_tibble(rownames = "var") |>
      mutate(
        benchmark = str_match(var, bench_pat)[, 2],
        model_name = str_split_i(var, str_c(benchmark, "_"), 2)
      ) |>
      clean_names(),
    by = c("benchmark", "model_name")
  )

blupdf_ranked <- blupdf |>
  mutate(
    old_rank = rank(-score, ties.method = "first"),
    new_rank = rank(-estimate_intercept, ties.method = "first"),
    delta = new_rank - old_rank
  )

plot_df <- blupdf_ranked |>
  select(benchmark, model_name, old_rank, new_rank, delta) |>
  pivot_longer(c(old_rank, new_rank), names_to = "type", values_to = "rank") |>
  mutate(type = recode(type, old_rank = "Old", new_rank = "New"))

focus_b <- c("taubench_airline", "swebench_verified_mini")

rp <- plot_df |>
  filter(benchmark %in% focus_b) |>
  ggplot(aes(x = type, y = rank, group = model_name, color = delta)) +
  geom_line(alpha = 1) +
  geom_point(size = 2) +
  geom_text_repel(
    data = filter(plot_df, benchmark %in% focus_b, type == "New"),
    aes(label = model_name),
    direction = "y", hjust = 0, nudge_x = 0.1, size = 3, segment.alpha = 0.3
  ) +
  scale_y_reverse() +
  xlim("Old", "New") +
  scale_color_gradient2_tableau() +
  facet_wrap(~benchmark, ncol = 2) +
  ggpubr::theme_pubr() +
  theme(
    text = element_text(family = "times", size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "Rank (1 = best)", color = "Rank Change\n(New - Old)")

ggsave(here("figures", "benchmark_ranking_changes_tau_swe.pdf"), rp, width = 10.5, height = 4.57, device = cairo_pdf)

plogis_blup <- function(x) stats::plogis(x)

pointp <- blupdf |>
  filter(benchmark %in% focus_b) |>
  mutate(
    estimate_intercept = plogis_blup(estimate_intercept),
    est_error_intercept = est_error_intercept * estimate_intercept * (1 - estimate_intercept)
  ) |>
  ggplot(aes(x = score, y = estimate_intercept, color = model_name)) +
  geom_pointrange(
    aes(
      ymin = estimate_intercept - est_error_intercept,
      ymax = estimate_intercept + est_error_intercept
    ),
    position = "jitter"
  ) +
  facet_wrap(~benchmark, scales = "free_x") +
  labs(y = "Estimated Score", x = "Original Mean Score") +
  paper_theme(13) +
  theme(legend.position = "none", panel.border = element_rect(color = "grey80", fill = NA, linewidth = 1))

ggsave(here("figures", "blups_vs_raw_scoresSE_tau_swe_nolegend.pdf"), pointp, width = 8.74, height = 5.32, device = cairo_pdf)

invisible(list(mvar = mvar, gbms = gbms, combdf = combdf, blupdf = blupdf))
}
} # RUN_ESTIMATION
