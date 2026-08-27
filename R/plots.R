# R/plots.R ----------------------------------------------------------------
#
# Shared theme, palette, and one builder per figure in the paper.
#
# Each `plot_*()` takes a tidy data frame produced by the analysis modules and
# returns a ggplot object; nothing here computes statistics. Analysis scripts
# therefore read as "compute, then draw", and a figure can be re-styled without
# touching an estimate.

# ---- Theme and palette -----------------------------------------------------

#' Shared figure theme.
#'
#' Uses `PLOT_FONT` from config.R, which defaults to R's device-independent
#' "serif" family so the figures build identically on a machine without any
#' particular font installed.
paper_theme <- function(base_size = 15, legend = "right") {
  family <- if (nzchar(PLOT_FONT)) PLOT_FONT else ""
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = family),
      legend.position = legend,
      panel.border = ggplot2::element_rect(colour = "grey85", fill = NA, linewidth = 0.4),
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.7))
    )
}

# Okabe-Ito palette, keyed by variance component. Interaction terms involving
# the benchmark are drawn at reduced alpha so that the pooled decomposition
# reads as "within-benchmark structure" versus "across-benchmark structure".
COMPONENT_COLOURS <- c(
  M = "#56B4E9", A = "#009E73", I = "#E69F00", B = "#CC79A7",
  IM = "#D55E00", IA = "#F0E442", MA = "#0072B2",
  BM = "#56B4E9", BA = "#009E73", BMA = "#0072B2", IMA = "#000000",
  e = "#000000"
)

COMPONENT_ALPHAS <- c(
  M = 1, A = 1, I = 1, B = 1, IM = 1, IA = 1, MA = 1,
  BM = 0.55, BA = 0.55, BMA = 0.55, IMA = 0.55, e = 1
)

# Two-colour scale for the two objects of measurement.
OBJECT_COLOURS <- c(model = "#D55E00", system = "#0072B2")
OBJECT_LABELS  <- c(model = "Model", system = "Scaffold-Model")

#' Replace benchmark codes with display names in a facetting variable.
#'
#' Defaults to the abbreviated names: a nine-panel row is too narrow for the
#' full ones, and a silently clipped strip label is worse than a short one.
label_benchmarks <- function(short = TRUE) {
  ggplot2::as_labeller(if (short) BENCHMARK_SHORT else BENCHMARK_LABELS)
}

#' Save a figure at the size used in the paper.
save_figure <- function(plot, name, width, height) {
  path <- file.path(PATHS$figures, paste0(name, ".pdf"))
  ggplot2::ggsave(path, plot, width = width, height = height,
                  device = grDevices::cairo_pdf)
  message("wrote ", path)
  invisible(path)
}

#' Save a results table as both CSV and a booktabs LaTeX fragment.
#'
#' The CSV is the machine-readable record; the `.tex` fragment is what the
#' manuscript includes, so a number in the paper can be traced to the script
#' that produced it without hand-transcription.
save_table <- function(dat, name, digits = 3, caption = NULL) {
  csv <- file.path(PATHS$tables, paste0(name, ".csv"))
  readr::write_csv(dat, csv)

  tex <- file.path(PATHS$tables, paste0(name, ".tex"))
  kableExtra::kable(dat, digits = digits, format = "latex",
                    booktabs = TRUE, caption = caption) |>
    kableExtra::kable_styling(latex_options = c("hold_position", "scale_down")) |>
    writeLines(tex)

  message("wrote ", csv, " and ", basename(tex))
  invisible(dat)
}


# ---- Figure 1: inter-scaffold reliability ----------------------------------

#' Inter-scaffold reliability by benchmark (paper Figure 1).
#'
#' Each bar is rho_AA' from Eq. 5: how much two independently sampled scaffolds
#' would agree about the model ordering, having averaged over all of the
#' benchmark's tasks. Low bars mean scaffold choice changes who looks best.
plot_inter_scaffold <- function(dat) {
  ggplot(dat, aes(x = benchmark, y = median)) +
    geom_col(fill = OBJECT_COLOURS[["model"]], alpha = 0.55, width = 0.7) +
    geom_pointrange(aes(y = mean, ymin = low, ymax = high), linewidth = 0.5) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_discrete(labels = BENCHMARK_LABELS) +
    labs(x = "Benchmark", y = "Inter-Scaffold Reliability",
         title = "Inter-Scaffold Reliability Across Benchmark Tasks") +
    paper_theme(legend = "none") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.x = element_text(angle = 30, hjust = 1))
}


# ---- Figure 2: per-benchmark D-study, two objects --------------------------

#' Reliability against task count, by object of measurement (paper Figure 2).
#'
#' The dashed line is the task-only ceiling of Eq. 8: where the curve would
#' land with unlimited tasks of the same construction. The gap between the two
#' curves is the variance that scaffold choice contributes -- error when
#' ranking models, signal when ranking deployable systems.
plot_bench_reliability <- function(dat, ceilings = NULL) {
  p <- ggplot(dat, aes(x = n_tasks, y = median, colour = object, fill = object)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, linewidth = 0) +
    geom_line(linewidth = 0.8)

  if (!is.null(ceilings)) {
    p <- p + geom_hline(aes(yintercept = median), data = ceilings,
                        linetype = "dashed", colour = OBJECT_COLOURS[["model"]],
                        linewidth = 0.4)
  }

  p +
    scale_x_log10() +
    scale_colour_manual(values = OBJECT_COLOURS, labels = OBJECT_LABELS, name = "Object") +
    scale_fill_manual(values = OBJECT_COLOURS, labels = OBJECT_LABELS, name = "Object") +
    facet_wrap(~benchmark, nrow = 1, scales = "free_x", labeller = label_benchmarks()) +
    labs(x = "Number of Agentic Tasks",
         y = expression("Est. Rank Reliability: E" * hat(rho)^2)) +
    paper_theme()
}


# ---- Figure 3: per-benchmark signal-to-noise -------------------------------

#' Signal-to-noise against task count (paper Figure 3).
#'
#' The shaded band is the external limit-of-detection reference region. It is a
#' measurement-design diagnostic rather than a pass/fail rule: decisions with
#' larger consequences warrant stronger evidence than its lower edge.
plot_bench_snr <- function(dat, y_cap = 6) {
  dat <- dat |> mutate(across(c(median, low, high, q95), \(x) pmin(x, y_cap)))
  # Finite bounds rather than +/-Inf: the x scale is logarithmic, and log10(-Inf)
  # is a NaN warning on every render.
  xr <- range(dat$n_tasks, na.rm = TRUE)

  ggplot(dat, aes(x = n_tasks, y = median)) +
    annotate("rect", xmin = xr[1], xmax = xr[2],
             ymin = LOD_BAND[["low"]], ymax = LOD_BAND[["high"]],
             fill = "#F0E442", alpha = 0.45) +
    geom_ribbon(aes(ymin = low, ymax = high),
                fill = OBJECT_COLOURS[["model"]], alpha = 0.25, linewidth = 0) +
    geom_line(colour = OBJECT_COLOURS[["model"]], linewidth = 0.8) +
    scale_x_log10() +
    facet_wrap(~benchmark, nrow = 1, scales = "free", labeller = label_benchmarks()) +
    ylim(0, NA) +
    labs(x = "Number of Agentic Tasks",
         y = expression("Signal-to-Noise Ratio: S/N(" * delta * ")")) +
    paper_theme(legend = "none")
}


# ---- Figure 4: leaderboard D-study -----------------------------------------

#' Leaderboard reliability against total task budget (paper Figure 4a).
#'
#' The x-axis is the evaluation budget actually spent, so the two curves answer
#' a design question directly: given N tasks, is it better to draw them all
#' from one benchmark or to spread them across nine?
plot_leaderboard_reliability <- function(dat, which_nb = c(1, 9)) {
  dat |>
    filter(n_benchmarks %in% which_nb, total_tasks >= 5) |>
    mutate(benches = factor(n_benchmarks,
                            labels = paste(which_nb, ifelse(which_nb == 1, "bench", "benches")))) |>
    ggplot(aes(x = total_tasks, y = median, colour = benches, fill = benches)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, linewidth = 0) +
    geom_line(linewidth = 1) +
    scale_x_log10() +
    labs(x = "Total Number of Agentic Tasks",
         y = expression("Estimated Overall Rank Reliability: E" * hat(rho)^2),
         colour = "Benchmarks\nSampled", fill = "Benchmarks\nSampled") +
    paper_theme(legend = "bottom")
}

#' Leaderboard signal-to-noise against total task budget (paper Figure 4b).
#'
#' One curve per battery size. Where a curve crosses the reference band is the
#' number of benchmarks a battery needs before model differences are detectable
#' at all.
plot_leaderboard_snr <- function(dat) {
  dat |>
    filter(total_tasks >= 5) |>
    mutate(benches = factor(n_benchmarks,
                            labels = paste(sort(unique(n_benchmarks)), "bench."))) |>
    ggplot(aes(x = total_tasks, y = median, colour = benches)) +
    annotate("rect",
             xmin = min(dat$total_tasks[dat$total_tasks >= 5], na.rm = TRUE),
             xmax = max(dat$total_tasks, na.rm = TRUE),
             ymin = LOD_BAND[["low"]], ymax = LOD_BAND[["high"]],
             fill = "#F0E442", alpha = 0.35) +
    geom_hline(yintercept = LOD_BAND[["mid"]], linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.9) +
    scale_colour_viridis_d(option = "viridis", direction = -1) +
    labs(x = "Total Number of Agentic Tasks",
         y = expression("Estimated Signal-to-Noise Ratio: S/N(" * delta * ")"),
         colour = "Benchmarks\nSampled") +
    paper_theme()
}


# ---- Figure 5: variance decomposition bar ----------------------------------

#' Stacked variance-decomposition bar (paper Figure 5).
#'
#' Reads left to right as the share of latent variance each component carries.
#' The point of the figure is the size of the model bar relative to the
#' condition-specific ones beside it.
#'
#' Uses `median_normalized` rather than `median`: posterior medians of shares
#' do not sum to one, and a stacked bar has to. Quote `median` for any single
#' component you report in text.
plot_variance_decomposition <- function(shares, label_size = 3.5) {
  shares |>
    mutate(component = factor(component, levels = names(COMPONENT_COLOURS)),
           median = median_normalized) |>
    arrange(component) |>
    ggplot(aes(x = "", y = median, fill = component, alpha = component)) +
    geom_col(width = 0.6) +
    ggrepel::geom_label_repel(
      aes(label = label), position = position_stack(vjust = 0.5),
      direction = "y", size = label_size, alpha = 1, fill = "white",
      label.size = 0.15, min.segment.length = 0, max.overlaps = Inf
    ) +
    coord_flip() +
    scale_fill_manual(values = COMPONENT_COLOURS) +
    scale_alpha_manual(values = COMPONENT_ALPHAS) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(y = "Proportion of Variance Explained") +
    paper_theme(legend = "none") +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.border = element_blank())
}


# ---- Figure 6: published rank vs posterior rank ----------------------------

#' Published rank against posterior rank, with credible intervals (Figure 6).
#'
#' Points far from the diagonal are models whose position depends on whether
#' task and scaffold variation is credited to the model. Intervals that overlap
#' the diagonal are ranks the data do not resolve.
plot_rank_comparison <- function(dat, x = "published_pctrank",
                                 lo = "rank_q2_5", hi = "rank_q97_5",
                                 y = "rank_median") {
  ggplot(dat, aes(x = .data[[x]], y = .data[[y]], colour = model_name)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    geom_pointrange(aes(ymin = .data[[lo]], ymax = .data[[hi]]),
                    position = position_jitter(width = 0.01, height = 0),
                    size = 0.35, linewidth = 0.4) +
    facet_wrap(~benchmark, labeller = label_benchmarks()) +
    guides(colour = guide_legend(ncol = 4)) +
    labs(x = "Original/Published Rank", y = "Posterior Rank (Median, 95% CI)") +
    paper_theme(legend = "bottom") +
    theme(legend.title = element_blank(),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 7))
}


# ---- Figure 7: rank-shift slope chart --------------------------------------

#' Published versus adjusted ordering as a slope chart (paper Figure 7).
#'
#' Colour encodes the size and direction of the move. Models paired with
#' unusually favourable scaffolds tend to fall once scaffold compatibility is
#' separated from model capability.
plot_rank_shift <- function(dat, ncol = 2) {
  long <- dat |>
    select(benchmark, model_name, delta,
           Published = published_rank, Adjusted = rank_median) |>
    tidyr::pivot_longer(c(Published, Adjusted),
                        names_to = "ordering", values_to = "rank") |>
    mutate(ordering = factor(ordering, levels = c("Published", "Adjusted")))

  ggplot(long, aes(x = ordering, y = rank, group = model_name, colour = delta)) +
    geom_line(alpha = 0.9) +
    geom_point(size = 1.6) +
    ggrepel::geom_text_repel(
      data = filter(long, ordering == "Adjusted"),
      aes(label = model_name), direction = "y", hjust = 0,
      nudge_x = 0.12, size = 2.6, segment.alpha = 0.3, max.overlaps = Inf
    ) +
    scale_y_reverse() +
    scale_x_discrete(expand = expansion(mult = c(0.15, 0.9))) +
    ggthemes::scale_color_gradient2_tableau() +
    facet_wrap(~benchmark, ncol = ncol, scales = "free_y", labeller = label_benchmarks()) +
    labs(y = "Rank (1 = best)", colour = "Rank change") +
    paper_theme(legend = "none") +
    theme(axis.title.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.border = element_blank())
}


# ---- Model vs. scaffold contrast -------------------------------------------

#' Posterior of the model-minus-scaffold variance contrast (Section 5.1).
#'
#' Mass to the left of zero is posterior probability that scaffold-related
#' variation exceeds model-related variation.
plot_model_scaffold_contrast <- function(contrast_row) {
  d <- tibble(x = contrast_row$draws[[1]])
  dens <- stats::density(d$x)
  dd <- tibble(x = dens$x, y = dens$y)

  ggplot(dd, aes(x, y)) +
    geom_area(data = filter(dd, x < 0), fill = OBJECT_COLOURS[["system"]], alpha = 0.6) +
    geom_area(data = filter(dd, x > 0), fill = OBJECT_COLOURS[["model"]], alpha = 0.5) +
    geom_line(linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    labs(
      x = "Proportion of Variance: Var(Model) - Var(Scaffold)",
      y = "Posterior Density",
      subtitle = glue::glue(
        "P(scaffold > model) = {scales::percent(contrast_row$p_scaffold_exceeds_model, 0.1)}"
      )
    ) +
    paper_theme(legend = "none")
}


# ---- Ablations -------------------------------------------------------------

#' Leave-one-benchmark-out D-study curves (paper Appendix D.4).
#'
#' One line per withheld benchmark. A benchmark whose removal moves the curve
#' is one the battery depends on.
plot_loo_curves <- function(dat, statistic = "ep2") {
  dat |>
    filter(statistic == .env$statistic, total_tasks >= 5) |>
    ggplot(aes(x = total_tasks, y = median,
               colour = removed, group = interaction(removed, n_benchmarks))) +
    geom_line(linewidth = 0.7) +
    scale_x_log10() +
    scale_colour_discrete(labels = c(BENCHMARK_LABELS, none = "None (full data)")) +
    facet_wrap(~n_benchmarks, labeller = as_labeller(\(x) paste(x, "benchmarks sampled"))) +
    labs(
      x = "Total Number of Agentic Tasks",
      y = if (statistic == "ep2") {
        expression("Estimated Reliability: E" * hat(rho)^2)
      } else {
        expression("Signal-to-Noise Ratio: S/N(" * delta * ")")
      },
      colour = "Benchmark removed"
    ) +
    paper_theme()
}


#' Stability of variance shares across an ablation (Tables 6 and 7).
plot_ablation_stability <- function(dat) {
  ggplot(dat, aes(x = stats::reorder(label, share_mean), y = share)) +
    geom_boxplot(outlier.size = 0.6, linewidth = 0.3, fill = "grey92") +
    coord_flip() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Proportion of variance across refits") +
    paper_theme(legend = "none")
}
