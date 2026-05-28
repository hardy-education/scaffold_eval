#' Fit or load a cached brms model.
fit_or_load_brm <- function(
    data,
    formula,
    file,
    fit_models = FALSE,
    family = NULL,
    save_pars = brms::save_pars(all = TRUE),
    ...) {
  args <- c(
    list(
      formula = formula,
      data = data,
      file = file,
      file_refit = if (fit_models) "always" else "on_change",
      save_pars = save_pars
    ),
    brm_defaults(),
    list(...)
  )
  if (!is.null(family)) args$family <- family

  if (file.exists(file) && !fit_models) {
    return(readRDS(file))
  }
  if (!fit_models) {
    stop(
      "Cached brms fit not found at ", file,
      ". Set FIT_MODELS <- TRUE to estimate (slow), or copy a saved .rds into data/."
    )
  }
  do.call(brms::brm, args)
}

#' Loop per-benchmark mixed models.
fit_per_benchmark <- function(
    df,
    benches,
    fit_fun,
    summary_fun = varcorr_summary) {
  fits <- stats::setNames(vector("list", length(benches)), benches)
  sums <- stats::setNames(vector("list", length(benches)), benches)
  for (b in benches) {
    message("fitting: ", b)
    fit <- df |> dplyr::filter(benchmark == b) |> fit_fun()
    fits[[b]] <- fit
    sums[[b]] <- summary_fun(fit)
  }
  list(fits = fits, summaries = dplyr::bind_rows(sums, .id = "benchmark"))
}
