#' Summarize a posterior draw vector (mean, HDI, Bayes factors, ROPE).
make_draws_stats <- function(drawstat, roperange = c(0, 0.005)) {
  summary(drawstat) |>
    dplyr::bind_cols(
      bayestestR::p_significance(drawstat, threshold = 0.05) |>
        dplyr::select(-Parameter),
      bayestestR::p_map(drawstat) |>
        dplyr::rename(pval = p_MAP) |>
        dplyr::select(-Parameter),
      bayestestR::p_map(drawstat, method = "KernSmooth") |>
        dplyr::select(-Parameter),
      bayestestR::map_estimate(drawstat) |>
        dplyr::select(-Parameter),
      bayestestR::equivalence_test(drawstat, range = roperange) |>
        dplyr::select(-Parameter)
    )
}

#' Posterior variance shares from brms SD parameters.
summarize_brm_variance_draws <- function(draws_rvars) {
  sd_draws <- draws_rvars[stringr::str_starts(names(draws_rvars), "sd|sig")]
  denom <- Reduce(function(acc, nxt) acc + nxt^2, sd_draws)
  lapply(sd_draws, function(x) make_draws_stats(x^2 / denom)) |>
    dplyr::bind_rows(.id = "var")
}

convert_to_percent_of_var_draws <- function(draws, term) {
  draws[[term]]^2 / sum(draws[[term]]^2)
}
