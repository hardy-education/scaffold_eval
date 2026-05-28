#' Harmonize facet names across lme4, brms, and disco outputs.
rename_lme_facet <- function(grp) {
  grp |>
    stringr::str_replace_all("task_id", "item") |>
    stringr::str_remove_all(stringr::regex("sd_|__Intercept|_name|_id")) |>
    stringr::str_replace_all(":", "_") |>
    stringr::str_replace_all("Residual", "sigma")
}

rename_brm_facet <- function(var) {
  var |>
    stringr::str_replace_all("task_id", "item") |>
    stringr::str_remove_all(stringr::regex("sd_|__Intercept|_name|_id")) |>
    stringr::str_replace_all(":", "_")
}

rename_disco_facet <- function(facet) {
  facet |>
    stringr::str_replace_all("item", "bench_item") |>
    stringr::str_replace_all("bench_", "benchmark_")
}

finalize_facet <- function(facet) {
  dplyr::if_else(facet == "benchmark_item_model_agent", "sigma", facet)
}
