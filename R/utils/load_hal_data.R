#' Load and preprocess the HAL agent evaluation response matrix.
#'
#' @param path Path to CSV (default: project data file).
#' @param collapse_model_agent If TRUE, add a `model` column as paste of model and agent.
#' @return A tibble with normalized `model_name` and optional `model` facet.
load_hal_data <- function(
    path = here::here("data", "hal_parsed_response_matrix7.csv"),
    collapse_model_agent = FALSE) {
  out <- data.table::fread(path) |>
    dplyr::mutate(
      model_name = stringr::str_split_i(model_name, ":", 1),
      model_name = dplyr::case_when(
        model_name == "DeepSeek-R1" ~ "deepseek-r1",
        model_name == "claude-sonnet-4.5" ~ "claude-sonnet-4-5",
        TRUE ~ model_name
      ),
      model_name = stringr::str_c(model_name, reasoning_effort)
    )

  if (collapse_model_agent) {
    out <- dplyr::mutate(out, model = stringr::str_c(model_name, agent_name))
  }

  tibble::as_tibble(out)
}

bench_item_counts <- function(df) {
  df |>
    dplyr::select(benchmark, task_id) |>
    dplyr::distinct() |>
    dplyr::group_by(benchmark) |>
    dplyr::summarize(n = dplyr::n_distinct(task_id), .groups = "drop")
}
