#' Nonparametric variance decomposition via energy::disco.
#'
#' @param df Data with score, benchmark, task_id, model_name, agent_name.
#' @param benchmark If set, subset to this benchmark; if NULL, use all rows (pooled).
#' @param collapse_model_agent Use collapsed `model` column instead of separate facets.
#' @param alpha Disco index parameter.
disco_decompose <- function(
    df,
    benchmark = NULL,
    collapse_model_agent = FALSE,
    alpha = 1) {
  dat <- df
  if (!is.null(benchmark)) {
    dat <- dplyr::filter(dat, benchmark == .env$benchmark)
  }
  dat <- dat |>
    dtplyr::lazy_dt() |>
    dplyr::filter(!is.na(score)) |>
    dplyr::mutate(resp = as.numeric(score))

  if (collapse_model_agent) {
    dat <- dat |>
      dplyr::mutate(
        item = as.factor(interaction(benchmark, task_id)),
        model = as.factor(model)
      ) |>
      dplyr::select(-dplyr::any_of(c(
        "run_id", "model_name", "agent_name", "score", "task_id",
        "reasoning_effort", "benchmark"
      )))
  } else if (!is.null(benchmark)) {
    dat <- dat |>
      dplyr::mutate(
        item = as.factor(interaction(benchmark, task_id)),
        agent = as.factor(agent_name),
        model = as.factor(model_name),
        item_agent = forcats::fct_drop(interaction(benchmark, item, agent)),
        item_model = forcats::fct_drop(interaction(benchmark, item, model)),
        model_agent = forcats::fct_drop(interaction(model, agent))
      ) |>
      dplyr::select(-dplyr::any_of(c(
        "run_id", "model_name", "agent_name", "score", "task_id",
        "reasoning_effort", "benchmark"
      )))
  } else {
    dat <- dat |>
      dplyr::mutate(
        item = as.factor(interaction(benchmark, task_id)),
        agent = as.factor(agent_name),
        model = as.factor(model_name),
        benchmark_f = as.factor(benchmark),
        item_agent = forcats::fct_drop(interaction(benchmark, item, agent)),
        item_model = forcats::fct_drop(interaction(benchmark, item, model)),
        model_agent = forcats::fct_drop(interaction(model, agent)),
        bench_agent = forcats::fct_drop(interaction(benchmark, agent)),
        bench_model = forcats::fct_drop(interaction(benchmark, model)),
        bench_model_agent = forcats::fct_drop(interaction(benchmark, model, agent)),
        item_model_agent = forcats::fct_drop(interaction(benchmark, item, model, agent))
      ) |>
      dplyr::select(-dplyr::any_of(c(
        "run_id", "model_name", "agent_name", "score", "task_id", "reasoning_effort"
      )))
  }

  dat_dt <- dplyr::collect(dat)
  resp <- dat_dt[, "resp", drop = FALSE]
  factors <- dat_dt[, setdiff(names(dat_dt), "resp"), drop = FALSE]

  disco_fit <- energy::disco(
    resp |> data.table::as.data.table(),
    factors = factors |> data.table::as.data.table(),
    distance = FALSE,
    index = alpha,
    R = 0
  )

  disco_fit$stats |>
    tibble::as_tibble() |>
    dplyr::select(-dplyr::any_of("p-value")) |>
    dplyr::mutate(
      facet = disco_fit$factor.names,
      Withins = Within,
      Within = disco_fit$within,
      Total = disco_fit$total,
      Ep2 = Trt / (Trt + Withins),
      SNR = Ep2 / (1 - Ep2),
      Ep2hat = (Trt / df1) / ((Trt / df1) + Withins / df2),
      SNRhat = Ep2hat / (1 - Ep2hat),
      pct = Trt / (Total - Within),
      mean = Trt / df1,
      var = Withins / df2,
      benchmark = benchmark,
      alpha = alpha
    ) |>
    dplyr::select(facet, Ep2, SNR, dplyr::everything())
}

run_disco_by_benchmark <- function(df, benches, alpha = 1) {
  purrr::map(benches, \(b) {
    message("disco: ", b)
    disco_decompose(df, benchmark = b, alpha = alpha) |>
      dplyr::mutate(
        benchmark = b,
        n_items = dplyr::n_distinct(interaction(df$benchmark, df$task_id)[df$benchmark == b])
      )
  }) |>
    dplyr::bind_rows()
}
