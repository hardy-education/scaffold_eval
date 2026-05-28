# Colors and labels for variance-component plots (paper figures).

component_cols <- c(
  agent = "#009E73",
  item_agent = "#F0E442",
  model_agent = "#0072B2",
  item_model = "#D55E00",
  item = "#E69F00",
  model = "#56B4E9",
  model_name = "#56B4E9",
  sigma = "#000000",
  benchmark = "#CC79A7",
  benchmark_agent = "#009E73",
  benchmark_item_agent = "#F0E442",
  benchmark_model_agent = "#0072B2",
  benchmark_item_model = "#D55E00",
  benchmark_item = "#E69F00",
  benchmark_model = "#56B4E9",
  benchmark_item_model_agent = "#000000"
)

component_labs <- c(
  agent = "Agent",
  benchmark_agent = "Benchmark-Agent",
  benchmark = "Benchmark",
  benchmark_model_agent = "Benchmark-Model-Agent",
  item_agent = "Item-Agent",
  benchmark_item_agent = "Item-Agent",
  item = "Item",
  benchmark_item = "Item",
  model_agent = "Model-Agent",
  item_model = "Model-Item",
  benchmark_item_model = "Model-Item",
  model = "Model",
  benchmark_model = "Benchmark-Model",
  benchmark_item_model_agent = "Item-Model-Agent",
  sigma = "Residual"
)

component_alphas <- c(
  agent = 1,
  item_agent = 1,
  model_agent = 1,
  item_model = 1,
  item = 1,
  model = 1,
  sigma = 1,
  benchmark = 1,
  benchmark_agent = 0.5,
  benchmark_item_agent = 1,
  benchmark_model_agent = 0.5,
  benchmark_item_model = 1,
  benchmark_item = 1,
  benchmark_model = 0.5,
  benchmark_item_model_agent = 0.5
)

paper_theme <- function(base_size = 15) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(text = ggplot2::element_text(family = "times"))
}
