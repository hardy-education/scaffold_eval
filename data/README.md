# Data directory

## `hal_parsed_response_matrix7.csv`

Run-level outcomes from the HAL agent evaluation matrix (7 benchmarks).

| Column | Description |
|--------|-------------|
| `run_id` | Unique run identifier |
| `benchmark` | Benchmark name (e.g. `swebench_verified_mini`) |
| `task_id` | Item / task within benchmark |
| `agent_name` | Agent scaffold or system |
| `model_name` | Base model (may include provider prefix before `:`) |
| `reasoning_effort` | Reasoning-effort suffix concatenated to `model_name` in code |
| `score` | Binary or proportion correct (0–1) |

Loaded via `load_hal_data()` in [`R/utils/load_hal_data.R`](../R/utils/load_hal_data.R).

## Cached outputs (generated locally)

| Path | Contents |
|------|----------|
| `data/brmfit_*.rds`, `data/gbrmfit_*.rds` | Cached `brms` fits (`FIT_MODELS=FALSE` reads these) |
| `data/results/disco_*.csv` | Nonparametric variance decomposition tables |

These artifacts are gitignored by default; regenerate with `FIT_MODELS=TRUE` / `RUN_DISCO=TRUE` in the analysis scripts.
