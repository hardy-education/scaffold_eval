# Agent leaderboard reliability (random-effects models)

Code accompanying the NeurIPS draft on **agent evaluation leaderboard reliability**: mixed-effects and Bayesian variance decompositions of HAL benchmark scores, plus nonparametric (`energy::disco`) baselines.

- Paper: [`NeurIPS_Agent_Leaderboard_Reliability_draft.pdf`](NeurIPS_Agent_Leaderboard_Reliability_draft.pdf)
- Narrative walkthrough: [`docs/reliability_models.qmd`](docs/reliability_models.qmd) (render with Quarto)

## Repository layout

| Path | Purpose |
|------|---------|
| [`data/hal_parsed_response_matrix7.csv`](data/hal_parsed_response_matrix7.csv) | Run-level scores (see [`data/README.md`](data/README.md)) |
| [`R/utils/`](R/utils/) | Shared loaders, formulas, Ep² helpers, `disco`, plotting |
| [`R/analysis/01_pooled_all_benchmarks.R`](R/analysis/01_pooled_all_benchmarks.R) | Full hierarchy pooled across benchmarks |
| [`R/analysis/02_per_benchmark.R`](R/analysis/02_per_benchmark.R) | Separate model per benchmark |
| [`R/analysis/03_per_benchmark_model_agent.R`](R/analysis/03_per_benchmark_model_agent.R) | Collapsed model×agent facet per benchmark |
| [`figures/`](figures/) | Generated PDF figures |

## Quick start

1. Install R ≥ 4.2 and dependencies:

   ```r
   source("install.R")
   ```

   Bayesian models require [CmdStan](https://mc-stan.org/cmdstanr/) via `cmdstanr::install_cmdstan()`.

2. Explore without refitting (fast):

   ```r
   source(here::here("R", "setup.R"))
   source(here::here("R", "utils", "source_utils.R"))
   df <- load_hal_data()
   ```

3. Render the Quarto walkthrough (does **not** run analysis scripts):

   ```bash
   quarto render docs/reliability_models.qmd
   ```

4. To refit models, set flags **before** sourcing an analysis script:

   ```r
   RUN_ESTIMATION <- TRUE   # required; pooled script ~30+ min for LME/GLME alone
   FIT_MODELS <- FALSE      # TRUE refits brms (very slow)
   MAKE_PLOTS <- TRUE
   source(here::here("R", "analysis", "02_per_benchmark.R"))
   ```

## Configuration flags

Each analysis script defaults to **no estimation** (`RUN_ESTIMATION = FALSE`):

- `RUN_ESTIMATION` — `TRUE` runs lme4/glmer/brms fits; `FALSE` exits immediately (safe default).
- `FIT_MODELS` — when estimating, `TRUE` refits brms; `FALSE` loads cached `data/*.rds` if present.
- `RUN_DISCO` — `FALSE` loads cached `data/results/disco_*.csv`; `TRUE` runs `energy::disco` (~1 hr).
- `MAKE_PLOTS` — `FALSE` skips figure generation.

## Scripts ↔ figures

| Script | Example outputs |
|--------|-----------------|
| `01_pooled_all_benchmarks.R` | `figures/bayes_variance_decomp_full.pdf`, `benchmark_ranking_changes_tau_swe.pdf`, BLUP diagnostics |
| `02_per_benchmark.R` | `figures/bayes_CI_per_bench.pdf`, `bayes_CI_per_bench_mvamp.pdf` |
| `03_per_benchmark_model_agent.R` | `figures/mod_agent_reliab_bar.pdf`, `bayes_CI_mod-age_per_bench.pdf` |

## License

MIT — see [`LICENSE`](LICENSE).
