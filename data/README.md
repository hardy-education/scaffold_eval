# Data

## `hal_response_matrix.csv`

The analysis dataset: **29,923 agent rollouts** drawn from the nine benchmarks
published on the Holistic Agent Leaderboard (HAL). One row per rollout, with a
binary outcome.

| Column | Description |
|---|---|
| `run_id` | Identifier of the HAL evaluation run the rollout belongs to |
| `benchmark` | *b* — one of the nine benchmarks listed below |
| `task_id` | *i* — task (item), unique within a benchmark, not across benchmarks |
| `agent_name` | *a* — agent scaffold |
| `model_name` | *m* — base LLM with reasoning effort appended (see below) |
| `reasoning_effort` | The effort setting, retained for reference |
| `score` | *y* ∈ {0, 1} — whether the rollout solved the task |

### What counts as a "model"

Throughout the paper and this code, a model is a **base LLM paired with a
reasoning-effort configuration**: `gpt_5_2025_08_07high` and
`gpt_5_2025_08_07minimal` are distinct levels of the model facet, because they
are distinct systems a leaderboard would list separately. That pairing is
already applied in the shipped file, and `normalize_hal_columns()` in
[`../R/data.R`](../R/data.R) applies it idempotently to a fresh HAL export.

### Design

| Benchmark | Models | Tasks | Scaffolds | Model–scaffold pairs | Rollouts | Mean score |
|---|---:|---:|---:|---:|---:|---:|
| AssistantBench | 18 | 33 | 2 | 30 | 1,056 | 0.04 |
| CORE-Bench Hard | 34 | 45 | 3 | 56 | 2,880 | 0.26 |
| GAIA | 21 | 165 | 2 | 34 | 5,926 | 0.46 |
| Online-Mind2Web | 13 | 300 | 2 | 23 | 6,900 | 0.33 |
| SciCode | 17 | 65 | 3 | 37 | 2,464 | 0.03 |
| ScienceAgentBench | 19 | 102 | 2 | 25 | 2,550 | 0.21 |
| SWE-bench Verified Mini | 24 | 50 | 2 | 26 | 1,444 | 0.31 |
| τ-bench Airline | 23 | 50 | 3 | 42 | 2,099 | 0.42 |
| USACO | 13 | 307 | 2 | 14 | 4,604 | 0.41 |

Across the whole matrix: 9 benchmarks, 1,117 benchmark-specific tasks, 54
models, 13 scaffolds. Regenerate this table with
`Rscript scripts/00_design_summary.R`.

### Sparsity, and why the design still identifies the components

Only about **3.7%** of the fully crossed benchmark × task × model × scaffold
design is observed. The decomposition does not need complete crossing; it needs
enough overlap to tell persistent model variation apart from variation
attributable to benchmarks, scaffolds, and their interactions. The relevant
anchors, all checkable from the data:

- **9 models are evaluated on every benchmark.** These separate the model main
  effect from benchmark-conditioned model performance: observations that share
  *M* while differing in *BM*.
- **One scaffold (`hal_generalist`) spans 8 of the 9 benchmarks**, and a second
  shared scaffold connects AssistantBench to the rest. Without these, a
  scaffold used on a single benchmark would be indistinguishable from a
  benchmark–scaffold interaction.
- **12 of 54 models are only ever seen under one scaffold**, and roughly half
  the scaffolds appear on a single benchmark. These levels are estimated
  through the hierarchical model's exchangeability assumptions, with the
  resulting uncertainty carried into the posterior rather than hidden.
- **Only 4.3% of (benchmark, task, model, scaffold) cells are replicated.** The
  four-way interaction therefore cannot be separated from cell-level execution
  noise, and both are folded into the terminal component σ²_{BIMA,e}.

`scripts/00_design_summary.R` prints all of these, and the paper's Appendix E
develops the identifiability argument they support.

## Cached outputs

Model fits and long-running decomposition results are written to `outputs/`
and are **not** tracked in git — the pooled `brms` fit alone is over 100 MB.
Regenerate them with `REFIT=TRUE Rscript scripts/01_fit_models.R`, or drop
previously saved `.rds` files into `outputs/fits/` to skip re-estimation.

## Provenance

Rollout outcomes come from the public HAL release (Kapoor et al., 2024). This
repository redistributes only the parsed outcome matrix — task identifiers and
binary scores — not the rollout traces themselves. For AssistantBench, GAIA,
SciCode, and SWE-bench Verified Mini, HAL evaluates a subset of the full
benchmark, so task counts here are smaller than in the original releases.
