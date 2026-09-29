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

---

## `harbor_index_data.csv`

The **corroboration dataset**: 1,476 trial-level outcomes from the Harbor
Index, a curated meta-benchmark distilled from roughly 6,000 candidate tasks
across 54 benchmarks down to 82 tasks spanning 29 benchmarks and eight
domains, filtered for difficulty and audited for correctness.

| Column | Description |
|---|---|
| `benchmark` | Source benchmark the task came from |
| `task_id` | Task identifier |
| `model` | LLM under test (9 models) |
| `agent` | Scaffold / harness (4 scaffolds) |
| `outcome` | Verifier–judge classification: `TP`, `TN`, `FP`, `FN` |
| `resp` | Raw verifier reward (may be fractional or missing) |
| `domain` | Topical grouping of the benchmark |
| `judged_score` | **Adjudicated binary outcome used in the analysis** |

### Why this dataset is here

Its sparsity runs the other way from HAL's. HAL has many tasks, but which
models and scaffolds appear varies from benchmark to benchmark, and the pooled
fit leans on a few shared anchors to connect them. Harbor has few tasks, but
runs the **same 9 models and 4 scaffolds on every benchmark**, so
benchmark × model and benchmark × scaffold coverage is complete. If the HAL
design conclusions were an artefact of its uneven coverage, they should not
reappear here.


### `judged_score`, and what it encodes

Harbor scores tasks with an automated verifier and has a judge LLM review
the verifier's output on Harbor Index. These judged scores are the adjudicated 
binary outcomes used in the analysis and represent instances where the Harbor
Index judge changed the verifier's classification (e.g., "True Solve").

### Inclusion threshold

Task and benchmark effects cannot be separated for a benchmark represented by
one or two items. The analysis uses benchmarks with **at least three tasks**,
which leaves 13 of 29 benchmarks and 1,098 of 1,476 rows.
`scripts/07_harbor_corroboration.R` re-runs the decomposition at thresholds of
three, four, and five so the sensitivity is visible rather than asserted; the
estimated model share falls as the threshold rises, which the paper attributes
to the non-random way tasks were selected to represent each benchmark.

---

## `hal_external_validation.csv`

Model-level scores from **four external benchmarks**, used to test whether the
latent model effect estimated on HAL transports outside the panel it was
fitted on. 101 rows, 28 models.

| Column | Description |
|---|---|
| `leaderboard`, `leaderboard_url` | Source and its URL |
| `benchmark` | `bfcl`, `terminal_bench_2_0`, `swebench`, `tau_2_core` |
| `benchmark_task_type` | What the benchmark measures |
| `evaluator`, `score_provenance` | Who ran it, and whether independently |
| `agent_scaffold_harness` | Harness the external score was obtained under |
| `task_set_size`, `metric` | Reported denominator and metric |
| `model_label_on_source` | Model name as printed by the source |
| `model_name` | Normalised to match the HAL response matrix |
| `score` | Reported accuracy |
| `rank_on_leaderboard`, `rank_denominator` | Position as published |

These are **scraped third-party leaderboard results, not reruns.** Several
sources report the same benchmark, so scores are averaged per
(benchmark, model) before use; different sources use different harnesses and
task subsets, which is part of why agreement is imperfect and why the
comparison is framed as convergent evidence rather than validation against
ground truth.

### Contamination control

Two of the four overlap a HAL benchmark in content:

| External benchmark | Overlaps | Handling |
|---|---|---|
| SWE-bench Verified | SWE-bench Verified Mini (a subset) | HAL benchmark removed from both estimators |
| τ²-bench Core | τ-bench Airline (a subset) | HAL benchmark removed from both estimators |
| BFCL v4 | — | full panel used |
| Terminal-Bench 2.0 | — | full panel used |

For the two overlapping comparisons, the latent effect is taken from the
matching leave-one-benchmark-out refit and the mean-score baseline is
recomputed on the remaining eight benchmarks. Correlating against a benchmark
the fit had already seen would credit the estimate for information it was
handed. `R/external.R` enforces this: it stops with an explanation if the
required refit is missing rather than falling back to the full-data fit.

---

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
