# Agent Evaluation Reliability: More Tasks Won't (Always) Fix an Agent Leaderboard

Code and data for the ICLR submission. The paper applies generalizability
theory to 22 agent benchmarks — nine from the Holistic Agent Leaderboard (HAL)
and 13 from the Harbor Index — decomposing score variance across model, task,
and scaffold, and asking **which conclusions an evaluation actually supports,
and what additional evaluation would improve them.**

The organising idea: reliability is not a property of a benchmark. It is a
property of a benchmark *and* a claim. The same data can rank deployable
systems precisely while leaving the underlying models unresolved.

Four results, and where to find each in the code:

| Result | Script |
|---|---|
| **Reliability depends on the measurement goal.** Fixed model–scaffold systems rank reliably (Eρ² ≈ 0.94–0.99); the underlying models do not (Eρ² ≈ 0.15–0.84) | [`02_benchmark_reliability.R`](scripts/02_benchmark_reliability.R) |
| **Changing the scaffold can change conclusions.** Inter-scaffold reliability spans 0.21–0.83; the scaffold changes *which tasks* get solved more than the model does | [`02`](scripts/02_benchmark_reliability.R), [`03`](scripts/03_leaderboard_reliability.R) |
| **More tasks cannot resolve all uncertainty.** Under observed scaffold coverage, unlimited tasks add at most ≈ 0.10 to model-ranking reliability | [`02`](scripts/02_benchmark_reliability.R), [`03`](scripts/03_leaderboard_reliability.R) |
| **Pooling diverse benchmarks helps, at lower cost.** Breadth lifts projected reliability from ≈ 0.44 to ≈ 0.75 at a fixed task budget; the latent model effect also transports better to held-out benchmarks | [`07`](scripts/07_harbor_corroboration.R), [`08`](scripts/08_external_validation.R), [`09`](scripts/09_cost_and_allocation.R) |

---

## Quick start

```bash
Rscript install.R              # packages (CmdStan only needed to re-fit)
Rscript tests/test_gtheory.R   # check the reliability engine against the paper
Rscript tests/test_ranks.R     # check the ranking and rank-agreement machinery
Rscript tests/test_design.R    # check the datasets, cost model, and subsampling
Rscript scripts/00_design_summary.R
```

The last command takes seconds, fits nothing, and reproduces Table 1 plus the
coverage statistics. If it prints 29,923 rollouts across 9 benchmarks, 54
models, and 13 scaffolds, the environment is working.

To reproduce the figures you need the model fits. Either drop previously saved
fits into `outputs/fits/`, or re-estimate:

```bash
REFIT=TRUE Rscript scripts/01_fit_models.R   # hours; see Cost below
make analysis                                # then the figures, in minutes
```

---

## The idea

A leaderboard reports one number per system, and readers treat the induced
ordering as a fact about models. Whether it is depends on how much of the
observed variation would persist under an admissible change to the evaluation —
different tasks of the same kind, a different scaffold, a different benchmark
in the battery.

Generalizability theory answers that by splitting score variance across the
facets of the design and asking what fraction supports the ordering of a chosen
**object of measurement**. Because the outcomes are binary and the design is
sparse and unbalanced, the split is estimated with Bayesian Bernoulli-logit
mixed models — random-item, many-facet Rasch models, where model capability,
task difficulty, and scaffold effects live on one latent log-odds scale.

Three things follow from choosing an object, and the code makes each of them
one argument rather than one function:

| Object | What counts as signal | What counts as error |
|---|---|---|
| `"model"` | persistent model differences | scaffold, task, benchmark, and their interactions with the model |
| `"scaffold"` | persistent scaffold differences | model, task, benchmark, and their interactions |
| `"system"` | the model–scaffold pair together, including their compatibility | task and benchmark variation only |

Two fits (`R/formulas.R`):

- **Benchmark level** (paper Eq. 2) — fitted within each benchmark:
  `score ~ 1 + (1|task_id) + (1|model_name) + (1|agent_name) + (1|task_id:model_name) + (1|task_id:agent_name) + (1|model_name:agent_name)`
- **Leaderboard level** (paper Eq. 6) — fitted jointly to all nine, with tasks
  nested in benchmarks and models and scaffolds crossed with benchmarks as far
  as the observed incidence graph supports.

Which object you choose changes the answer from the *same* posterior. When the
object is the **model**, scaffold effects are error, and no number of tasks
removes them. When the object is the **model–scaffold system**, the same
variation is signal. This is the paper's central inferential point, and it is
also the design principle of the code.

Two further datasets carry the argument beyond HAL:

- **Harbor Index** (`R/harbor.R`) — a second meta-benchmark with a different
  sparsity pattern: far fewer tasks, but the same models and scaffolds on every
  benchmark, where HAL's coverage varies benchmark to benchmark. It uses the
  *same* decomposition, so agreement is a genuine check rather than a
  restatement. (Harbor is not fully crossed on model × scaffold — 18 of 36
  pairs — so it corroborates the benchmark-breadth result, not
  scaffold-specific ones.)
- **External validation** (`R/external.R`) — four contemporaneous benchmarks
  outside the panel, used to test whether the estimated latent model effect
  transports, with explicit contamination control where content overlaps.

---

## Repository layout

```
.here                    project-root marker for here::here(); keep it
config.R                 paths, sampler settings, D-study grids — change things here
R/
  setup.R                sources everything below, in dependency order
  data.R                 loading, validation, design and coverage descriptives
  formulas.R             the random-effects structures of Eq. 2 and Eq. 6
  fit.R                  fit-or-load wrappers for brms and lme4, diagnostics
  variance.R             fits -> variance components (Bayesian or frequentist)
  gtheory.R              G-coefficients, signal-to-noise, D-studies
  ranks.R                latent capability, posterior ranks, rank agreement
  cost.R                 design cost, allocation, repeated task subsampling
  harbor.R               the Harbor Index corroboration dataset
  external.R             out-of-panel validation of the latent model effect
  disco.R                nonparametric distance-components decomposition
  plots.R                shared theme and one builder per figure
scripts/
  00_design_summary.R           design summary, coverage, connectivity (seconds)
  01_fit_models.R               estimates every decomposition        (hours)
  02_benchmark_reliability.R    per-benchmark reliability by object
  03_leaderboard_reliability.R  pooled D-study, variance decomposition
  04_rank_analysis.R            rank shifts, agreement, indistinguishability
  05_ablations.R                leave-one-out sensitivity analyses
  06_method_comparison.R        LME / GLME / Bayes / DISCO contrast
  07_harbor_corroboration.R     the same design question on a second dataset
  08_external_validation.R      does the latent model effect transport?
  09_cost_and_allocation.R      cost frontier and task subsampling
tests/
  test_gtheory.R         checks the engine reproduces the published equations
  test_ranks.R           checks rank direction, ties, and rank agreement
  test_design.R          checks the datasets, cost model, and subsampling
data/                    three datasets, documented in data/README.md
outputs/                 figures, tables, results, cached fits (all gitignored)
```

Analysis scripts are pure consumers: they load cached fits, compute, and draw.
Only `01_fit_models.R` estimates anything, and it refuses to run without
`REFIT=TRUE`, so re-running a figure can never silently reproduce it from a
different fit than the one it reports.

---

## The reliability engine

Nearly all of the paper's numbers are ratios of variance components, and the
usual way to write them — one hand-coded function per coefficient — is how such
code drifts out of agreement with its paper. Here there is one engine
(`R/gtheory.R`), driven by three declarations: the object of measurement, which
facets are sampled versus fixed, and how many levels of each the design has.
Everything else follows from how each component is indexed:

| | classification |
|---|---|
| universe score σ²_τ | components indexed only by object facets (and fixed facets) |
| relative error σ²_δ | components indexed by an object facet **and** a sampled facet — these can reorder objects |
| absolute error σ²_Δ | σ²_δ plus components carrying no object facet — these shift scores but not ranks |

Each error term is divided by the number of levels of the sampled facets it is
indexed by, because the design averages over them. Every published coefficient
is then one line:

```r
vc <- vcomp_from_brms(fit, level = "leaderboard")

ep2(vc, design(n_tasks = 50, n_benchmarks = 9))     # model ranking, pooled
ep2(vc_b, design(n_tasks = 50), object = "system")  # model-scaffold systems
ep2(vc_b, design(n_tasks = 50), object = "scaffold")# ranking scaffolds
inter_scaffold_reliability(vc_b, n_tasks = 50)      # rho_AA'
reliability_ceiling(vc, n_benchmarks = 9)           # the task-only limit
snr(vc, design(50, 9))                              # signal-to-noise
phi(vc_gaussian, design(50, 9))                     # absolute, not relative
```

Which facets count as error is derived from the object, not hardcoded: ranking
models makes the scaffold an error facet, ranking scaffolds makes the model
one. The two are exact mirrors, and a test asserts it.

Design economics sit on top of the same components:

```r
frontier <- cost_reliability_frontier(vc)            # reliability vs dollars
cheapest_design_reaching(frontier, 0.75)             # cheapest way to a target
subsample_rankings(df, k_values = c(5, 15, 30))      # model-free cross-check
```

Inter-scaffold reliability is not a special case in the implementation: it is
Eρ² for the model with the *task* facet held fixed, which is exactly what "would
two scaffolds rank these models the same way, on the same tasks?" means.

`tests/test_gtheory.R` writes Eq. 3, 4, 5, 7, and 8 out longhand and checks the
engine against them, along with the monotonicity claims of Proposition 3.1 and
the bound of Corollary 3.2. A change to the engine that breaks a paper formula
fails there.

The same engine consumes frequentist variance components
(`vcomp_from_lme4()`), which is why the ablations and the estimator comparison
reuse it rather than duplicating the algebra.

---

## Reproducing the paper

| Output | Command |
|---|---|
| Table 1, coverage, connectivity | `Rscript scripts/00_design_summary.R` |
| Figure 1 (inter-scaffold reliability) | `Rscript scripts/02_benchmark_reliability.R` |
| Figures 2–3 (per-benchmark D-studies, ceilings) | same |
| Figures 4–5 (leaderboard D-study, variance decomposition) | `Rscript scripts/03_leaderboard_reliability.R` |
| Section 5.1 scaffold-vs-model probabilities | same |
| Figures 6–7, Table 2 (rank shifts and agreement) | `Rscript scripts/04_rank_analysis.R` |
| Appendix D.5–D.6 ablations | `Rscript scripts/05_ablations.R` |
| Appendix D.4 leave-one-benchmark-out | `RUN_LOBO_BAYES=TRUE REFIT=TRUE Rscript scripts/05_ablations.R` |
| Estimator comparison (LME/GLME/Bayes) | `Rscript scripts/06_method_comparison.R` |
| DISCO nonparametric decomposition | `RUN_DISCO=TRUE Rscript scripts/06_method_comparison.R` |
| Harbor Index corroboration, HAL vs Harbor D-study | `Rscript scripts/07_harbor_corroboration.R` |
| Harbor task-threshold sensitivity | `RUN_HARBOR_THRESHOLDS=TRUE REFIT=TRUE Rscript scripts/07_harbor_corroboration.R` |
| Absolute reliability Φ, observed scale | `RUN_ABSOLUTE=TRUE REFIT=TRUE Rscript scripts/02_benchmark_reliability.R` |
| External validation, Kendall's τ table and Δτ | `Rscript scripts/08_external_validation.R` |
| Cost frontier, allocation, task subsampling | `Rscript scripts/09_cost_and_allocation.R` |

Figures go to `outputs/figures/`; tables are written as both CSV and a
booktabs LaTeX fragment to `outputs/tables/`.

### Cost

Wall-clock on an Apple M1 Max, matching Appendix C.6:

| Step | Time |
|---|---|
| Design summary; scripts 02–04, 07–08 from cached fits | minutes |
| Harbor Index fit (1,098 rollouts) | minutes |
| Cost frontier and task subsampling (script 09) | ~5 min |
| Estimator comparison (script 06): 20 lme4 refits | ~30–60 min |
| Pooled leaderboard fit (6 chains × 9,000 iter) | ~5.5 h |
| Nine benchmark-level fits | ~6.5 min each |
| Leave-one-benchmark-out Bayesian refits | ~57 h |
| Leave-one-model-out / leave-one-benchmark-scaffold-out (lme4) | tens of minutes |
| DISCO, nine benchmarks | ~17 h |
| PSIS-LOO with moment matching | ~38 h |

The expensive steps are all opt-in through environment flags. Sampler settings
live in `SAMPLER` in [`config.R`](config.R).

### Environment flags

| Flag | Effect |
|---|---|
| `REFIT=TRUE` | Re-estimate rather than load cached fits |
| `FIT_SENSITIVITY=TRUE` | Also fit the Gaussian-link and four-way-interaction variants |
| `RUN_LOO=TRUE` | PSIS-LOO diagnostics (Appendix D.2) |
| `RUN_LOBO_BAYES=TRUE` | Bayesian leave-one-benchmark-out refits |
| `RUN_DISCO=TRUE` | Nonparametric decomposition |
| `RUN_HARBOR_THRESHOLDS=TRUE` | Refit Harbor at task thresholds of 3, 4, 5 |
| `RUN_ABSOLUTE=TRUE` | Absolute reliability (Φ) from Gaussian per-benchmark fits |
| `N_SUBSAMPLES=n` | Replicates per k in the task-subsampling check (default 500) |
| `N_CORES=n` | Parallel workers for the frequentist ablations |

---

## Notes on reading the output

**Reliability is not a property of a benchmark.** It is a property of a
benchmark, a scaffold sample, a task sample, and the competitor set being
ranked. A clustered frontier-model pool depresses σ²_M regardless of benchmark
quality. The estimates here characterise the studied HAL designs and model
pool, and should be re-estimated as the evaluated population changes — not read
as fixed properties of benchmark names.

**Latent scale.** Components are estimated on the logit scale while
leaderboards report proportions; the two do not translate one-to-one, so the
coefficients are approximate bounds on observed-scale rank stability.

**Medians, not means.** Ratios of variance components are skewed, sometimes
severely. Figures report posterior medians with 0.68 HDI intervals
(`CI_LEVEL` in `config.R`), and `summarize_quantity()` returns both so either
can be checked.

**Variance shares.** `variance_shares()` returns `median` and
`median_normalized`. Quote `median` for a single component; `median_normalized`
exists only so stacked bars sum to one, because the median of a sum is not the
sum of medians.

**Dollar figures are assumptions, not measurements.** The response matrices
record outcomes, not spend, so the per-rollout price in `COST` (config.R) is an
input. The default is a single flat rate back-calculated so the full observed
battery reproduces its reported total, which makes every projection a
transparent linear function of one number. Real cost varies substantially by
benchmark, model, and episode length, so **the ordering of designs on the cost
frontier is far more trustworthy than the levels.** If you have per-benchmark
prices, set `COST$per_trial_by_benchmark` and every cost figure updates
consistently. Rank-agreement columns in the subsampling table do not depend on
the price and reproduce regardless.

**External validation is convergent evidence, not ground truth.** The external
scores are scraped third-party leaderboard results obtained under their own
harnesses, not reruns under controlled conditions. Agreement shows the latent
estimate carries information that transports; it does not establish that
either ranking is correct, nor that a single latent dimension captures agentic
capability. With six to twenty overlapping models per benchmark, single models
can move a rank correlation noticeably, which is why
`scripts/08_external_validation.R` reports leave-one-model-out ranges and flags
sign reversals alongside the headline numbers.

**Wide intervals are a finding.** Most benchmarks contain two or three
scaffolds, so they carry little direct information about the distribution of
plausible harnesses. The posterior intervals reflect that honestly, and
narrowing them would mean hiding it.

---

## License

MIT — see [LICENSE](LICENSE). Rollout outcomes derive from the public HAL
release; see [`data/README.md`](data/README.md) for provenance.
