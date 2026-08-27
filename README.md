# How Reliable Are Agent Leaderboards? A Variance-Decomposition Analysis

Code and data for the ICLR submission. The paper applies generalizability
theory to nine benchmarks on the Holistic Agent Leaderboard (HAL),
decomposing agent score variance across model, task, and scaffold, and asks how
much of a published ranking would survive a redraw of the evaluation
conditions.

Three results, and where to find each in the code:

| Result | Script |
|---|---|
| Scaffold is a non-negligible measurement facet: for a fixed task, scaffold variation rivals or exceeds model variation | [`scripts/03_leaderboard_reliability.R`](scripts/03_leaderboard_reliability.R) |
| Benchmarks do not carry equal signal: only two of nine reach Eρ² > 0.75 for model ranking even with unlimited same-construction tasks | [`scripts/02_benchmark_reliability.R`](scripts/02_benchmark_reliability.R) |
| Reliability is design-conditional: task-only scaling asymptotes near Eρ² ≈ 0.44, broadening across benchmarks raises it to ≈ 0.75 | [`scripts/03_leaderboard_reliability.R`](scripts/03_leaderboard_reliability.R) |

---

## Quick start

```bash
Rscript install.R              # packages (CmdStan only needed to re-fit)
Rscript tests/test_gtheory.R   # check the reliability engine against the paper
Rscript tests/test_ranks.R     # check the ranking and rank-agreement machinery
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

## The idea, in one page

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
  disco.R                nonparametric distance-components decomposition
  plots.R                shared theme and one builder per figure
scripts/
  00_design_summary.R           Table 1, coverage, connectivity  (seconds)
  01_fit_models.R               estimates both decompositions    (hours)
  02_benchmark_reliability.R    Figures 1-3, ceilings table
  03_leaderboard_reliability.R  Figures 4-5, model vs. scaffold
  04_rank_analysis.R            Figures 6-7, Table 2, indistinguishability
  05_ablations.R                leave-one-out sensitivity analyses
  06_method_comparison.R        LME / GLME / Bayes / DISCO contrast
tests/
  test_gtheory.R         checks the engine reproduces Eq. 3-8 longhand
  test_ranks.R           checks rank direction, ties, and rank agreement
data/                    the response matrix, documented in data/README.md
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

ep2(vc, design(n_tasks = 50, n_benchmarks = 9))   # Eq. 7  model ranking
ep2(vc_b, design(n_tasks = 50), object = "system")# Eq. 4  model-scaffold system
inter_scaffold_reliability(vc_b, n_tasks = 50)    # Eq. 5  rho_AA'
reliability_ceiling(vc, n_benchmarks = 9)         # Eq. 8  task-only limit
snr(vc, design(50, 9))                            # Eq. 1  S/N_delta
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
| Appendix G.3 estimator comparison | `Rscript scripts/06_method_comparison.R` |
| Appendix C.5 DISCO | `RUN_DISCO=TRUE Rscript scripts/06_method_comparison.R` |

Figures go to `outputs/figures/`; tables are written as both CSV and a
booktabs LaTeX fragment to `outputs/tables/`.

### Cost

Wall-clock on an Apple M1 Max, matching Appendix C.6:

| Step | Time |
|---|---|
| Design summary; scripts 02–04 from cached fits | minutes |
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

**Wide intervals are a finding.** Most benchmarks contain two or three
scaffolds, so they carry little direct information about the distribution of
plausible harnesses. The posterior intervals reflect that honestly, and
narrowing them would mean hiding it.

---

## License

MIT — see [LICENSE](LICENSE). Rollout outcomes derive from the public HAL
release; see [`data/README.md`](data/README.md) for provenance.
