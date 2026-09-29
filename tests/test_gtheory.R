# tests/test_gtheory.R -----------------------------------------------------
#
# Checks that the declarative G-theory engine in R/gtheory.R reproduces the
# reliability equations as they are written in the paper.
#
# The engine derives every coefficient from a component's facet indexing. That
# is a good deal more compact than one function per formula, but the compaction
# is only safe if it is pinned to the published algebra -- which is what this
# file does. Each test writes an equation out longhand and compares.
#
# Run with:  Rscript tests/test_gtheory.R
#
# No test framework required: failures raise, and the script exits non-zero.

source(here::here("R", "setup.R"))

failures <- 0L
check <- function(label, actual, expected, tol = 1e-10) {
  ok <- if (is.character(expected) || is.character(actual)) {
    identical(as.character(actual), as.character(expected))
  } else {
    isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = tol))
  }
  cat(if (ok) "  ok   " else "  FAIL ", label, "\n", sep = "")
  if (!ok) {
    cat("         actual:   ", format(actual), "\n",
        "         expected: ", format(expected), "\n", sep = "")
    failures <<- failures + 1L
  }
  invisible(ok)
}

# Arbitrary but distinct component variances, so that any mis-assignment of a
# component to numerator, denominator, or divisor changes the answer.
bench_v <- c(I = 0.90, M = 0.20, A = 0.15, IM = 0.35, IA = 0.25, MA = 0.30, e = pi^2 / 3)
pool_v  <- c(B = 1.10, I = 0.90, M = 0.20, A = 0.15, BM = 0.22, BA = 0.18,
             MA = 0.30, IM = 0.35, IA = 0.25, BMA = 0.12, IMA = 0, e = pi^2 / 3)

vcb <- new_vcomp(as.list(bench_v), "benchmark",   source = "test fixture")
vcp <- new_vcomp(as.list(pool_v),  "leaderboard", source = "test fixture")

ni <- 40; na <- 3; nb <- 6


cat("\nBenchmark-level decomposition (Eq. 2)\n")

# Eq. 3: Erho^2_M(b) = sM / (sM + sIM/ni + sMA/na + sIMA,e/(ni na))
with(as.list(bench_v), {
  expected <- M / (M + IM / ni + MA / na + e / (ni * na))
  check("Eq. 3  model-ranking reliability",
        ep2(vcb, design(n_tasks = ni, n_scaffolds = na)), expected)

  check("Eq. 1  S/N is Erho^2 / (1 - Erho^2)",
        snr(vcb, design(n_tasks = ni, n_scaffolds = na)),
        expected / (1 - expected))
})

# Eq. 4: Erho^2_MA(b) = (sM+sA+sMA) / (sM+sA+sMA + (sIM+sIA+sIMA,e)/ni)
with(as.list(bench_v), {
  num <- M + A + MA
  check("Eq. 4  model-scaffold system reliability",
        ep2(vcb, design(n_tasks = ni), object = "system"),
        num / (num + (IM + IA + e) / ni))
})

# Eq. 5: rho_AA'(ni) = (sM + sIM/ni) / (sM + sMA + sIM/ni + sIMA,e/ni)
with(as.list(bench_v), {
  check("Eq. 5  inter-scaffold reliability",
        inter_scaffold_reliability(vcb, n_tasks = ni),
        (M + IM / ni) / (M + MA + IM / ni + e / ni))
})

# Phi adds error terms that shift all models together: the task main effect,
# the scaffold main effect, and the task-scaffold interaction.
with(as.list(bench_v), {
  check("Phi  dependability includes object-free components",
        phi(vcb, design(n_tasks = ni, n_scaffolds = na)),
        M / (M + IM / ni + MA / na + e / (ni * na) +
               I / ni + A / na + IA / (ni * na)))
})


cat("\nLeaderboard-level decomposition (Eq. 6)\n")

# Eq. 7: Erho^2_M = sM / (sM + sBM/nb + sMA/na + sBMA/(nb na)
#                          + sIM/(nb ni) + sBIMA,e/(nb ni na))
with(as.list(pool_v), {
  check("Eq. 7  pooled model-ranking reliability",
        ep2(vcp, design(n_tasks = ni, n_benchmarks = nb, n_scaffolds = na)),
        M / (M + BM / nb + MA / na + BMA / (nb * na) +
               IM / (nb * ni) + e / (nb * ni * na)))
})

# Eq. 8: the task-only ceilings. All item-indexed error vanishes; benchmark-
# and scaffold-indexed error does not.
with(as.list(bench_v), {
  check("Eq. 8  benchmark-level ceiling as ni -> Inf",
        reliability_ceiling(vcb, n_scaffolds = na),
        M / (M + MA / na))
})
with(as.list(pool_v), {
  check("Eq. 8  leaderboard ceiling as ni -> Inf",
        reliability_ceiling(vcp, n_benchmarks = nb, n_scaffolds = na),
        M / (M + BM / nb + MA / na + BMA / (nb * na)))
})

# Corollary 3.2: with sBM > 0 or sMA > 0 the ceiling is bounded away from 1.
check("Corollary 3.2  ceiling is strictly below 1",
      reliability_ceiling(vcp, n_benchmarks = 1e6, n_scaffolds = 1e6) < 1, TRUE)


cat("\nStructural properties\n")

# Monotonicity in each design dimension (Proposition 3.1).
grid <- dstudy(vcp, n_tasks = c(1, 10, 100, 1000, Inf), n_benchmarks = 9)
check("Erho^2 is non-decreasing in the number of tasks",
      all(diff(grid$estimate[grid$statistic == "ep2"]) >= -1e-12), TRUE)

check("Erho^2 is non-decreasing in the number of benchmarks",
      all(diff(purrr::map_dbl(1:9, \(nb) ep2(vcp, design(50, nb, 1)))) >= -1e-12), TRUE)

# Treating the model-scaffold pair as the object moves model-scaffold
# compatibility from error into signal, so it cannot lower reliability.
check("system reliability >= model reliability at the same task count",
      ep2(vcb, design(n_tasks = ni), object = "system") >=
        ep2(vcb, design(n_tasks = ni)), TRUE)

# Every component must be classified exactly once, with a divisor that is the
# product of the sampled facets it is indexed by.
part <- gt_partition(vcp, design(ni, nb, na))
check("all leaderboard components are classified",
      nrow(part$terms), nrow(component_registry("leaderboard")))
check("model main effect is the universe score",
      part$terms$role[part$terms$component == "M"], "universe")
check("benchmark main effect carries no model signal",
      part$terms$role[part$terms$component == "B"], "absolute_only")
check("task-model interaction divisor is nb * ni",
      part$terms$divisor[part$terms$component == "IM"], nb * ni)
check("residual divisor is nb * ni * na",
      part$terms$divisor[part$terms$component == "e"], nb * ni * na)

# A component absent from a fit contributes nothing, so dropping the four-way
# interaction term must not change a coefficient computed without it.
vcp_no4 <- new_vcomp(as.list(pool_v[setdiff(names(pool_v), "IMA")]), "leaderboard")
check("absent components are treated as zero variance",
      ep2(vcp_no4, design(ni, nb, na)), ep2(vcp, design(ni, nb, na)))


cat("\nVariance shares\n")
shares <- variance_shares(vcp)
check("variance shares sum to one", sum(shares$estimate), 1)
check("shares are reported for every component",
      nrow(shares), nrow(component_registry("leaderboard")))


cat("\nModel vs. scaffold, as the paper reports it\n")

# The revised framing compares single components pairwise rather than summing
# model-side against scaffold-side terms; `component_contrast()` below checks
# the arithmetic, and these check the wrapper wires it up the right way round.
ct <- model_vs_scaffold_contrast(vcp, design(1, 1, 1), scope = "main")
with(as.list(pool_v), {
  total <- B + I + M + A + BM + BA + MA + IM + IA + BMA + IMA + e
  check("main-effect wrapper contrasts M against A",
        ct$estimate, (M - A) / total)
})
check("scaffold and model exceedance probabilities are complementary",
      ct$p_model_exceeds_scaffold + ct$p_scaffold_exceeds_model, 1)

ct_task <- model_vs_scaffold_contrast(vcp, design(10, 1, 1), scope = "task")
with(as.list(pool_v), {
  ni <- 10
  total <- B + I / ni + M + A + BM + BA + MA + IM / ni + IA / ni + BMA + IMA / ni + e / ni
  check("task-interaction wrapper contrasts IM against IA",
        ct_task$estimate, (IM / ni - IA / ni) / total)
})

# The wrapper must also work on a decomposition with no benchmark facet.
ct_bench <- model_vs_scaffold_contrast(vcb, design(10, 1, 1), scope = "task")
with(as.list(bench_v), {
  ni <- 10
  total <- I / ni + M + A + MA + IM / ni + IA / ni + e / ni
  check("contrast works on a per-benchmark decomposition",
        ct_bench$estimate, (IM / ni - IA / ni) / total)
})


cat("\nScaffold as the object of measurement\n")

# Ranking scaffolds is the mirror image of ranking models: the model facet
# becomes the error term.
# The mirror is exact: whichever facet is not the object supplies the
# compatibility error term, so n_models plays the role n_scaffolds plays when
# ranking models.
nm <- 7
with(as.list(bench_v), {
  check("scaffold-ranking reliability mirrors model-ranking",
        ep2(vcb, design(n_tasks = ni, n_models = nm), object = "scaffold"),
        A / (A + IA / ni + MA / nm + e / (ni * nm)))
})
with(as.list(bench_v), {
  check("the two objects are exact mirrors under swapped counts",
        ep2(vcb, design(n_tasks = ni, n_scaffolds = nm), object = "model"),
        M / (M + IM / ni + MA / nm + e / (ni * nm)))
})

part_a <- gt_partition(vcp, design(ni, nb, na), object = "scaffold")
check("scaffold main effect is the universe score when ranking scaffolds",
      part_a$terms$role[part_a$terms$component == "A"], "universe")
check("model main effect carries no scaffold signal",
      part_a$terms$role[part_a$terms$component == "M"], "absolute_only")
check("model-scaffold interaction is relative error for both objects",
      part_a$terms$role[part_a$terms$component == "MA"], "relative_error")

# With the scaffold as object, the scaffold count is no longer an error
# divisor -- averaging over more scaffolds cannot sharpen a scaffold ranking.
check("scaffold count does not affect scaffold-ranking reliability",
      ep2(vcb, design(n_tasks = ni, n_scaffolds = 1), object = "scaffold"),
      ep2(vcb, design(n_tasks = ni, n_scaffolds = 9), object = "scaffold"))
check("model count does not affect model-ranking reliability",
      ep2(vcb, design(n_tasks = ni, n_models = 1), object = "model"),
      ep2(vcb, design(n_tasks = ni, n_models = 9), object = "model"))


cat("\nTask reliability\n")

with(as.list(bench_v), {
  check("Erho^2_M(i) is the model share against task-model interaction",
        task_reliability(vcb, "model"), M / (M + IM))
  check("Erho^2_A(i) is the scaffold share against task-scaffold interaction",
        task_reliability(vcb, "scaffold"), A / (A + IA))
  check("averaging over tasks raises task reliability",
        task_reliability(vcb, "model", n_tasks = 10), M / (M + IM / 10))
})


cat("\nComponent contrasts\n")

# The two contrasts the paper reports are pairwise on single components.
main <- component_contrast(vcp, "M", "A", design(1, 1, 1))
with(as.list(pool_v), {
  total <- B + I + M + A + BM + BA + MA + IM + IA + BMA + IMA + e
  check("main-effect contrast is (M - A) over total variance",
        main$estimate, (M - A) / total)
})

task_c <- component_contrast(vcp, "IM", "IA", design(1, 1, 1))
with(as.list(pool_v), {
  total <- B + I + M + A + BM + BA + MA + IM + IA + BMA + IMA + e
  check("task-interaction contrast is (IM - IA) over total variance",
        task_c$estimate, (IM - IA) / total)
})

check("a contrast of a component with itself is exactly zero",
      component_contrast(vcp, "M", "M", design(1, 1, 1))$estimate, 0)

check("reversing the sides negates the contrast",
      component_contrast(vcp, "A", "M", design(1, 1, 1))$estimate, -main$estimate)

# Components absent from a decomposition must contribute nothing rather than
# collapsing the arithmetic.
check("contrasts work on a decomposition lacking the named components",
      is.finite(component_contrast(vcb, c("M", "BM"), c("A", "BA"),
                                   design(1, 1, 1))$estimate), TRUE)

both <- model_vs_scaffold_contrast(vcp, design(1, 1, 1), scope = "both")
check("the paper reports two contrasts, not one", nrow(both), 2)
check("both contrasts are labelled", all(!is.na(both$label)), TRUE)


cat("\nDesign cost\n")

check("trials are tasks x benchmarks x scaffolds",
      evaluation_trials(10, 9, 2), 180)
check("cost is linear in trials and models",
      design_cost(10, 9, 2, n_models = 5, cost_per_trial = 2), 180 * 5 * 2)
check("adding scaffolds multiplies cost without adding tasks",
      design_cost(10, 9, 2, n_models = 1, cost_per_trial = 1) /
        design_cost(10, 9, 1, n_models = 1, cost_per_trial = 1), 2)


cat("\n", if (failures == 0L) "All checks passed.\n" else
    sprintf("%d check(s) FAILED.\n", failures), sep = "")
quit(status = if (failures == 0L) 0L else 1L)
