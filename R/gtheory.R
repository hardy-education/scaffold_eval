# R/gtheory.R --------------------------------------------------------------
#
# Generalizability coefficients from variance components.
#
# This is the analytical core of the repository. Rather than hand-coding one
# function per coefficient -- which is how these formulas are usually written,
# and how they usually drift out of agreement with the paper -- everything is
# derived from three declarations:
#
#   1. the object of measurement (what is supposed to stay put when evaluation
#      conditions change),
#   2. which facets are sampled (random) and which are held fixed,
#   3. how many levels of each sampled facet the design contains.
#
# Given those, the classification of every variance component follows
# mechanically:
#
#   universe score   sigma^2_tau : components indexed only by object facets
#                                  (and fixed facets)
#   relative error   sigma^2_d   : components indexed by an object facet AND a
#                                  sampled facet -- these can reorder objects
#   absolute error   sigma^2_D   : sigma^2_d plus components carrying no object
#                                  facet -- these shift scores but not ranks
#
# and each error term is divided by the number of levels of the sampled facets
# it is indexed by, because the design averages over them.
#
# Every published coefficient is then a two-line specification. Concretely:
#
#   Eq. 3  Erho^2_M(b)   = ep2(vc_bench, design(n_tasks = ni, n_scaffolds = na))
#   Eq. 4  Erho^2_MA(b)  = ep2(vc_bench, design(n_tasks = ni), object = "system")
#   Eq. 5  rho_AA'(b)    = inter_scaffold_reliability(vc_bench, ni)
#   Eq. 7  Erho^2_M      = ep2(vc_pooled, design(ni, nb, na))
#   Eq. 8  the ni -> Inf limits are `n_tasks = Inf`
#
# `tests/test_gtheory.R` checks the engine against those equations written out
# by hand, so a change here that breaks a paper formula fails loudly.

OBJECT_FACETS <- list(
  model  = "M",        # a base LLM x reasoning-effort configuration
  system = c("M", "A") # the deployable model-scaffold pair
)

SAMPLED_FACETS <- c("B", "I", "A")


#' Describe a D-study design.
#'
#' @param n_tasks Tasks sampled per benchmark. `Inf` gives the task-only
#'   asymptote of Proposition 3.1.
#' @param n_benchmarks Benchmarks in the battery (leaderboard level only).
#' @param n_scaffolds Scaffolds each system is evaluated under. The observed
#'   leaderboard reports one scaffold per reported score, so `1` is the design
#'   current practice corresponds to.
design <- function(n_tasks = 1, n_benchmarks = 1, n_scaffolds = 1) {
  stopifnot(n_tasks > 0, n_benchmarks > 0, n_scaffolds > 0)
  c(I = n_tasks, B = n_benchmarks, A = n_scaffolds)
}


#' Split a facet string such as "BIM" into its individual facets.
facet_set <- function(s) strsplit(s, "", fixed = TRUE)[[1]]


#' Averaging divisor for each component.
#'
#' The number of levels of the sampled facets a component is indexed by. Tasks
#' are nested within benchmarks, so an item-indexed component at the
#' leaderboard level picks up both n_B and n_I.
#'
#' @param vc A `vcomp` object.
#' @param d A design from `design()`.
#' @param averaged Facets the design averages over. Defaults to every sampled
#'   facet; `gt_partition()` narrows it when a facet belongs to the object.
component_divisors <- function(vc, d, averaged = SAMPLED_FACETS) {
  divisors <- purrr::map_dbl(vc$registry$facets,
                             \(f) prod(d[intersect(averaged, facet_set(f))]))
  stats::setNames(divisors, vc$registry$component)
}

#' Design-scaled variance components: each divided by its averaging divisor.
#'
#' @return Named list of scaled components, in the registry's order.
scaled_components <- function(vc, d, averaged = SAMPLED_FACETS) {
  divisors <- component_divisors(vc, d, averaged)
  purrr::imap(vc$values, \(x, k) x / divisors[[k]])
}


#' Classify components and compute their design-scaled contributions.
#'
#' The workhorse behind every coefficient. Returns the three sums that make up
#' a G-coefficient, so callers only choose which to divide by which.
#'
#' @param vc A `vcomp` object.
#' @param d A design from `design()`.
#' @param object "model" or "system".
#' @param fixed Facets held fixed rather than resampled. Setting `fixed = "I"`
#'   asks how well two scaffolds agree *on the same tasks*, which is the
#'   inter-scaffold reliability of Eq. 5.
#' @return A list with `universe`, `relative_error`, `absolute_error`, and the
#'   per-component `terms` used to build them.
gt_partition <- function(vc, d, object = c("model", "system"), fixed = character()) {
  object <- match.arg(object)
  obj_facets <- OBJECT_FACETS[[object]]
  stopifnot(all(fixed %in% SAMPLED_FACETS))

  # A facet that is part of the object is not a source of error, and a facet
  # declared fixed is averaged over rather than generalized to.
  random <- setdiff(SAMPLED_FACETS, c(obj_facets, fixed))
  averaged <- setdiff(SAMPLED_FACETS, obj_facets) # random + fixed

  divisors <- component_divisors(vc, d, averaged)

  terms <- vc$registry |>
    mutate(
      has_object = purrr::map_lgl(facets, \(f) any(obj_facets %in% facet_set(f))),
      has_random = purrr::map_lgl(facets, \(f) any(random %in% facet_set(f))),
      divisor    = unname(divisors[component]),
      role = case_when(
        has_object & !has_random ~ "universe",
        has_object &  has_random ~ "relative_error",
        TRUE                     ~ "absolute_only"
      ),
      contribution = purrr::map2(component, divisor, \(k, div) vc$values[[k]] / div)
    )

  sum_role <- function(roles) {
    keep <- terms$contribution[terms$role %in% roles]
    if (!length(keep)) return(0)
    Reduce(`+`, keep)
  }

  list(
    universe       = sum_role("universe"),
    relative_error = sum_role("relative_error"),
    absolute_error = sum_role(c("relative_error", "absolute_only")),
    terms          = select(terms, component, label, facets, role, divisor)
  )
}


#' Relative generalizability coefficient, Erho^2 (paper Eq. 1).
#'
#' The squared correlation between an object's score under design `d` and its
#' universe score: how informative the observed ordering is about the ordering
#' a fresh draw from the universe would produce.
ep2 <- function(vc, d, object = "model", fixed = character()) {
  p <- gt_partition(vc, d, object, fixed)
  p$universe / (p$universe + p$relative_error)
}

#' Signal-to-noise ratio for relative decisions, S/N_delta (paper Eq. 1).
#'
#' Equal to Erho^2 / (1 - Erho^2). Reported alongside Erho^2 because the
#' limit-of-detection literature states minimum-detection criteria on this
#' scale.
snr <- function(vc, d, object = "model", fixed = character()) {
  p <- gt_partition(vc, d, object, fixed)
  p$universe / p$relative_error
}

#' Dependability coefficient Phi, for absolute rather than relative decisions.
#'
#' Includes error terms that shift all objects together -- a uniformly hard
#' task set, a systematically permissive scaffold. Those cancel from rankings
#' but not from "does this system exceed 60% accuracy" claims.
phi <- function(vc, d, object = "model", fixed = character()) {
  p <- gt_partition(vc, d, object, fixed)
  p$universe / (p$universe + p$absolute_error)
}


#' Inter-scaffold reliability, rho_AA' (paper Eq. 5).
#'
#' Would two independently sampled scaffolds rank the models the same way,
#' given the same `n_tasks` tasks? This is Erho^2 for the model with the task
#' facet held fixed: task-model interactions become part of the signal shared
#' by the two scaffolds, and only scaffold-indexed terms remain as error.
#'
#' Low values mean scaffold choice changes which models look strongest, even
#' when each scaffold on its own produces internally consistent scores.
inter_scaffold_reliability <- function(vc, n_tasks = 1, n_scaffolds = 1, n_benchmarks = 1) {
  ep2(vc,
      design(n_tasks = n_tasks, n_benchmarks = n_benchmarks, n_scaffolds = n_scaffolds),
      object = "model", fixed = "I")
}


#' Task-only reliability ceiling (paper Proposition 3.1 / Eq. 8).
#'
#' The limit of Erho^2 as tasks of the same construction accumulate without
#' bound. Benchmark- and scaffold-indexed error terms do not shrink with tasks,
#' so this ceiling is strictly below 1 whenever sigma^2_BM or sigma^2_MA > 0
#' (Corollary 3.2).
reliability_ceiling <- function(vc, n_benchmarks = 1, n_scaffolds = 1, object = "model") {
  ep2(vc, design(Inf, n_benchmarks, n_scaffolds), object = object)
}


# ---- D-studies ------------------------------------------------------------

#' Sweep reliability over a grid of evaluation designs.
#'
#' Produces the tidy table behind the paper's D-study figures. Each row is one
#' hypothetical design; `total_tasks` is the evaluation budget actually spent
#' (benchmarks x tasks per benchmark), which is what makes the "500 tasks from
#' one benchmark vs. 55 from each of nine" comparison a fair one.
#'
#' @param vc A `vcomp`.
#' @param n_tasks Vector of task counts per benchmark.
#' @param n_benchmarks Vector of benchmark counts.
#' @param n_scaffolds Vector of scaffold counts.
#' @param objects Which objects of measurement to evaluate.
#' @param statistics Which coefficients to compute.
#' @param ci Credible mass for interval summaries.
#' @param max_total_tasks Optional cap on `n_benchmarks * n_tasks`, to keep the
#'   sweep to designs within a plausible budget.
dstudy <- function(vc,
                   n_tasks,
                   n_benchmarks = 1,
                   n_scaffolds = DSTUDY$n_scaffolds,
                   objects = "model",
                   statistics = c("ep2", "snr"),
                   ci = CI_LEVEL,
                   max_total_tasks = NULL) {
  stat_fun <- list(ep2 = ep2, snr = snr, phi = phi)
  statistics <- match.arg(statistics, names(stat_fun), several.ok = TRUE)

  grid <- tidyr::expand_grid(
    n_benchmarks = n_benchmarks,
    n_tasks      = n_tasks,
    n_scaffolds  = n_scaffolds,
    object       = objects,
    statistic    = statistics
  ) |>
    mutate(total_tasks = n_benchmarks * n_tasks)

  if (!is.null(max_total_tasks)) {
    grid <- filter(grid, is.infinite(total_tasks) | total_tasks <= max_total_tasks)
  }

  grid |>
    mutate(
      .summary = purrr::pmap(
        list(n_tasks, n_benchmarks, n_scaffolds, object, statistic),
        function(ni, nb, na, obj, st) {
          summarize_quantity(
            stat_fun[[st]](vc, design(ni, nb, na), object = obj),
            ci = ci
          )
        }
      )
    ) |>
    tidyr::unnest(.summary)
}


#' D-study across every benchmark's own fit.
#'
#' @param vcs Named list of `vcomp` objects, one per benchmark.
#' @param n_tasks Task grid, or a named list of per-benchmark grids.
dstudy_by_benchmark <- function(vcs, n_tasks, ...) {
  purrr::imap(vcs, function(vc, b) {
    grid <- if (is.list(n_tasks)) n_tasks[[b]] else n_tasks
    dstudy(vc, n_tasks = grid, ...) |> mutate(benchmark = b, .before = 1)
  }) |>
    bind_rows()
}


# ---- Model vs. scaffold variance -------------------------------------------

#' Posterior contrast between model-related and scaffold-related variation.
#'
#' Section 5.1 asks whether the scaffold contributes as much rank-relevant
#' variation as the model does. The contrast is taken as a share of total
#' design-scaled variance, so it is on a common footing across designs, and its
#' `p_scaffold_exceeds_model` is the probability quoted in the text.
#'
#' @param scope "benchmark" compares the variation seen when ranking on a
#'   randomly chosen benchmark (model and benchmark-model terms against
#'   scaffold and benchmark-scaffold terms). "task" additionally includes the
#'   task-indexed interactions, i.e. the comparison at the level of a single
#'   task within a benchmark.
#' @return A one-row tibble of posterior summaries, plus the draws themselves
#'   in the `draws` list-column for plotting.
model_vs_scaffold_contrast <- function(vc,
                                       d = design(1, 1, 1),
                                       scope = c("benchmark", "task"),
                                       ci = CI_LEVEL) {
  scope <- match.arg(scope)

  # Design-scaled components, so that the contrast is on the same footing as
  # the reliability coefficients. Components a decomposition does not contain
  # (there is no benchmark facet in a per-benchmark fit) contribute nothing.
  sc <- scaled_components(vc, d)
  term <- function(k) if (is.null(sc[[k]])) 0 else sc[[k]]

  model_side    <- term("M") + term("BM")
  scaffold_side <- term("A") + term("BA")

  if (scope == "task") {
    model_side    <- model_side + term("IM")
    scaffold_side <- scaffold_side + term("IA")
  }

  total <- Reduce(`+`, sc)
  contrast <- (model_side - scaffold_side) / total

  summarize_quantity(contrast, ci = ci) |>
    mutate(
      scope = scope,
      n_tasks = d[["I"]], n_benchmarks = d[["B"]], n_scaffolds = d[["A"]],
      p_model_exceeds_scaffold = prob_positive(contrast),
      p_scaffold_exceeds_model = 1 - prob_positive(contrast),
      draws = list(if (inherits(contrast, "rvar")) as.numeric(posterior::draws_of(contrast)) else contrast),
      .before = 1
    )
}
