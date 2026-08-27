# tests/test_ranks.R -------------------------------------------------------
#
# Checks the rank machinery on synthetic posterior draws with a known answer.
#
# Two things here are easy to get subtly wrong and hard to notice in a figure:
# the direction of a percentile rank (Figure 6 plots one against the other with
# a diagonal, so a reversed axis looks like a dramatic finding), and the
# mapping of brms random-effect level labels back to models when model names
# themselves contain underscores.
#
# Run with:  Rscript tests/test_ranks.R

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

# Three models with unambiguous ordering: c is best, a is worst, on every draw.
draws <- tidyr::expand_grid(
  .draw = 1:200,
  model_name = c("a_default", "b_high", "c_low")
) |>
  mutate(
    benchmark = factor("gaia", levels = BENCHMARKS),
    theta = c(a_default = -1, b_high = 0, c_low = 1)[model_name] + 1e-6 * .draw
  )

cat("\nRank direction\n")

int_ranks <- posterior_ranks(draws)
check("integer rank 1 is the best model",
      int_ranks$model_name[which.min(int_ranks$rank_median)], "c_low")
check("integer rank is largest for the worst model",
      int_ranks$model_name[which.max(int_ranks$rank_median)], "a_default")

pct_ranks <- posterior_ranks(draws, percentile = TRUE)
check("percentile rank 1 is the best model",
      pct_ranks$model_name[which.max(pct_ranks$rank_median)], "c_low")
check("percentile ranks span [0, 1]",
      range(pct_ranks$rank_median), c(0, 1))

# The published-rank axis must run the same way as the percentile axis, or the
# diagonal in Figure 6 is meaningless.
obs <- tibble(
  run_id = as.character(1:3), benchmark = factor("gaia", levels = BENCHMARKS),
  task_id = "t1", agent_name = "s1",
  model_name = c("a_default", "b_high", "c_low"),
  reasoning_effort = "default", score = c(0, 0.5, 1)
)
pub <- published_scores(obs)
check("published percentile rank agrees in direction with posterior percentile",
      stats::cor(pub$published_pctrank[order(pub$model_name)],
                 pct_ranks$rank_median[order(pct_ranks$model_name)]), 1)
check("published integer rank 1 is the best model",
      pub$model_name[which.min(pub$published_rank)], "c_low")

cat("\nOrdering probabilities\n")

pw <- pairwise_ordering(draws, threshold = 0.95)
check("a strictly dominated model never wins its pairs",
      all(pw$p_a_beats_b %in% c(0, 1)), TRUE)
check("every pair is resolved when the ordering is deterministic",
      all(pw$resolved), TRUE)
check("share of indistinguishable pairs is zero here",
      indistinguishable_share(draws)$prop_indistinguishable, 0)

# With no signal at all, nothing should be resolvable.
set.seed(7)
noise <- draws |> mutate(theta = stats::rnorm(dplyr::n()))
check("pure noise resolves almost nothing",
      indistinguishable_share(noise)$prop_indistinguishable > 0.5, TRUE)

cat("\nRank agreement\n")

check("a perfectly matching ordering gives Spearman 1",
      rank_agreement(c(1, 2, 3), c(10, 20, 30))$spearman, 1)
check("a reversed ordering gives Spearman -1",
      rank_agreement(c(1, 2, 3), c(30, 20, 10))$spearman, -1)

cat("\nRank shifts\n")

shifts <- rank_shifts(int_ranks, pub)
check("no rank moves when the two orderings agree", shifts$delta, c(0, 0, 0))

cat("\n", if (failures == 0L) "All checks passed.\n" else
    sprintf("%d check(s) FAILED.\n", failures), sep = "")
quit(status = if (failures == 0L) 0L else 1L)
